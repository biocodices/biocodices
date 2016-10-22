import re
import time
import requests
from multiprocessing import Pool

import pandas as pd
from bs4 import BeautifulSoup

from biocodices.helpers import randomize_sleep_time
from biocodices.annotation import AnnotatorWithCache


class Omim(AnnotatorWithCache):
    @staticmethod
    def _key(mim_id):
        return 'omim:%s' % mim_id

    @staticmethod
    def _url(mim_id):
        return 'http://omim.org/entry/{0}'.format(mim_id)

    @classmethod
    def _query(cls, mim_id):
        user_agent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) ' + \
                     'AppleWebKit/537.36 (KHTML, like Gecko) ' + \
                     'Chrome/39.0.2171.95 Safari/537.36'

        response = requests.get(cls._url(mim_id),
                                headers={'user-agent': user_agent})
        if not response.ok:
            print(response.status_code, 'status code for id "{0}"'.format(mim_id))

        return response.text

    def _batch_query(self, mim_ids, parallel, sleep_time):
        html_dict = {}

        # OMIM seems to be very strict against web crawlers and bans IPs.
        # That's why we override the _batch_query method and avoid
        # parallelization here.
        # Never ommit the sleeping step between queries, and never let it
        # be less than a couple of seconds, just in case:
        if sleep_time < 3:
            sleep_time = 3

        for i, mim_id in enumerate(mim_ids):
            if i > 0:
                # To further simulate real human behavior when visiting the page,
                # randomize the sleeping time:
                random_sleep_time = randomize_sleep_time(sleep_time)
                print('  Random sleep: {0:.2f} seconds'.format(random_sleep_time))
                time.sleep(random_sleep_time)

            print('[%s/%s] Visit %s' % (i+1, len(mim_ids), self._url(mim_id)))
            html = self._query(mim_id)
            self._cache_set({mim_id: html})
            html_dict[mim_id] = html

        return html_dict

    @classmethod
    def parse_html_dict(cls, html_dict):
        """
        Parses a dict like { '60555': '<html ...>' }. Meant for the result
        of #annotate() for this class.
        Returns a DataFrame of the OMIM variants per OMIM entry.
        """
        entries = cls._html_dict_to_entries_list(html_dict)
        references = cls._html_dict_to_references_list(html_dict)

        for entry in entries:
            if 'pubmeds_summary' not in entry:
                continue
            pmids = [pmid for pmid in entry['pubmeds_summary'].values() if pmid]
            entry['pubmeds'] = [ref for ref in references
                                if 'pmid' in ref and ref['pmid'] in pmids]

        df = pd.DataFrame(entries)

        df['sub_id'] = df['pheno'].str.extract(r'\.(\d+) ', expand=False)
        df['phenotype'] = df['pheno'].str.extract(r'\.\d+ (.*)', expand=False)
        df.drop('pheno', axis=1, inplace=True)

        df['gene'] = df['variant'].str.extract(r'^(\w+),', expand=False)
        df['rs'] = df['variant'].str.findall(r'dbSNP:(rs\d+)').str.join('|')

        prot_regex = r'\w+, (.+?)(?:,| \[| -| \()'
        df['prot_change'] = df['variant'].str.extract(prot_regex, expand=False)
        return df

    @classmethod
    def _html_dict_to_references_list(cls, html_dict):
        references = []

        with Pool(7) as pool:
            results = [pool.apply_async(cls._parse_html_references, (html, ))
                       for html in html_dict.values()]
            for result in results:
                references += (result.get() or [])

        return references

    @classmethod
    def _html_dict_to_entries_list(cls, html_dict):
        entries = []

        with Pool(7) as pool:
            results = [pool.apply_async(cls._parse_html, (html, omim_id))
                       for omim_id, html in html_dict.items()]
            for result in results:
                entries += (result.get() or [])

        return entries

    @classmethod
    def _parse_html_references(cls, html):
        soup = BeautifulSoup(html, 'html.parser')
        references = []

        in_the_references_section = False
        for td in soup.select('td#floatingEntryContainer .wrapper-table td'):
            if not in_the_references_section:
                # Check again if the section started:
                in_the_references_section = (td.get('id') == 'references')
                continue

            # Get in the references section
            for td2 in td.select('.wrapper-table td.text'):
                if td2.get('id') and 'reference' in td2.get('id'):
                    continue

                fields = ['authors', 'title', 'publication', 'pubmed_str']
                values = [span.text.strip() for span in td2.select('span')]
                reference = dict(zip(fields, values))

                if 'pubmed_str' in reference:
                    pubmed_match = re.search(r'PubMed: (\d+)',
                                             reference['pubmed_str'])
                    if pubmed_match:
                        reference['pmid'] = pubmed_match.group(1)
                    del(reference['pubmed_str'])

                references.append(reference)

            break  # Leaving the references section, ignore the rest of <td>s

        return references

    @classmethod
    def _parse_html(cls, html, omim_id):
        soup = BeautifulSoup(html, 'html.parser')
        entry_type = ' '.join([e['title'] for e in soup.select('.title.definition')])
        if 'Phenotype description' in entry_type:
            return None

        entries = []
        current_entry = {}

        in_the_variants_section = False
        table_rows = soup.select('td#floatingEntryContainer .wrapper-table td')
        for td in [td for td in table_rows if td.text]:
            inner_text = re.sub(r'\s+', ' ', td.text).strip()

            # Detect start and end of the ALLELIC VARIANTS section
            if not in_the_variants_section:
                in_the_variants_section = re.search(r'Table View', inner_text)
                continue
            elif re.search(r'REFERENCES', inner_text):
                entries.append(current_entry) if current_entry else None
                break

            # Title of new entry
            if re.match(r'\.\d{4} ', inner_text):
                entries.append(current_entry) if current_entry else None
                current_entry = {}
                if not re.search(r'REMOVED FROM|MOVED TO', inner_text):
                    current_entry['pheno'] = inner_text
                continue

            if 'variant' not in current_entry:
                if '[ClinVar]' in inner_text:
                    current_entry['variant'] = inner_text.replace(' [ClinVar]', '')
                else:
                    current_entry['pheno'] += ' {0}'.format(inner_text)
                continue

            # Rest of the entry will be the review
            # Extract the PubMed references
            if 'pubmed' not in current_entry:
                current_entry['pubmeds_summary'] = {}

            for anchor in td.select('a.entry-reference'):
                current_entry['pubmeds_summary'][anchor.text] = anchor.get('pmid')

            # Get the review text
            if 'review' not in current_entry:
                current_entry['review'] = td.text
                continue

            current_entry['review'] += td.text

        for entry in entries:
            entry['omim_id'] = omim_id

        return entries

