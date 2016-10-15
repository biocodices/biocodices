import re
import time
import requests
import json
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
        # Never ommit the sleeping step between queries, and never let it
        # be less than a couple of seconds, just in case:
        if sleep_time < 3:
            sleep_time = 3

        for mim_id in mim_ids:
            print('Visit %s' % self._url(mim_id))
            html = self._query(mim_id)
            self._cache_set({mim_id: html})
            html_dict[mim_id] = html

            if len(mim_ids) > 1:
                # To further simulate real human behavior when visiting the page,
                # randomize the sleeping time:
                random_sleep_time = randomize_sleep_time(sleep_time)
                print('Random sleep: {0:.2f} seconds'.format(random_sleep_time))
                time.sleep(random_sleep_time)

        return html_dict

    @classmethod
    def _entries_from_htmls(cls, html_dict):
        entries = []

        with Pool(7) as pool:
            results = [pool.apply_async(cls._entries_from_omim_html, (html, omim_id))
                       for omim_id, html in html_dict.items()]
            for result in results:
                entries += (result.get() or [])

        return entries


    @classmethod
    def dataframe_from_htmls(cls, html_dict):
        """
        Parses a dict like { '60555': '<html ...>' }. Meant for the result
        of #annotate() for this class.
        Returns a DataFrame of the OMIM variants per OMIM entry.
        """
        entries = cls._entries_from_htmls(html_dict)
        for entry in entries:
            if 'pubmeds' in entry:
                entry['pubmeds'] = json.dumps(entry['pubmeds'])
        df = pd.DataFrame(entries)

        df['sub_id'] = df['pheno'].str.extract(r'\.(\d+) ', expand=False)
        df['phenotype'] = df['pheno'].str.extract(r'\.\d+ (.*)', expand=False)
        df.drop('pheno', axis=1, inplace=True)

        df['gene'] = df['variant'].str.extract(r'^(\w+),', expand=False)
        df['rs'] = df['variant'].str.extract(r'dbSNP:(rs\d+)', expand=False)
        prot_regex = r'\w+, (.+?)(?:,| \[| -| \()'
        df['prot_change'] = df['variant'].str.extract(prot_regex, expand=False)
        return df

    @classmethod
    def _entries_from_omim_html(cls, html, omim_id):
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
                current_entry['pubmeds'] = {}

            for anchor in td.select('a.entry-reference'):
                current_entry['pubmeds'][anchor.text] = anchor.get('pmid')

            # Get the review text
            if 'review' not in current_entry:
                current_entry['review'] = td.text
                continue

            current_entry['review'] += td.text

        for entry in entries:
            entry['omim_id'] = omim_id

        return entries

