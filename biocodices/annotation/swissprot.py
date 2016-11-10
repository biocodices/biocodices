import re
from collections import namedtuple
from functools import lru_cache
from urllib.error import HTTPError
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from Bio import SwissProt as SwissProtReader
from Bio import ExPASy
from Bio.SeqUtils import seq3
from tqdm import tqdm


class SwissProt():
    def __init__(self, swissprot_id):
        self.id = swissprot_id

    def __repr__(self):
        return '{}({!r})'.format(self.__class__.__name__, self.id)

    @classmethod
    def variants_from_gene(cls, swissprot_id):
        return cls(swissprot_id).variants

    @classmethod
    def variants_from_genes(cls, swissprot_ids):
        with ThreadPoolExecutor(max_workers=10) as pool:
            futures = [pool.submit(cls.variants_from_gene, swissprot_id)
                       for swissprot_id in swissprot_ids]

            futures_iter = tqdm(as_completed(futures),
                                total=len(swissprot_ids), position=0)

        frames = [f.result() for f in futures_iter]
        return pd.concat(frames, ignore_index=True)

    @staticmethod
    @lru_cache(maxsize=1000)
    def fetch_swissprot_record(swissprot_id):
        try:
            handle = ExPASy.get_sprot_raw(swissprot_id)
        except ValueError:
            print('  "%s" not found in SwissProt' % swissprot_id)
        else:
            record = SwissProtReader.read(handle)
            return record

    @property
    def record(self):
        return self.__class__.fetch_swissprot_record(self.id)

    @property
    def comments(self):
        comments = {}
        for comment in self.record.comments:
            title, text = comment.split(': ', maxsplit=1)
            comments[title.capitalize()] = text
        return comments

    @staticmethod
    def _extract_first_group(regex, text):
        match = re.search(regex, text)
        return match and match.group(1)

    @property
    def description(self):
        return self._extract_first_group(r'RecName: Full=(.+?);( AltName)?',
                                         self.record.description)

    @property
    def alt_name(self):
        return self._extract_first_group(r'AltName: Full=(.+?);',
                                         self.record.description)

    @property
    def gene_symbol(self):
        return self._extract_first_group(r'Name=(.+?)( \{.+\})?;( Synonyms)?',
                                         self.record.gene_name)

    @property
    @lru_cache()
    def variants(self):
        variant_dicts = [self._parse_variant(feature)
                         for feature in self.record.features
                         if feature[0] == 'VARIANT']
        return pd.DataFrame(variant_dicts)

    def _parse_variant(self, feature):
        variant = {
            'id': self.id,
            'url': 'http://www.uniprot.org/uniprot/{}'.format(self.id),
            'variant_id': feature[4],
            'desc': feature[3],
            'pos': feature[1],
            'pos_stop': feature[2],
            'gene_symbol': self.gene_symbol,
        }

        variant_url = 'http://web.expasy.org/variant_pages/{}.html'
        variant['variant_url'] = variant_url.format(variant['variant_id'])

        regex = r'(?P<old_aa>[A-Z]+) -> (?P<new_aa>[A-Z]+)'
        match = re.search(regex, variant['desc'])
        if match:
            variant['old_aa'], variant['new_aa'] = match.groups()
            variant['prot_change'] = 'p.{}{}{}'.format(seq3(variant['old_aa']),
                    variant['pos'], seq3(variant['new_aa']))

        variant['pmids'] = re.findall(r'PubMed:(\d+)', variant['desc'])

        matches = re.search(r'dbSNP:(rs\d+)', variant['desc'])
        if matches:
            assert len(matches.groups()) == 1
            variant['dbSNP_id'] = matches.group(1)

        match = re.search(r'\((.+)\)', variant['desc'])
        if match:
            review = re.sub(r'; dbSNP:rs\d+', '', match.group(1))
            variant['review'] = review

        return variant

