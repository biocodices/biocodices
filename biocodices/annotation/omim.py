import time
import requests

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
            response.raise_for_status()

        return response.text

    def _batch_query(self, mim_ids, parallel, sleep_time):
        html_dict = {}

        for i, mim_id in enumerate(mim_ids):
            print('Visit %s' % self._url(mim_id))
            html_dict[mim_id] = self._query(mim_id)

            if i > 0:
                # OMIM seems to be very strict against web crawlers and bans IPs
                # Never ommit this sleeping step between queries, and never let it
                # be less than a couple of seconds, just in case:
                if sleep_time < 5:
                    sleep_time = 5
                # To further simulate real human behavior when visiting the page,
                # randomize the sleeping time:
                sleep_time = randomize_sleep_time(sleep_time)
                print('Sleep random time: %s seconds' % sleep_time)
                time.sleep(sleep_time)

        self._cache_set(html_dict)
        return self._cache_get([mim_ids])

    @classmethod
    def parse_html_dict(cls, html_dict):
        """
        Parses a dict like { '60555': '<html ...>' }. Meant for the result
        of #annotate() for this class.
        Returns a dict like { '60555': { data about the OMIM ID } }
        """
        return {key: cls.parse_html(html) for key, html in html_dict.items()}

    @staticmethod
    def parse_html(html):
        soup = BeautifulSoup(html)

        info = {
        }
        return info
