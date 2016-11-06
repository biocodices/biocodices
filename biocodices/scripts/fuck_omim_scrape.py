#!/usr/bin/env python

import time
from os import sys

import stem
from stem import Signal
from stem.control import Controller
import requests
import pandas as pd

from biocodices.annotation import Omim
from biocodices.helpers import in_groups_of

# 1. Start n tor instances with SocksPort ControlPort and PidFile
# 2. Keep the info about the ports + pidfile
# 3. Start n python processes, each will use one of the available Tor instances.
# 4. After annotating, kill the tor instances using the stored PID.






#  proxies = {'http': 'socks5://127.0.0.1:9050'}
#  proxies = {'http': 'socks5://127.0.0.1:9052'}
headers = {'User-Agent': ('Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) '
                          'AppleWebKit/537.36 (KHTML, like Gecko) '
                          'Chrome/39.0.2171.95 Safari/537.36')}



#  def start_tor_instance():
    #  return socks_port, control_port

#  def my_ip():
    #  r = requests.get('http://canihazip.com/s', proxies=proxies)
    #  return r.text.strip()

#  def change_identity():
    #  with Controller.from_port(port=9053) as controller:
        #  controller.authenticate(password='elhobbit')
        #  controller.signal(Signal.NEWNYM)
        #  controller.close()

#  def main():
    #  omim_genes = pd.read_csv('~/omim_genes.csv')
    #  omim_gene_ids = list(omim_genes.mim.dropna())

    #  omim_annotator = Omim()
    #  omim_annotator.PROXIES = proxies

    #  for omim_ids in in_groups_of(100, omim_gene_ids):
        #  print('My IP:', my_ip())
        #  html_dict = omim_annotator.annotate(omim_ids)
        #  change_identity()


if __name__ == '__main__':
    main()
