#!/usr/bin/env python

from beeprint import pp
from glob import glob
from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import ceil

import requests
import pandas as pd
#  from stem import Signal
#  from stem.control import Controller

from biocodices.annotation import Omim
from biocodices.helpers import in_groups_of


headers = {'User-Agent': ('Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) '
                          'AppleWebKit/537.36 (KHTML, like Gecko) '
                          'Chrome/39.0.2171.95 Safari/537.36')}

TorNode = namedtuple('TorNode', 'socks_port control_port pid')

def get_running_tor_nodes():
    tor_nodes = []
    for pidfile in glob('*.pid'):
        socks_port, control_port = pidfile.split('.')[0].split('-')[1:]
        with open(pidfile) as f:
            pid = f.read().strip()
        tor_nodes.append(TorNode(socks_port, control_port, pid))
    return tor_nodes

def tor_node_to_proxy(tor_node):
    return {'http': 'socks5://localhost:{}'.format(tor_node.socks_port)}

#  def change_identity(tor_node):
    #  with Controller.from_port(port=int(tor_node.control_port)) as controller:
        #  controller.authenticate()
        #  controller.signal(Signal.NEWNYM)
        #  controller.close()

def my_ip(proxy_dict):
    response = requests.get('http://canihazip.com/s', proxies=proxy_dict)
    return response.text.strip()

def get_mim_ids():
    omim_genes = pd.read_csv('~/omim_genes.csv')
    return list(omim_genes.mim.dropna())

def annotate_mim_ids_in_tor_node(mim_ids, tor_node):
    proxies = tor_node_to_proxy(tor_node)
    omim_annotator = Omim()
    omim_annotator.PROXIES = proxies
    node_ip = my_ip(proxies)
    omim_annotator.TQDM_PREFIX = 'TOR @ %s' % node_ip
    html_dict = omim_annotator.annotate(mim_ids)
    return html_dict


def main():
    tor_nodes = get_running_tor_nodes()
    mim_ids = get_mim_ids()
    group_size = ceil(len(mim_ids) / len(tor_nodes))
    grouped_mim_ids = in_groups_of(group_size, mim_ids)

    with ThreadPoolExecutor(max_workers=len(tor_nodes)) as executor:
        for tor_node, mim_ids_group in zip(tor_nodes, grouped_mim_ids):
            executor.submit(annotate_mim_ids_in_tor_node, mim_ids_group, tor_node)


if __name__ == '__main__':
    main()
