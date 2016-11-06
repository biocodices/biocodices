#!/usr/bin/env python
"""
Start N TOR instances!

Usage:
    start_tor_instances.py <number-of-nodes> [--start-port 9052]
"""

import os
from sys import argv
from concurrent.futures import ThreadPoolExecutor
from collections import namedtuple

from docopt import docopt


TorNode = namedtuple('TorNode', 'socks_port control_port')


def define_tor_nodes(number_of_nodes, start_port):
    start_port = (start_port and int(start_port)) or 9050
    nodes = [TorNode(start_port, start_port + 1)]
    for _ in range(number_of_nodes - 1):
        last_socks_port, last_control_port = nodes[-1]
        nodes.append(TorNode(last_socks_port + 2, last_control_port + 2))
    return nodes

def start_tor_node(node):
    data_dir = os.path.join('/tmp', 'tor-{0}-{1}'.format(*node))
    command = ('tor --SocksPort {0} --ControlPort {1} '
               '--PidFile tor-{0}-{1}.pid --DataDirectory {2}')
    command = command.format(node.socks_port, node.control_port, data_dir)
    os.system(command)

def start_tor_nodes(number_of_nodes, start_port):
    nodes = define_tor_nodes(number_of_nodes, start_port=start_port)
    with ThreadPoolExecutor(max_workers=number_of_nodes) as thread_pool:
        thread_pool.map(start_tor_node, nodes)

if __name__ == '__main__':
    arguments = docopt(__doc__)
    start_tor_nodes(number_of_nodes=int(arguments['<number-of-nodes>']),
                    start_port=arguments['--start-port'])
