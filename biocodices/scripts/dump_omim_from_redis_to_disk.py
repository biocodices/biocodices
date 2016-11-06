#!/usr/bin/env python

from os import makedirs
from os.path import expanduser, join, isfile
import time
from concurrent.futures import ProcessPoolExecutor

import redis
from tqdm import tqdm as with_progressbar


redis_client = redis.StrictRedis(host='localhost', port='6379', db=0)
omim_keys = [key.decode('utf-8') for key in redis_client.keys('omim*')]

dest_dir = expanduser('~/new_500gb/omim_dump')
makedirs(dest_dir, exist_ok=True)

for key in with_progressbar(omim_keys):
    html = redis_client.get(key).decode('utf-8')

    if 'blocked because was identified as a crawler' in html[:200]:
        print('Delete', key)
        redis_client.delete(key)
        continue

    dest_file = join(dest_dir, key.replace('omim:', ''))
    if not isfile(dest_file):
        with open(dest_file, 'w') as f:
            f.write(html)

