import requests
import redis
from urllib.error import HTTPError

class dbSNP:
    def dbsnp_url(rs):
        return 'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_gene.cgi?rs=%s' % rs

    def query_dbsnp(rs):
        url = dbsnp_url(rs)
        headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
        print('%s : Query %s' % (rs, url))
        response = requests.get(url, headers)

        if response.ok:
            print(' -> OK!')
            expire_time = 60 * 60 * 24 * 30  # One month
            redis_client.setex(rs, expire_time, json.dumps(response.json()))
            return response.json()
        else:
            print(' -> %s %s' % (response.status_code, response.reason))
            return


    def dbsnp(rs, use_cache=True):
        if not redis_client.get(rs):
            query_dbsnp(rs)

        cache = redis_client.get(rs)

        if cache:
            return json.loads(cache.decode('utf8'))
