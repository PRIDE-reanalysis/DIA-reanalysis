import pandas as pd
import string 
import itertools
import requests as r
import time 

from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

df = pd.read_csv("phl004_canonical_sall_osw.tsv.gz", compression='gzip', sep='\t')
ids = list(df.UniprotID.unique())
uni = list(filter(lambda x: '_' not in x,  # only proper accessions, otherwise query will fail for `accession` field 
				set(filter(lambda x: not x.startswith(tuple(string.digits)), 
			list(itertools.chain.from_iterable([id.split('/') for id in ids]))))))

# Peptides
# METHOD 0
what = df.drop_duplicates(subset='PeptideSequence', keep='first', inplace=False).apply(lambda x: SeqIO.SeqRecord(id=x['UniprotID'], seq=x['PeptideSequence']) , axis = 1)
pep_fasta = dict()
for so in what:
	nid = so.id + '-1' 
	c=1
	while nid in pep_fasta: 
		c+=1
		nid = so.id + '-' + str(c)
	so.id = nid
	so.seq = Seq(so.seq)
	so.description = "pan-human target peptides"
	pep_fasta[nid] = so
	
with open('/tmp/pan-human_peptides.fasta','w') as f:
	SeqIO.write(pep_fasta.values(), f, "fasta")


# Proteins
# METHOD 1
# this werks but dead slow, undocumented limit, brute-force explored, 10 is max
fasta = dict()
up_url = "https://rest.uniprot.org/uniprotkb/accessions?accessions={accs}&format=fasta"
#r.get(up_url.format(accs=','.join(uni))).text  # accessions endpoint has undocumented limits 
step = 10
for i in range(0, len(uni), step):
	qf = r.get(up_url.format(accs=','.join(uni[i:i+step]))).text
	fasta.update(SeqIO.to_dict(SeqIO.parse(StringIO(qf), "fasta")))
	print(len(fasta))
	time.sleep(.125)

with open('/tmp/pan-human_proteins.fasta','w') as f:
	SeqIO.write(fasta.values(), f, "fasta")
	


# METHOD 2 
# this uses a SOLR query like this
#"https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+AND+(accession:P0DSE2+OR+accession:Q61508)"
# has two sub-variants due to the size restrictions in place. According to [swagger doc](https://rest.uniprot.org/docs/?urls.primaryName=idmapping#/uniprotkb/searchCursor) limit is 500

# METHOD2 1st variant - query all, page results only works 
up_url = "https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+AND+format=fasta+AND+({accs})"
fasta = dict()

#1 queries too big will just go poof because of URL length limit
#2 issues invalid accessions will break the whole query
q = ['accession:'+ac for ac in uni]
qf = r.get(up_url.format(accs='+OR+'.join(q))).json()
num_pages = qf['last_page']

for page in range(2, num_pages + 1):
    r_sanfran = requests.get("https://api.angel.co/1/tags/1664/jobs", params={'page': page}).json()
    print r_sanfran['page']
    # TODO: extract the data

	fasta.update(SeqIO.to_dict(SeqIO.parse(StringIO(qf), "fasta")))
	print(len(fasta))


# METHOD2 2nd variant - query chunks
up_url = "https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+AND+format=fasta+AND+({accs})"
fasta = dict()
step = 500 
for i in range(0, len(uni), step):
	q = ['accession:'+ac for ac in uni[i:i+step]]
	qf = r.get(up_url.format(accs='+OR+'.join(q))).text
	fasta.update(SeqIO.to_dict(SeqIO.parse(StringIO(qf), "fasta")))
	print(len(fasta))

# METHOD2 3rd variant - 'stream'
# The stream endpoint can handle at most result sets with 5,000,000 entries
# But it craps out far earlier due to: requests.exceptions.HTTPError: 413 Client Error: Request Entity Too Large 
up_url = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=organism_id:9606+AND+({accs})"
q = ['accession:'+ac for ac in uni]
with r.get(up_url.format(accs='+OR+'.join(q)), stream=True) as rq:
    rq.raise_for_status()
    with open('stream.fasta.gz', 'wb') as f:
        for chunk in rq.iter_content(chunk_size=2**20):
            f.write(chunk)

#3rd.2
para = {
	'query':'{accs}'.format(accs='+OR+'.join(q)),
}
response = r.post('https://rest.uniprot.org/uniprotkb/stream', data=para)
#no POST!

#3rd.3
session = r.Session()
session.post('https://rest.uniprot.org/uniprotkb/stream', data=para)
# no POST!



with open('/tmp/fa.sta','w') as f:
	f.write(fasta)

with open('/tmp/fa.sta') as o:
    print(len(list(SeqIO.parse(o, "fasta"))))

with open('/tmp/accessions','w') as f:
	f.write(' '.join(uni))
	


#METHOD 4 Response chunks

import re
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = r.adapters.Retry(total=5, backoff_factor=0.125, status_forcelist=[500, 502, 503, 504])
session = r.Session()
session.mount("https://", r.adapters.HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)  # from https://www.uniprot.org/help/api_queries - ARGH hardcoded object defined outside used inside function scope
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

from itertools import zip_longest
up_url = "https://rest.uniprot.org/uniprotkb/search?fields=accession&format=fasta&query=organism_id:9606+AND+({accs})"
q = ['accession:'+ac for ac in uni]
chunk_size = 1000
fastas = list()
for q_chunk in zip_longest(*[iter(q)]*chunk_size, fillvalue=''):  # for numpy version of list chunking see QC tool code
	with r.get(up_url.format(accs='+OR+'.join(filter(lambda x: x,q_chunk)), stream=True)) as rq:
		for chunk in rq.iter_content(chunk_size=2**20):  # this is simply to chunk the response, paginationstill needs to be adressed separately
			fastas.append(chunk)  
#so this does not work

up_url = "https://rest.uniprot.org/uniprotkb/search?fields=accession&format=fasta&size={re_chunk_size}&query=organism_id:9606+AND+({accs})"
q = ['accession:'+ac for ac in uni]
query_chunk_size = 100  # uniprot is really stupid in that sense as it cant handle a query size of 1000, so that this craps out with urllib3.exceptions.MaxRetryError: HTTPSConnectionPool (Caused by ResponseError('too many 500 error responses',))
fastas = dict()
for q_chunk in zip_longest(*[iter(q)]*query_chunk_size, fillvalue=''):  # for numpy version of list chunking see QC tool code
	for batch, total in get_batch(up_url.format(accs='+OR+'.join(filter(lambda x: x,q_chunk)))): 
		fastas.update(SeqIO.to_dict(SeqIO.parse(StringIO(batch.text), "fasta")))
		print(f'{len(fastas)} / {total}')	
		
# https://rest.uniprot.org/docs/?urls.primaryName=idmapping#/
# https://biopython.org/wiki/SeqIO
# https://www.uniprot.org/help/api_queries
# https://www.uniprot.org/help/programmatic_access
# https://www.uniprot.org/help/query-fields
