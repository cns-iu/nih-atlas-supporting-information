#!/usr/bin/env python3
"""Collect a versioned, updated OpenAlex citing-work snapshot by PMID.

This is intentionally a *new* snapshot, authorized by the user on 2026-10-05;
it must not be represented as the manuscript's 2026-09-24 snapshot.
"""
from __future__ import annotations
import concurrent.futures as cf
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path
import time
import os
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SHARD=os.environ.get('OPENALEX_SHARD')
TAG=f'_shard{SHARD}' if SHARD is not None else ''
OUT=ROOT/'data'/'citations'/f'openalex_updated_20261005_citing_works{TAG}.jsonl.gz'
META=ROOT/'data'/'citations'/f'openalex_updated_20261005_snapshot_metadata{TAG}.json'
DONE=ROOT/'data'/'citations'/f'openalex_updated_20261005_collection_status{TAG}.jsonl'
WORKERS=int(os.environ.get('OPENALEX_WORKERS','12')); PER_PAGE=200; RETRIES=6
LIMIT=int(os.environ.get('OPENALEX_BATCH_SIZE','0'))
N_SHARDS=int(os.environ.get('OPENALEX_N_SHARDS','1'))
API_KEY=os.environ.get('OPENALEX_API_KEY')
HEADERS={'User-Agent':'NIH-map-supporting-information/round3 (contact: research@example.org)'}
SELECT='id,publication_year,primary_topic'

def get(url):
    if API_KEY:
        url += ('&' if '?' in url else '?') + urlencode({'api_key':API_KEY})
    for attempt in range(RETRIES):
        try:
            req=Request(url,headers=HEADERS)
            with urlopen(req,timeout=60) as r: return json.loads(r.read().decode())
        except HTTPError as e:
            if e.code not in (429,500,502,503,504) or attempt==RETRIES-1: raise
        except URLError:
            if attempt==RETRIES-1: raise
        time.sleep(min(30,2**attempt))

def topic_fields(topic):
    topic=topic or {}; field=topic.get('field') or {}; domain=topic.get('domain') or {}
    return {'primary_field_id':field.get('id'),'primary_field_name':field.get('display_name'),'primary_domain_id':domain.get('id'),'primary_domain_name':domain.get('display_name')}

def one(pmid):
    try:
        focal=get(f'https://api.openalex.org/works/pmid:{pmid}?select={SELECT}')
    except HTTPError as e:
        if e.code==404: return [],{'pmid':str(pmid),'status':'no_openalex_work','http_status':404}
        raise
    f={'focal_pmid':str(pmid),'focal_work_id':focal.get('id'),'focal_pub_year':focal.get('publication_year'),**{f'focal_{k}':v for k,v in topic_fields(focal.get('primary_topic')).items()}}
    if not f['focal_work_id']: return [],{'pmid':str(pmid),'status':'no_openalex_work'}
    out=[]; cursor='*'
    while cursor:
        params={'filter':f"cites:{f['focal_work_id'].rsplit('/',1)[-1]}",'per-page':PER_PAGE,'cursor':cursor,'select':SELECT}
        page=get('https://api.openalex.org/works?'+urlencode(params));
        for w in page.get('results',[]): out.append({**f,'citing_work_id':w.get('id'),'citing_year':w.get('publication_year'),**{f'citing_{k}':v for k,v in topic_fields(w.get('primary_topic')).items()}})
        cursor=page.get('meta',{}).get('next_cursor')
        if not page.get('results'): break
    return out,{'pmid':str(pmid),'status':'ok','citing_works':len(out)}

def main():
    pmids=sorted(pd.read_csv(ROOT/'data'/'results'/'publication_functional_classification.csv',dtype={'pmid':str}).pmid.unique())
    if SHARD is not None: pmids=pmids[int(SHARD)::N_SHARDS]
    done=set()
    if DONE.exists():
        with DONE.open() as h:
            for line in h:
                try:
                    s=json.loads(line)
                    if s.get('status') in ('ok','no_openalex_work'): done.add(s['pmid'])
                except json.JSONDecodeError: pass
    todo=[x for x in pmids if x not in done]
    if LIMIT: todo=todo[:LIMIT]
    statuses=[]
    with gzip.open(OUT, 'at', encoding='utf-8') as h, DONE.open('a') as dh, cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futures={ex.submit(one,p):p for p in todo}
        for i,future in enumerate(cf.as_completed(futures),1):
            pmid=futures[future]
            try:
                rows,status=future.result(); statuses.append(status); dh.write(json.dumps(status)+'\n'); dh.flush()
                for row in rows: h.write(json.dumps(row,ensure_ascii=False)+'\n')
            except Exception as e:
                status={'pmid':pmid,'status':'failed','error':str(e)}; statuses.append(status); dh.write(json.dumps(status)+'\n'); dh.flush()
            if i%25==0: print(f'{i}/{len(todo)} focal PMIDs completed',flush=True)
    all_status=[]
    if DONE.exists():
        all_status=[json.loads(line) for line in DONE.open() if line.strip()]
    finished={s['pmid'] for s in all_status if s.get('status') in ('ok','no_openalex_work')}
    META.write_text(json.dumps({'snapshot_type':'updated_OpenAlex_snapshot','retrieved_at_utc':datetime.now(timezone.utc).isoformat(),'authorized_by_user':'2026-10-05 user request','workers':WORKERS,'source':'OpenAlex works API; citing works via filter=cites:<OpenAlex work ID>','n_target_pmids':len(pmids),'n_completed_pmids':len(finished),'n_failed_attempts':sum(s.get('status')=='failed' for s in all_status),'complete':len(finished)==len(pmids)},indent=2))
    print(f'Wrote {OUT} and {META}')
if __name__=='__main__': main()
