#!/usr/bin/env python3
"""Later portable reconstruction, NOT the unrecovered historical repair executable."""
import argparse,csv,json
from pathlib import Path
REPAIRS={'s91124.1':'ucs91124.1','s22768.1':'ucs22768.1','p81541.1':'ucp81541.1','i09541.1':'uci09541.1'}
def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--input',required=True,type=Path);a.add_argument('--output',required=True,type=Path);args=a.parse_args()
    if args.input.resolve()==args.output.resolve() or args.output.exists():a.error('Output must be a new file, distinct from input')
    changed=0;rows=0
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.input.open() as src,args.output.open('w',newline='') as dst:
        first=src.readline();header=first.rstrip('\r\n').split('\t')[:2]==['query','target']
        if not header:src.seek(0)
        for line in src:
            values=line.rstrip('\r\n').split('\t')
            if len(values)!=17:raise ValueError('Expected17 alignment columns')
            old=values[:2];values[:2]=[REPAIRS.get(v,v) for v in old];changed+=values[:2]!=old;rows+=1;dst.write('\t'.join(values)+'\n')
    print(json.dumps(dict(rows=rows,changed_rows=changed,original_header_removed=header)))
if __name__=='__main__':main()
