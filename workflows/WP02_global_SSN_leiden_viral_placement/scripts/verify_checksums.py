#!/usr/bin/env python3
"""Verify a sealed package against its SHA256 inventory (stdlib only)."""
import argparse,hashlib
from pathlib import Path
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True,type=Path);a=p.parse_args();root=a.root.resolve();bad=[];n=0
    for line in (root/'SHA256SUMS').read_text().splitlines():
        expected,rel=line.split('  ',1);target=(root/rel).resolve()
        if not target.is_relative_to(root):raise ValueError('Unsafe checksum path')
        h=hashlib.sha256()
        if not target.is_file():bad.append(rel);continue
        with target.open('rb') as f:
            for block in iter(lambda:f.read(1048576),b''):h.update(block)
        n+=1
        if h.hexdigest()!=expected:bad.append(rel)
    print(f'Checked {n} files; mismatches/missing: {len(bad)}')
    for rel in bad:print(rel)
    raise SystemExit(1 if bad else 0)
if __name__=='__main__':main()
