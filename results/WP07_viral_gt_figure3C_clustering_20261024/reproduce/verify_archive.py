"""Verify archive bytes only. Does not execute scientific code."""
import csv,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
with (root/'ARCHIVE_SHA256.tsv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
for row in rows:
 p=root/row['path'];assert p.resolve().is_relative_to(root.resolve()),row['path'];h=hashlib.sha256()
 with p.open('rb') as stream:
  for block in iter(lambda:stream.read(1048576),b''):h.update(block)
 assert p.stat().st_size==int(row['size_bytes']) and h.hexdigest()==row['sha256'],row['path']
print('PASS',len(rows),'archived file sizes and SHA256 values; no science executed.')
