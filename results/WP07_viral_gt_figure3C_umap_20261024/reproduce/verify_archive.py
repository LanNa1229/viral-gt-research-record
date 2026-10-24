"""Verify the archive SHA256 manifest using only the Python standard library."""
import csv,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
with (root/'ARCHIVE_SHA256.tsv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
for row in rows:
 p=root/row['path'];assert p.resolve().is_relative_to(root.resolve()),row['path']
 digest=hashlib.sha256()
 with p.open('rb') as stream:
  for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
 assert p.stat().st_size==int(row['size_bytes']) and digest.hexdigest()==row['sha256'],row['path']
print('PASS:',len(rows),'archived files match size and SHA256; no scientific code was executed.')
