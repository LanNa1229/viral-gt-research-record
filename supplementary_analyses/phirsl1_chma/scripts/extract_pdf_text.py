#!/usr/bin/env python3
"""Extract primary PDF tables with pypdf; optional root path and --bootstrap."""
from pathlib import Path
import sys,json,hashlib,datetime
from pypdf import PdfReader,__version__
args=[a for a in sys.argv[1:] if a!='--bootstrap'];bootstrap='--bootstrap' in sys.argv
root=Path(args[0]) if args else Path(__file__).resolve().parents[1]
paths=['sources/PMC10299810/publisher_1.pdf','sources/PMC9365700/41586_2022_5013_MOESM1_ESM.pdf'];checks=[]
for rel in paths:
 p=root/rel
 if not p.exists() and bootstrap:continue
 reader=PdfReader(p);text='\n\n'.join(pg.extract_text() or '' for pg in reader.pages);out=p.with_suffix('.txt');out.write_text(text,encoding='utf-8')
 checks.append(dict(input=rel,output=str(out.relative_to(root)),pages=len(reader.pages),sha256_input=hashlib.sha256(p.read_bytes()).hexdigest(),sha256_text=hashlib.sha256(text.encode()).hexdigest()))
(root/'logs/pdf_extraction.json').write_text(json.dumps(dict(pypdf_version=__version__,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),bootstrap=bootstrap,files=checks),indent=2))
print('Extracted primary table text',len(checks))
