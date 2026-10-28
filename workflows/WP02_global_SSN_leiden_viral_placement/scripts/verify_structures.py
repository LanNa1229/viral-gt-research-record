#!/usr/bin/env python3
"""Verify original coordinate ZIP members against the full per-file manifest."""
import argparse,csv,hashlib,json,zipfile
from pathlib import Path
def ensure_case_sensitive_directory(folder):
    """Reject case-insensitive extraction before any original member is written."""
    folder.mkdir(parents=True)
    lower=folder/'.structure_case_probe_a';upper=folder/'.STRUCTURE_CASE_PROBE_A'
    created=[]
    try:
        with lower.open('x') as f:f.write('lower')
        created.append(lower)
        try:
            with upper.open('x') as f:f.write('upper')
            created.append(upper)
        except FileExistsError:
            raise ValueError('Extraction requires a CASE-SENSITIVE filesystem: 27 original filename pairs differ only by case. Use native WSL/Linux storage, not a default Windows-mounted directory.')
    finally:
        for path in created:
            if path.exists():path.unlink()
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-root',required=True,type=Path);p.add_argument('--extract-dir',type=Path);p.add_argument('--report',type=Path);a=p.parse_args()
    folder=a.data_root/'inputs/structures';rows=read(folder/'structures_manifest.tsv');expected={r['relative_filepath']:r for r in rows};ids=[r['protein_id'] for r in rows];metadata=read(a.data_root/'inputs/metadata/005_GTA_node_metadata_completed_20260727.tsv')
    if len(rows)!=4854 or len(expected)!=4854 or len(set(ids))!=4854:raise ValueError('Coordinate count/ID uniqueness failure')
    if set(ids)!={r['name'] for r in metadata}:raise ValueError('Coordinate/node ID mismatch')
    if a.extract_dir and a.extract_dir.exists():p.error('Use a new extraction directory')
    if a.extract_dir:ensure_case_sensitive_directory(a.extract_dir)
    checked=0;total=0
    with zipfile.ZipFile(folder/'GT_A_4854_input_structures.zip') as z:
        files=[i for i in z.infolist() if not i.is_dir()]
        if len(files)!=4854 or len({i.filename for i in files})!=4854 or set(expected)!={i.filename for i in files}:raise ValueError('Coordinate archive member mismatch')
        if z.testzip() is not None:raise ValueError('ZIP CRC failure')
        for info in files:
            r=expected[info.filename];h=hashlib.sha256();n=0;output=None
            if a.extract_dir:
                target=(a.extract_dir/info.filename).resolve()
                if not target.is_relative_to(a.extract_dir.resolve()):raise ValueError('Unsafe archive member')
                target.parent.mkdir(parents=True,exist_ok=True);output=target.open('wb')
            try:
                with z.open(info) as src:
                    for block in iter(lambda:src.read(1048576),b''):
                        h.update(block);n+=len(block)
                        if output:output.write(block)
            finally:
                if output:output.close()
            if n!=int(r['file_size_bytes']) or h.hexdigest()!=r['sha256']:raise ValueError('Coordinate byte mismatch: '+info.filename)
            checked+=1;total+=n
    result=dict(coordinate_files_verified=checked,total_uncompressed_bytes=total,node_metadata_IDs_equal=True,CRC_pass=True,all_SHA256_equal=True,extracted_to=str(a.extract_dir) if a.extract_dir else None)
    if a.report:a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
