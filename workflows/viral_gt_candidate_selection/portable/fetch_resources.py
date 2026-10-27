#!/usr/bin/env python3
"""Fetch checksum-pinned public Step 1-3 resources; no inference dependencies."""
import argparse,csv,hashlib,json,urllib.request,zipfile,time
from pathlib import Path
def sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def read(p):
    with Path(p).open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def obtain(url,dest,expected):
    if dest.exists() and sha(dest)==expected:return "existing_verified_cache"
    if not url.startswith(("https://zenodo.org/records/","https://raw.githubusercontent.com/LanNa1229/viral-gt-research-record/")):raise ValueError("Unexpected public source URL")
    req=urllib.request.Request(url,headers={"User-Agent":"Viral-GT-Step4-reproducibility/1.0.0"})
    # Access errors are reported; there is no authentication or access-control bypass.
    with urllib.request.urlopen(req,timeout=90) as response,dest.open("wb") as f:
        while chunk:=response.read(1024*1024):f.write(chunk)
    if sha(dest)!=expected:raise ValueError("Downloaded SHA256 mismatch "+url)
    return "downloaded_public_bytes"
def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument("--out",required=True);a.add_argument("--archive-cache");o=a.parse_args()
    root=Path(__file__).resolve().parents[1];out=Path(o.out).resolve();out.mkdir(parents=True,exist_ok=True)
    archive_out=out/"archives";archive_out.mkdir(exist_ok=True);archives={};receipt=[]
    for r in read(root/"references/public_archives.tsv"):
        src=Path(o.archive_cache)/r["filename"] if o.archive_cache else archive_out/r["filename"]
        if src.exists() and sha(src)==r["sha256"]:method="reused_checksum_verified_public_archive"
        else:src=archive_out/r["filename"];method=obtain(r["url"],src,r["sha256"])
        archives[r["archive_id"]]=src;receipt.append(dict(key=r["archive_id"],url=r["url"],sha256=sha(src),method=method))
    for r in read(root/"references/public_resources.tsv"):
        dest=out/r["key"]
        if r["archive_id"]:
            with zipfile.ZipFile(archives[r["archive_id"]]) as z:data=z.read(r["archive_member"])
            if hashlib.sha256(data).hexdigest()!=r["sha256"]:raise ValueError("Public archive member SHA256 mismatch "+r["key"])
            dest.write_bytes(data);method="extracted_verified_public_archive_member"
        else:method=obtain(r["url"],dest,r["sha256"])
        receipt.append(dict(key=r["key"],url=r["url"],archive_member=r["archive_member"],sha256=sha(dest),method=method))
    (out/"PUBLIC_FETCH_RECEIPT.json").write_text(json.dumps(dict(status="PASS",resources=len(receipt)-len(archives),archives=len(archives),files=receipt),indent=2)+"\n")
    print(json.dumps(dict(status="PASS",public_files=len(receipt)-len(archives),public_archives=len(archives),cache=str(out))))
if __name__=="__main__":main()

