"""Preserve installed tool; guard absent realm/domain thresholds in a local copy."""
import pathlib,shutil,hashlib,json,difflib
OUT=pathlib.Path(__file__).resolve().parents[1]
source=OUT/'env/taxonomy/lib/python3.11/site-packages/vcontact3';dest=OUT/'compat/vendor/vcontact3'
if not dest.exists():
    dest.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(source,dest,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
p=dest/'gbga.py';original=(source/'gbga.py').read_text();text=p.read_text()
needle="                if id_realm in not_implemented:\n"
guard="                if id_realm != 'default' and id_realm not in distances[ref_domain]:\n                    logger.warning(f'No calibrated thresholds for realm {id_realm} in domain {ref_domain}; retaining it as unassigned')\n                    continue\n\n"
if guard not in text:
    assert text.count(needle)==1;text=text.replace(needle,guard+needle);p.write_text(text)
(OUT/'reports/vcontact_runtime_guard.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),fromfile='installed/vcontact3/gbga.py',tofile='local_overlay/vcontact3/gbga.py')))
(OUT/'reports/vcontact_runtime_guard.json').write_text(json.dumps(dict(version='3.1.6',installed_file_sha256=hashlib.sha256(original.encode()).hexdigest(),runtime_overlay_sha256=hashlib.sha256(text.encode()).hexdigest(),change='Skip unsupported realm/domain combinations when no calibrated threshold exists; no thresholds modified; installed package preserved'),indent=2))
print('Installed documented local missing-threshold guard; original environment unchanged')
