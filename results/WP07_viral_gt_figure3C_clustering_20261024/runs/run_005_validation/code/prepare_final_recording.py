from pathlib import Path
OUT=Path(__file__).resolve().parents[1]
for name in ['validate.py','validate_extended.py']:
 p=OUT/'scripts'/name;lines=p.read_text(encoding='utf-8-sig').splitlines()
 for i,line in enumerate(lines):
  if line.startswith('check(all(abs(float(r['):lines[i]='check(all(abs(float(r["weight"])-(1-d[lookup[r["source"]],lookup[r["target"]]]))<1e-12 for r in edges),"independent_edge_cosine_weights")'
 p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
s=(OUT/'scripts/run_recorded.py').read_text(encoding='utf-8-sig')
s=s.replace("inputs=c.read(c.OUT/'INPUT_MANIFEST.csv')","inputs=c.read(c.OUT/'INPUT_MANIFEST.csv')\nif (c.OUT/'DERIVED_INPUT_MANIFEST.csv').exists():inputs+=c.read(c.OUT/'DERIVED_INPUT_MANIFEST.csv')")
s=s.replace("p.suffix.lower() in ('.py','.pyd','.dll')","p.suffix.lower() in ('.pyd','.dll')")
s=s.replace("api=dict(HDBSCAN_signature", "entrypoints={name:dict(path=str(module.__file__),sha256=c.sha(Path(module.__file__))) for name,module in modules.items()}\napi=dict(HDBSCAN_signature")
s=s.replace("record=dict(run_id=run_id,", "record=dict(dependency_policy='Full Python/native dependency fingerprint preserved in run_001_science; postprocessing rehashes native binaries, METADATA/RECORD and loaded module entrypoints',module_entrypoint_hashes=entrypoints,derived_input_manifest_sha256=c.sha(c.OUT/'DERIVED_INPUT_MANIFEST.csv') if (c.OUT/'DERIVED_INPUT_MANIFEST.csv').exists() else None,run_id=run_id,")
s=s.replace("if p.returncode!=0:raise SystemExit(p.returncode)\nassert record['inputs_unchanged'] and source_unchanged", "if p.returncode!=0:raise SystemExit(p.returncode)\nif not (record['inputs_unchanged'] and source_unchanged):raise SystemExit('Post-execution integrity guard failed; preserved record reports exact discrepancy.')")
(OUT/'scripts/run_final_recorded.py').write_text(s,encoding='utf-8')
