"""Frozen-weight viral BGC pilot and matched MIBiG reference construction.

Run stages in order: prepare, esm, pilot, reference. No source files are edited.
"""
import argparse
import collections
import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import pickle
import sys
import tarfile
import types

OUT = Path(__file__).resolve().parents[1]
ROOT = OUT.parents[1]
RECOVERY = ROOT / 'analysis/bgcprophet_sequence_recovery_20260912'
LABELS = ['Alkaloid', 'Terpene', 'NRP', 'Polyketide', 'RiPP', 'Saccharide', 'Other']
CHECKPOINT_SHA = 'f45ca298980d4740dd43131a7e69c99cb233368c815470572cfc46f33abf4eb5'
CLASSIFIER_SHA = 'fa4e1d045ac2c94559a9ab41054fd6b28d991d82be63e98585d7bc014564dca8'
BGC = 'S-H68_QBP05878-QBP05897'
# Accession numbers have five digits after QBP.
CORE = [f'QBP{i:05d}.1' for i in range(5878, 5898)]


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def tsv(path, rows, fields=None):
    rows = list(rows)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def read(path, delimiter='\t'):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def strategies(ordered_ids, core_ids):
    require(core_ids and len(core_ids) == len(set(core_ids)), 'Empty/duplicate core IDs')
    require(len(ordered_ids) == len(set(ordered_ids)), 'Duplicate genome protein IDs')
    positions = [ordered_ids.index(p) for p in core_ids]
    require(positions == sorted(positions), 'Core IDs not in genomic order')
    span = positions[-1] - positions[0] + 1
    status = 'available'
    real = []
    if len(ordered_ids) < 128:
        status = 'unavailable_genome_fewer_than_128_protein_coding_genes'
    elif span > 128:
        status = 'unavailable_core_span_exceeds_128'
    else:
        start = max(0, min(positions[0] - (128 - span) // 2, len(ordered_ids) - 128))
        real = ordered_ids[start:start + 128]
        require(set(core_ids).issubset(real), 'REAL_128 lost core genes')
    repeated = (core_ids * ((128 + len(core_ids) - 1) // len(core_ids)))[:128]
    return [('REAL_128', real, status), ('REPEAT_128', repeated, 'available'),
            ('CORE_ONLY', core_ids, 'available')]


def prepare():
    from Bio import SeqIO
    source = OUT / 'sources/MK016663.1.gb'
    records = list(SeqIO.parse(source, 'genbank'))
    require(len(records) == 1 and records[0].id == 'MK016663.1', 'Unexpected genome')
    record = records[0]
    genes, excluded = [], []
    for i, f in enumerate(record.features):
        if f.type != 'CDS':
            continue
        q = f.qualifiers
        if 'translation' not in q or 'protein_id' not in q:
            excluded.append(dict(feature_index=i, location=str(f.location), reason='no_translated_protein'))
            continue
        require(len(q['protein_id']) == len(q['translation']) == 1, 'Ambiguous CDS qualifiers')
        seq = q['translation'][0].replace(' ', '').replace('\n', '').upper()
        require(seq and set(seq) <= set('ACDEFGHIKLMNPQRSTVWYBXZJUO'), 'Invalid protein sequence')
        start, end = int(f.location.start), int(f.location.end)
        require(0 <= start < end <= len(record), 'Invalid coordinates')
        require(len(f.location.parts) == 1, 'Compound CDS: explicit handling required')
        translated = str(f.extract(record.seq).translate(table=int(q.get('transl_table', ['1'])[0]), cds=True))
        require(translated == seq, 'CDS translation mismatch: ' + q['protein_id'][0])
        genes.append(dict(genome_id=record.id, contig_id=record.id, organism=record.annotations['organism'],
                          protein_id=q['protein_id'][0], start_1based=start+1, end_1based=end,
                          strand=f.location.strand, source_feature_index=i, product=q.get('product', [''])[0],
                          aa_length=len(seq), sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                          sequence=seq, is_manual_core=int(q['protein_id'][0] in CORE),
                          translation_verified=True, esm_truncated=len(seq)>1022))
    genes.sort(key=lambda g: (g['start_1based'], g['end_1based'], g['source_feature_index']))
    ids = [g['protein_id'] for g in genes]
    require(set(CORE).issubset(ids), 'Missing requested proteins')
    require([p for p in ids if p in CORE] == CORE, 'Requested range not in genomic order')
    for i, g in enumerate(genes):
        g['gene_index_0based'] = i
    core = [g for g in genes if g['is_manual_core']]
    require(len(core) == 20, 'Expected exactly 20 manually specified genes')
    tsv(OUT/'genome_protein_manifest.tsv', genes)
    tsv(OUT/'core_protein_manifest.tsv', core)
    tsv(OUT/'excluded_cds.tsv', excluded, ['feature_index', 'location', 'reason'])
    tsv(OUT/'viral_bgc_manifest.tsv', [dict(bgc_id=BGC, genome_id=record.id, contig_id=record.id,
        organism=record.annotations['organism'], genome_description=record.description,
        genome_completeness='partial genome' if 'partial' in record.description.lower() else 'see source record',
        ordered_core_protein_ids=' '.join(CORE),
        core_gene_count=len(CORE), genome_gene_count=len(genes), genome_length_nt=len(record),
        core_start_1based=min(g['start_1based'] for g in core), core_end_1based=max(g['end_1based'] for g in core),
        manual_expected_function='glycan biosynthesis; expected Saccharide',
        manual_annotation_source='User task and explicit inclusive QBP05878.1-QBP05897.1 range, 2026-10-16',
        membership_source='user_manual_selection; no model predictions used',
        source_genbank=str(source.relative_to(OUT)), source_sha256=sha(source))])
    plan, positions = [], []
    for name, selected, status in strategies(ids, CORE):
        counts = collections.Counter(selected)
        mask = [int(p in CORE) for p in selected]
        plan.append(dict(bgc_id=BGC, strategy=name, status=status, synthetic=name=='REPEAT_128',
                         gene_count=len(selected), unique_gene_count=len(counts), core_gene_count=len(CORE),
                         core_positions=sum(mask), ordered_protein_ids=' '.join(selected),
                         membership_mask=' '.join(map(str, mask)), gene_copy_counts=json.dumps(dict(counts)),
                         attention_mask='float32 additive 1=core,0=flank; not boolean exclusion',
                         truncated_core_list=len(CORE)>128 and name=='REPEAT_128'))
        for i, p in enumerate(selected):
            positions.append(dict(bgc_id=BGC, strategy=name, position_0based=i, protein_id=p,
                                  membership=mask[i], copies_in_input=counts[p], genome_gene_index_0based=ids.index(p)))
    tsv(OUT/'input_strategy_manifest.tsv', plan)
    tsv(OUT/'input_positions.tsv', positions)
    with (OUT/'sources/genome_proteins.faa').open('w') as f:
        for g in genes:
            f.write('>'+g['protein_id']+'\n'+g['sequence']+'\n')
    save_json(OUT/'manifest_validation.json', dict(status='PASS', n_core=len(core), n_genome_proteins=len(genes),
        missing_core_ids=[], duplicate_protein_ids=len(ids)-len(set(ids)),
        duplicate_core_sequences=len(core)-len({g['sequence_sha256'] for g in core}),
        all_cds_translations_match=True, excluded_CDS=len(excluded),
        coordinate_convention='1-based inclusive; genomic ascending order regardless of strand'))
    print(json.dumps(read(OUT/'viral_bgc_manifest.tsv'), indent=2), flush=True)


def torch_setup():
    import torch
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    torch.manual_seed(1729)
    return torch


def esm_stage():
    import numpy as np
    torch = torch_setup()
    import esm
    genes = read(OUT/'genome_protein_manifest.tsv')
    needed = set(' '.join(r['ordered_protein_ids'] for r in read(OUT/'input_strategy_manifest.tsv')).split())
    selected = [g for g in genes if g['protein_id'] in needed]
    checkpoint = RECOVERY/'model_cache/hub/checkpoints/esm2_t6_8M_UR50D.pt'
    before = sha(checkpoint)
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(checkpoint))
    model.eval().float()
    converter = alphabet.get_batch_converter(truncation_seq_length=1022)
    cache, rows, vectors = {}, [], []
    with torch.inference_mode():
        for i, g in enumerate(selected):
            key = g['sequence_sha256']
            if key not in cache:
                seq = g['sequence']
                _, _, tokens = converter([(key, seq)])
                h = model(tokens, repr_layers=[6], return_contacts=False)['representations'][6]
                cache[key] = h[0, 1:min(len(seq), 1022)+1].mean(0).cpu().numpy()
            v = cache[key]
            require(v.shape == (320,) and np.isfinite(v).all(), 'Invalid ESM vector')
            vectors.append(v)
            rows.append(dict(embedding_row=i, protein_id=g['protein_id'], sequence_sha256=key,
                             aa_length=g['aa_length'], encoded_length=min(int(g['aa_length']),1022),
                             truncated=int(g['aa_length'])>1022))
            if i % 16 == 0:
                print('ESM', i, '/', len(selected), flush=True)
    np.save(OUT/'protein_esm2_embeddings.npy', np.stack(vectors))
    tsv(OUT/'protein_embedding_metadata.tsv', rows)
    require(sha(checkpoint)==before, 'ESM checkpoint changed')
    save_json(OUT/'esm_provenance.json', dict(checkpoint=str(checkpoint), sha256=before,
        model='esm2_t6_8M_UR50D', layer=6, pooling='residue mean excluding BOS/EOS', truncation=1022,
        unique_sequences=len(cache), n_proteins=len(rows), torch=torch.__version__, esm=esm.__version__,
        dtype='float32', device='cpu', checkpoint_unchanged=True))


def load_model(torch):
    checkpoint = OUT/'sources/classifier.pt'
    if not checkpoint.exists():
        with tarfile.open(OUT/'sources/model.tar.gz') as tf:
            members = [m for m in tf.getmembers() if m.isfile() and Path(m.name).name=='classifier.pt']
            require(len(members)==1, 'Ambiguous archive checkpoint')
            checkpoint.write_bytes(tf.extractfile(members[0]).read())
    require(sha(checkpoint)==CHECKPOINT_SHA, 'Checkpoint differs from validated HPC checkpoint')
    src = ROOT/'sources/upstream/bgc_prophet/train'
    require(sha(src/'classifier.py')==CLASSIFIER_SHA, 'Classifier source differs from validated HPC source')
    # Load exact source modules without unrelated package __init__ imports or package installation.
    package = types.ModuleType('frozen_bgc_train')
    package.__path__ = [str(src)]
    sys.modules['frozen_bgc_train'] = package
    spec = importlib.util.spec_from_file_location('frozen_bgc_train.classifier', src/'classifier.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    model = module.transformerClassifier(d_model=320, nhead=5, num_encoder_layers=2, max_len=128,
        dim_feedforward=1280, labels_num=7, transformer_dropout=.1, mlp_dropout=.5)
    model.load_state_dict(torch.load(checkpoint, map_location='cpu', weights_only=True), strict=True)
    model.eval()
    return model


def infer(torch, model, vectors, mask, check_hook=False):
    import numpy as np
    x = torch.tensor(np.asarray(vectors), dtype=torch.float32).unsqueeze(0)
    m = torch.tensor(mask, dtype=torch.float32).unsqueeze(0)
    require(x.shape == (1,len(mask),320) and torch.isfinite(x).all(), 'Invalid classifier input')
    captured = {}
    # Keep autograd enabled as in the validated original classify() call.
    if check_hook:
        baseline = model(src=x, src_key_padding_mask=m).detach().clone()
    h1 = model.encoder.register_forward_hook(lambda _, args, result: captured.update(h=result.detach().clone()))
    h2 = model.classifier.register_forward_pre_hook(lambda _, args: captured.update(pooled=args[0].detach().clone()))
    try:
        scores = model(src=x, src_key_padding_mask=m).detach()
    finally:
        h1.remove()
        h2.remove()
    h, pooled = captured['h'], captured['pooled']
    require(h.shape==x.shape and pooled.shape==(1,320), 'Unexpected hidden tensor shape')
    require(torch.equal(h.mean(1),pooled), 'Hook differs from published mean pooling')
    if check_hook:
        require(torch.equal(baseline,scores), 'Hook changed classification scores')
    core = h[0, m[0].bool()].mean(0)
    result = [scores[0].numpy(), pooled[0].numpy(), core.numpy()]
    require(all(np.isfinite(v).all() for v in result), 'Nonfinite classifier output')
    return result


def predicted(scores):
    threshold = max(float(scores[6]), .5)
    labels = [name for name, score in zip(LABELS, scores) if score >= threshold]
    if 'Other' in labels and len(labels)>1:
        labels.remove('Other')
    return ' '.join(labels) or 'Unknown', threshold


def pilot():
    import numpy as np
    torch = torch_setup()
    model = load_model(torch)
    original = {k:v.clone() for k,v in model.state_dict().items()}
    protein_meta = read(OUT/'protein_embedding_metadata.tsv')
    protein_vectors = np.load(OUT/'protein_esm2_embeddings.npy', allow_pickle=False)
    require(protein_vectors.shape==(len(protein_meta),320), 'ESM alignment failure')
    lookup = {r['protein_id']:protein_vectors[i] for i,r in enumerate(protein_meta)}
    for i,r in enumerate(protein_meta):
        require(int(r['embedding_row'])==i, 'ESM metadata row mismatch')
    rows, embeddings, core_embeddings, metadata, qc = [], [], [], [], []
    for strategy in read(OUT/'input_strategy_manifest.tsv'):
        if strategy['status'] != 'available':
            continue
        name = strategy['strategy']
        ids = strategy['ordered_protein_ids'].split()
        mask = list(map(int, strategy['membership_mask'].split()))
        require(name=='CORE_ONLY' or len(ids)==128, 'Expected 128 input positions')
        scores, pooled, core = infer(torch, model, [lookup[p] for p in ids], mask, True)
        labels, threshold = predicted(scores)
        row = dict(bgc_id=strategy['bgc_id'], strategy=name,
                   **{f'score_{k}':float(v) for k,v in zip(LABELS,scores)},
                   saccharide_score=float(scores[5]), saccharide_rank=1+int((scores>scores[5]).sum()),
                   applied_threshold=threshold, predicted_labels=labels,
                   saccharide_predicted='Saccharide' in labels.split(),
                   score_vector_json=json.dumps(scores.tolist()), score_type='raw model sigmoid outputs',
                   isBGC_reporting_gate='bypassed_for_manually_curated_locus')
        rows.append(row)
        metadata.append(dict(embedding_row=len(embeddings), bgc_id=strategy['bgc_id'],strategy=name,
                             gene_count=len(ids), unique_gene_count=len(set(ids)), core_gene_count=int(strategy['core_gene_count']),
                             core_position_count=sum(mask),
                             manual_expected_function='glycan biosynthesis; expected Saccharide',
                             pooling='published all-position mean immediately before classifier',
                             synthetic=name=='REPEAT_128'))
        embeddings.append(pooled)
        core_embeddings.append(core)
        qc.append(dict(bgc_id=strategy['bgc_id'],strategy=name,input_shape=str((1,len(ids),320)),
                       hook_scores_bitwise_equal=True,pooled_matches_classifier_input=True,finite=True))
    # Exact uniform repetition invariance: use actual pilot proteins, with 1,2,4,8,16-core controls.
    equivalence=[]
    for n in [1,2,4,8,16]:
        ids = CORE[:n]
        b = infer(torch,model,[lookup[p] for p in (ids*(128//n))],[1]*128)
        c = infer(torch,model,[lookup[p] for p in ids],[1]*n)
        passed = all(np.allclose(u,v,atol=1e-5,rtol=1e-5) for u,v in zip(b,c))
        equivalence.append(dict(n=n,control='prefix of pilot core; numerical control only',
            score_max_abs=float(np.max(np.abs(b[0]-c[0]))),embedding_max_abs=float(np.max(np.abs(b[1]-c[1]))),
            passed=passed,atol=1e-5,rtol=1e-5))
        require(passed, 'Uniform repetition numerical equivalence failed')
    require(all(torch.equal(v,model.state_dict()[k]) for k,v in original.items()), 'Model weights changed')
    require(sha(OUT/'sources/classifier.pt')==CHECKPOINT_SHA, 'Checkpoint changed')
    tsv(OUT/'classification_scores.tsv',rows)
    tsv(OUT/'embedding_metadata.tsv',metadata)
    np.save(OUT/'viral_bgc_embeddings.npy',np.stack(embeddings))
    np.save(OUT/'viral_bgc_core_position_embeddings.npy',np.stack(core_embeddings))
    tsv(OUT/'pilot_quality_controls.tsv',qc)
    tsv(OUT/'uniform_repetition_controls.tsv',equivalence)
    save_json(OUT/'pilot_validation.json',dict(status='PASS',n_bgcs=1,n_strategies=len(rows),
        checkpoint_sha256=CHECKPOINT_SHA,classifier_source_sha256=CLASSIFIER_SHA,weights_unchanged=True,
        all_hook_probabilities_bitwise_equal=True,catovirus_regression='NOT_AVAILABLE: historical query outputs not located',
        attention_mask='float32 0/1 additive key attention bias',positional_encoding='computed but unused, original forward',
        core_length=20,repeat_core_counts='first 8 genes have 7 copies; remaining 12 have 6',
        primary_pooling='all-position mean',supplementary_pooling='core-position mean',torch=torch.__version__))
    print(json.dumps(rows,indent=2),flush=True)


class SafeTensorUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        import torch
        if (module,name)==('torch._utils','_rebuild_tensor_v2'):
            return torch._utils._rebuild_tensor_v2
        if (module,name)==('torch.storage','_load_from_bytes'):
            return lambda b:torch.load(io.BytesIO(b),map_location='cpu',weights_only=True)
        if (module,name)==('collections','OrderedDict'):
            return collections.OrderedDict
        raise pickle.UnpicklingError('Forbidden global '+module+'.'+name)


def reference(limit=0):
    import numpy as np
    import lmdb
    require(json.loads((OUT/'pilot_validation.json').read_text())['status']=='PASS','Pilot must pass first')
    torch=torch_setup()
    model=load_model(torch)
    train=read(ROOT/'BGC_train_dataset_classify.csv',',')
    require(len(train)==12510,'Unexpected training population')
    grouped=collections.defaultdict(list)
    for r in train:
        grouped[r['ID'].rsplit('_',1)[0]].append(r)
    require(len(grouped)==2502,'Expected 2502 biological BGCs')
    metadata=read(ROOT/'analysis/bgcprophet_reference/results/consensus_metadata.tsv')
    require({r['BGC'] for r in metadata}==set(grouped),'Reference population mismatch')
    db=RECOVERY/'lmdb_snapshot/data.mdb'
    stat_before=(db.stat().st_size,db.stat().st_mtime_ns)
    env=lmdb.open(str(db.parent),readonly=True,lock=False,create=False,readahead=False)
    refs={'REPEAT_128':[],'CORE_ONLY':[]}
    reference_rows=[]
    checks=[]
    old=np.load(ROOT/'analysis/bgcprophet_reference/data/mibig_12510_embeddings.npy',mmap_mode='r')
    old_scores=np.load(ROOT/'analysis/bgcprophet_reference/data/mibig_12510_probabilities.npy',mmap_mode='r')
    oldmeta=read(ROOT/'analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv')
    require(len(oldmeta)==len(train),'Original row alignment failure')
    for i,(r,m) in enumerate(zip(train,oldmeta)):
        require(int(m['embedding_row'])==i and all(r[k]==m[k] for k in ['ID','labels','isBGC']), 'Original row identity mismatch')
    def fetch(txn,p):
        payload=txn.get(p.encode('ascii'))
        require(payload is not None,'Missing LMDB protein '+p)
        obj=SafeTensorUnpickler(io.BytesIO(payload)).load()
        require(obj['label']==p,'LMDB identity mismatch')
        v=obj['mean_representations'][6]
        require(tuple(v.shape)==(320,) and torch.isfinite(v).all(),'Invalid LMDB vector')
        return v.numpy()
    with env.begin() as txn:
        # Cross-runtime regression on 20 original 128-position inputs before any reference generation.
        test_indices=np.linspace(0,len(train)-1,20,dtype=int)
        for i in test_indices:
            r=train[i]
            scores,emb,_=infer(torch,model,[fetch(txn,p) for p in r['TDsentence'].split()],
                               [int(float(v)) for v in r['TDlabels'].split()],True)
            ep=np.max(np.abs(emb-old[i])); sp=np.max(np.abs(scores-old_scores[i]))
            passed=np.allclose(emb,old[i],atol=1e-4,rtol=1e-4) and np.allclose(scores,old_scores[i],atol=1e-4,rtol=1e-4)
            checks.append(dict(source_embedding_row=int(i),ID=r['ID'],embedding_max_abs=float(ep),
                               score_max_abs=float(sp),passed=bool(passed),atol=1e-4,rtol=1e-4))
        tsv(OUT/'mibig_original_regression.tsv',checks)
        require(all(r['passed'] for r in checks),'Original reference regression failed; do not scale')
        for i,m in enumerate(metadata[:limit or None]):
            group=grouped[m['BGC']]
            require(len(group)==5 and {r['ID'].rsplit('_',1)[1] for r in group}=={'-1','0','1','2','3'},'Augmentation mismatch')
            require(len({r['sentence'] for r in group})==len({r['labels'] for r in group})==1,'Core/label varies across augmentations')
            require(group[0]['labels']==m['label'],'Consensus labels differ from source')
            core=group[0]['sentence'].split()
            require(len(core)==len(set(core)) and len(core)>0,'Invalid reference core IDs')
            base=[fetch(txn,p) for p in core]
            for strategy in refs:
                vectors=(base*((128+len(base)-1)//len(base)))[:128] if strategy=='REPEAT_128' else base
                _,e,_=infer(torch,model,vectors,[1]*len(vectors))
                refs[strategy].append(e)
            reference_rows.append(dict(embedding_row=i,BGC=m['BGC'],label=m['label'],
                core_gene_count=len(core),n_source_augmentations=5,n_biological_reference_rows=1,
                repeated_core_truncated=len(core)>128))
            if i%100==0: print('Reference',i,'/2502',flush=True)
    env.close()
    require(stat_before==(db.stat().st_size,db.stat().st_mtime_ns),'LMDB modified')
    suffix='_pilot' if limit else ''
    for strategy,values in refs.items():
        np.save(OUT/f'mibig_{strategy}{suffix}_embeddings.npy',np.stack(values))
    tsv(OUT/f'mibig_matched_reference_metadata{suffix}.tsv',reference_rows)
    save_json(OUT/f'reference_validation{suffix}.json',dict(status='PASS',n_biological_bgcs=len(reference_rows),
        lmdb_readonly=True,lmdb_stat_unchanged=True,original_regression_passed=True,
        original_regression_atol=1e-4,original_regression_rtol=1e-4,
        repeated_core_truncation_count=sum(r['repeated_core_truncated'] for r in reference_rows)))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['prepare','esm','pilot','reference'])
    parser.add_argument('--limit',type=int,default=0)
    args=parser.parse_args()
    if args.stage=='reference': reference(args.limit)
    else: {'prepare':prepare,'esm':esm_stage,'pilot':pilot}[args.stage]()
