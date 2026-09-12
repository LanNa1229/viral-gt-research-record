from common import *
from scipy.spatial.distance import pdist,cdist
from PIL import Image
def main():
 checks=[]
 def check(name,test,detail=""):
  checks.append(dict(check=name,passed=bool(test),detail=detail))
  if not test:
   tab("provenance/independent_validation.tsv",checks);raise AssertionError(name)
 arrays=raw();ix=np.load(ROOT/"results/group_row_indices.npy")
 meta=read(ROOT/"results/noncore_context_metadata.tsv");cm=read(ROOT/"results/consensus_metadata.tsv")
 n=meta.core_gene_count.to_numpy()[:,None]
 alternate=arrays["B"]+128/(128-n)*(arrays["A"]-arrays["B"])
 check("Independent rearrangement for all N entries",np.allclose(alternate,arrays["N"],atol=1e-12,rtol=1e-12))
 st=read(ROOT/"tables/A_B_N_augmentation_stability.tsv")
 for rep,rawarray in arrays.items():
  pair=np.array([1-pdist(rawarray[g],metric="cosine") for g in ix])
  check(rep+" SciPy all 25020 pairwise cosine means",np.allclose(pair.mean(1),st[rep+"_mean_pairwise_cosine"],atol=2e-15,rtol=0))
  check(rep+" SciPy all pairwise minima",np.allclose(pair.min(1),st[rep+"_min_pairwise_cosine"],atol=2e-15,rtol=0))
  check(rep+" pairwise medians",np.allclose(np.median(pair,axis=1),st[rep+"_median_pairwise_cosine"],atol=2e-15,rtol=0))
  check(rep+" pairwise sample SD",np.allclose(pair.std(1,ddof=1),st[rep+"_std_pairwise_cosine"],atol=2e-15,rtol=0))
  cent=np.array([sum((rawarray[i] for i in group),np.zeros(320))/5 for group in ix])
  check(rep+" exact independent consensus sum",np.array_equal(cent,np.load(ROOT/f"results/{rep}_consensus_embeddings.npy")))
  d=np.array([cdist(rawarray[group],c[None,:],metric="cosine")[:,0] for group,c in zip(ix,cent)])
  check(rep+" independent centroid distances",np.allclose(d.mean(1),st[rep+"_mean_centroid_cosine_distance"],atol=2e-15,rtol=0))
  full=cdist(cent,cent,metric="cosine");np.fill_diagonal(full,np.inf)
  neighbors=np.argsort(full,axis=1,kind="stable")[:,:50]
  check(rep+" all top50 via SciPy original cosine",np.array_equal(neighbors,np.load(ROOT/f"results/{rep}_neighbor_indices_top50.npy")))
 cancel=read(ROOT/"tables/core_noncore_cancellation.tsv")
 saved=cancel.set_index("embedding_row").loc[np.arange(12510)]
 cross=np.empty(12510);ec=np.empty(12510);en=np.empty(12510);ea=np.empty(12510)
 for group in ix:
  aa,bb,nn=(arrays[r][group] for r in ["A","B","N"])
  da=aa-aa.mean(0);db=bb-bb.mean(0);dn=nn-nn.mean(0)
  w=float(n[group[0],0]/128);ct=w*db;nt=(1-w)*dn
  for q,i in enumerate(group):
   cross[i]=2*np.dot(ct[q],nt[q]);ec[i]=np.dot(ct[q],ct[q]);en[i]=np.dot(nt[q],nt[q]);ea[i]=np.dot(da[q],da[q])
 check("All cross terms via independent dot products",np.allclose(cross,saved.cross_term,atol=1e-12,rtol=1e-12))
 check("All Euclidean variance identities",np.allclose(ea,ec+en+cross,atol=1e-12,rtol=1e-12))
 check("All cancellation metadata row identities",np.array_equal(saved.ID,meta.ID))
 counts=cm.core_gene_count.to_numpy()[:,None]
 ac,bc,nc=[np.load(ROOT/f"results/{r}_consensus_embeddings.npy") for r in ["A","B","N"]]
 check("Consensus weighted-mean identity",np.allclose(ac,counts/128*bc+(128-counts)/128*nc,atol=1e-12,rtol=1e-12))
 nn={r:np.load(ROOT/f"results/{r}_neighbor_indices_top50.npy") for r in ["A","B","N"]}
 summary=read(ROOT/"tables/A_B_N_saccharide_knn_summary.tsv")
 for row in summary.itertuples():
  mask=cm[row.definition].to_numpy(bool)
  counts=np.array([sum(mask[j] for j in nn[row.representation][i,:row.k]) for i in np.flatnonzero(mask)])
  check(f"Precision {row.representation} {row.definition} {row.k}",abs(counts.mean()/row.k-row.precision_at_k)<1e-14)
 perm=np.load(CONS/"results/shared_label_permutations.npy",mmap_mode="r")
 null=np.load(ROOT/"results/N_permutation_null_precision.npy")
 for i in list(range(25))+list(range(9975,10000)):
  for j,definition in enumerate(DEFS):
   lab=cm[definition].to_numpy(bool)[perm[i]]
   for h,k in enumerate(KS):
    q=[sum(lab[x] for x in nn["N"][i,:k])/k for i in np.flatnonzero(lab)]
    check(f"Independent N permutation {i} {definition} {k}",abs(np.mean(q)-null[i,j,h])<1e-14)
 for row in summary[summary.representation=="N"].itertuples():
  nv=null[:,DEFS.index(row.definition),KS.index(row.k)]
  p=(1+(nv>=row.precision_at_k).sum())/10001
  check(f"N empirical p {row.definition} {row.k}",abs(p-row.empirical_p)<1e-14)
 ov=read(ROOT/"tables/A_B_N_neighbor_overlap.tsv")
 for (pair,k),g in ov.groupby(["comparison","k"],sort=False):
  a,b=pair.split(" vs ");actual=[np.isin(x[:k],y[:k]).sum()/k for x,y in zip(nn[a],nn[b])]
  check(f"All overlaps {pair} {k}",np.allclose(actual,g.overlap,atol=1e-15,rtol=0))
 for p in (ROOT/"figures").glob("*.png"):
  with Image.open(p) as im:im.verify()
  check(p.name+" readable",True);check(p.stem+" PDF",(p.with_suffix(".pdf")).stat().st_size>1000)
 protected=integrity();check("All protected scientific inputs unchanged",True,str(protected))
 tab("provenance/independent_validation.tsv",checks)
 js("results/final_validation.json",dict(passed=True,checks=len(checks),protected_files_unchanged=protected))
 print(json.dumps(json.loads((ROOT/"results/final_validation.json").read_text())),flush=True)
if __name__=="__main__":main()
