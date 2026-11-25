"""Exact cosine graphs, RB-configuration Leiden and independent precomputed HDBSCAN."""
import collections,itertools,time
import common as c
import numpy as np
import igraph as ig
import leidenalg as la
from scipy.stats import skew
from sklearn.cluster import HDBSCAN
from sklearn.metrics import adjusted_rand_score,silhouette_score

def canonical(membership,ids,prefix):
 groups=collections.defaultdict(list)
 for i,g in enumerate(membership):
  if int(g)>=0:groups[int(g)].append(ids[i])
 mapping={g:f'{prefix}{i+1:02d}' for i,(g,_) in enumerate(sorted(groups.items(),key=lambda kv:min(kv[1])))}
 return np.array([mapping.get(int(g),'Noise') for g in membership])

def distances(x):
 raw=np.asarray(x);u=raw.astype(np.float64);u/=np.linalg.norm(u,axis=1,keepdims=True)
 sim=np.clip(u@u.T,-1,1);sim=(sim+sim.T)/2;np.fill_diagonal(sim,1)
 _,inv=np.unique(raw,axis=0,return_inverse=True)
 identical=inv[:,None]==inv[None,:];sim[identical]=1
 d=np.clip(1-sim,0,2);np.fill_diagonal(d,0)
 assert np.isfinite(d).all() and np.array_equal(d,d.T) and np.all(d>=0)
 return sim,d,inv

def graph(sim,d,ids,k,inv):
 n=len(ids);ids=np.asarray(ids);neighbors=[];indegree=np.zeros(n,int);edges=set();ties=0
 for i in range(n):
  order=np.lexsort((ids,d[i]));order=order[order!=i];chosen=order[:k];neighbors.append(chosen);indegree[chosen]+=1
  ties+=int(len(order)>k and d[i,order[k-1]]==d[i,order[k]])
  edges.update((min(i,int(j)),max(i,int(j))) for j in chosen)
 pairs=sorted((i,j) for i,j in edges if sim[i,j]>0);weights=[float(sim[i,j]) for i,j in pairs]
 g=ig.Graph(n=n,edges=pairs,directed=False);g.es['weight']=weights
 components=g.connected_components();degrees=np.asarray(g.degree());counts=np.bincount(inv)
 nearest=np.array([d[i,ns[0]] for i,ns in enumerate(neighbors)]);kth=np.array([d[i,ns[-1]] for i,ns in enumerate(neighbors)])
 diagnostics=dict(n_nodes=n,k=k,n_edges=g.ecount(),density=g.density(loops=False),n_components=len(components),component_sizes=sorted(components.sizes(),reverse=True),isolated_nodes=int(np.sum(degrees==0)),min_degree=int(degrees.min()),max_degree=int(degrees.max()),mean_degree=float(degrees.mean()),max_knn_indegree=int(indegree.max()),indegree_skew=float(skew(indegree)) if np.std(indegree)>0 else 0.0,indegree_gini=float(np.abs(indegree[:,None]-indegree).sum()/(2*n*indegree.sum())),exact_vector_groups=int(len(counts)),duplicate_vector_loci=int(counts[counts>1].sum()),duplicate_pair_edges=sum(inv[i]==inv[j] for i,j in pairs),duplicate_edge_fraction=sum(inv[i]==inv[j] for i,j in pairs)/len(pairs) if pairs else 0.0,nonpositive_union_edges_omitted=len(edges)-len(pairs),nodes_tied_at_k_boundary=ties,minimum_retained_weight=min(weights) if weights else None,maximum_retained_weight=max(weights) if weights else None)
 for label,values in [('pairwise_cosine',sim[np.triu_indices(n,1)]),('nearest_distance',nearest),('kth_neighbor_distance',kth)]:
  diagnostics.update({label+'_'+q:float(v) for q,v in zip(['min','q05','q25','median','q75','q95','max'],np.quantile(values,[0,.05,.25,.5,.75,.95,1]))})
 nodes=[dict(sample_id=ids[i],degree=int(degrees[i]),knn_indegree=int(indegree[i]),component=int(components.membership[i]),nearest_distance=float(nearest[i]),kth_distance=float(kth[i]),exact_vector_multiplicity=int(counts[inv[i]])) for i in range(n)]
 return g,diagnostics,nodes,neighbors

def leiden(g,ids,resolution,seed,prefix):
 cfg=c.CONFIG['leiden'];p=la.RBConfigurationVertexPartition(g,initial_membership=list(range(g.vcount())),weights='weight',resolution_parameter=resolution)
 opt=la.Optimiser();opt.set_rng_seed(int(seed));opt.max_comm_size=cfg['max_comm_size'];improvements=[]
 for step in range(cfg['iteration_limit']):
  delta=float(opt.optimise_partition(p,n_iterations=cfg['iterations_per_call']));improvements.append(delta)
  if delta<=cfg['stop_improvement_at_most']:break
 labels=canonical(p.membership,ids,prefix)
 info=dict(n_communities=len(set(labels)),community_sizes=sorted(collections.Counter(labels).values(),reverse=True),rb_objective=float(p.quality()),weighted_modularity_gamma1=float(g.modularity(p.membership,weights='weight',resolution=1.0)),actual_optimizer_calls=len(improvements),iterations_per_call=cfg['iterations_per_call'],converged=improvements[-1]<=cfg['stop_improvement_at_most'],improvements=improvements)
 return labels,info

def run_cohort(name,x,ids,prefix,hdbscan=True):
 print('COHORT',name,len(ids),flush=True);sim,d,inv=distances(x);runs={};summaries=[];graphrows=[];allnodes=[];run_details=[]
 for k in c.CONFIG['grid']['k']:
  g,diag,nodes,neighbors=graph(sim,d,ids,k,inv);graphrows.append(dict(cohort=name,**diag))
  allnodes.extend(dict(cohort=name,k=k,**row) for row in nodes)
  if k==c.CONFIG['baseline']['k']:
   baseline_graph=g
   c.write(f'graphs/{name}_k20_edges.csv',[dict(source=ids[e.source],target=ids[e.target],weight=float(e['weight']),duplicate_vector_edge=bool(inv[e.source]==inv[e.target])) for e in g.es])
   if name=='selected287':
    c.write('graphs/selected287_directed_k20_neighbors.csv',[dict(sample_id=ids[i],rank=rank,neighbor_id=ids[j],cosine=float(sim[i,j]),distance=float(d[i,j]),exact_duplicate=bool(inv[i]==inv[j])) for i,ns in enumerate(neighbors) for rank,j in enumerate(ns,1)])
  for resolution in c.CONFIG['grid']['resolution']:
   seeds=list(c.CONFIG['grid']['seeds'])
   if k==20 and resolution==1.0:seeds.append(42)
   for seed in seeds:
    run=f'k{k}_r{resolution:g}_s{seed}';labels,info=leiden(g,ids,resolution,seed,prefix);runs[run]=labels
    run_details.append(dict(cohort=name,run_id=run,k=k,resolution=resolution,seed=seed,**info))
    summaries.append(dict(cohort=name,run_id=run,k=k,resolution=resolution,seed=seed,**{key:value for key,value in info.items() if key!='improvements'}))
 baseline=runs['k20_r1_s42']
 for row in summaries:row['ARI_to_cohort_baseline']=float(adjusted_rand_score(baseline,runs[row['run_id']]))
 c.write(f'assignments/{name}_leiden_grid.csv',[dict(sample_id=sid,**{run:lab[i] for run,lab in runs.items()}) for i,sid in enumerate(ids)])
 for row in summaries:row['community_sizes']=';'.join(map(str,row['community_sizes']))
 c.write(f'stability/{name}_grid_summary.csv',summaries);c.save(f'provenance/{name}_optimizer_iterations.json',run_details)
 for row in graphrows:row['component_sizes']=';'.join(map(str,row['component_sizes']))
 c.write(f'graphs/{name}_diagnostics.csv',graphrows);c.write(f'graphs/{name}_node_QC.csv',allnodes)
 grid=[r for r in runs if r!='k20_r1_s42'];seedpairs=[]
 for k in c.CONFIG['grid']['k']:
  for res in c.CONFIG['grid']['resolution']:
   rs=[f'k{k}_r{res:g}_s{s}' for s in c.CONFIG['grid']['seeds']];values=[adjusted_rand_score(runs[a],runs[b]) for a,b in itertools.combinations(rs,2)]
   seedpairs.append(dict(cohort=name,k=k,resolution=res,pairwise_seed_ARI_min=float(min(values)),pairwise_seed_ARI_median=float(np.median(values)),pairwise_seed_ARI_max=float(max(values))))
 c.write(f'stability/{name}_seed_ARI.csv',seedpairs)
 if name=='selected287':
  pairwise=np.array([[adjusted_rand_score(runs[a],runs[b]) for b in grid] for a in grid]);c.write('stability/selected287_all_grid_ARI.csv',[dict(run_id=run,**{other:float(pairwise[i,j]) for j,other in enumerate(grid)}) for i,run in enumerate(grid)])
  co=np.mean(np.stack([labels[:,None]==labels[None,:] for key,labels in runs.items() if key!='k20_r1_s42']),axis=0)
  c.write('stability/selected287_coassignment.csv',[dict(sample_id=sid,**{other:float(co[i,j]) for j,other in enumerate(ids)}) for i,sid in enumerate(ids)])
 hdb={};hrows=[]
 if hdbscan:
  for mcs in c.CONFIG['hdbscan_grid']['min_cluster_size']:
   for ms in c.CONFIG['hdbscan_grid']['min_samples']:
    parameters={k:v for k,v in c.CONFIG['hdbscan'].items() if k not in ('implementation','min_samples_convention')};parameters.update(min_cluster_size=mcs,min_samples=ms)
    before=d.copy();model=HDBSCAN(**parameters).fit(d);assert np.array_equal(d,before),'HDBSCAN mutated original distance matrix'
    labels=canonical(model.labels_,ids,'HC');probs=model.probabilities_;noise=labels=='Noise';keep=~noise
    hkey=f'mcs{mcs}_ms{ms}';hdb[hkey]=(labels,probs)
    hrows.append(dict(cohort=name,run_id=hkey,min_cluster_size=mcs,min_samples=ms,n_clusters=len(set(labels)-{'Noise'}),cluster_sizes=';'.join(str(v) for _,v in sorted(collections.Counter(labels[keep]).items())),noise_n=int(noise.sum()),noise_fraction=float(noise.mean()),ARI_to_Leiden_including_noise=float(adjusted_rand_score(baseline,labels)),ARI_to_Leiden_nonnoise=float(adjusted_rand_score(baseline[keep],labels[keep])) if keep.sum()>1 else None,nonnoise_coverage=float(keep.mean())))
    c.save(f'provenance/{name}_{hkey}_HDBSCAN_parameters.json',model.get_params())
  c.write(f'assignments/{name}_HDBSCAN_grid.csv',[dict(sample_id=sid,**{key:labels[i] for key,(labels,probs) in hdb.items()},**{key+'_probability':float(probs[i]) for key,(labels,probs) in hdb.items()}) for i,sid in enumerate(ids)])
  c.write(f'stability/{name}_HDBSCAN_summary.csv',hrows)
 silhouette=float(silhouette_score(d,baseline,metric='precomputed')) if 1<len(set(baseline))<len(ids) else None
 print('DONE',name,'baseline',collections.Counter(baseline),'silhouette',silhouette,flush=True)
 return dict(name=name,x=x,ids=list(ids),sim=sim,d=d,inv=inv,graph=baseline_graph,runs=runs,baseline=baseline,grid_summary=summaries,seed_ARI=seedpairs,hdb=hdb,hdb_summary=hrows,silhouette=silhouette,graph_diagnostics=graphrows)
