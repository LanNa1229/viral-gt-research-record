from common import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":10,"pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
palette=json.loads((PROJECT/"analysis/bgcprophet_reference/results/umap_all_bgc_classes_provenance.json").read_text())["colors"]
COLORS={"A":palette["NRP"],"B":palette["Polyketide"],"N":palette["Other"]}
LABELS={"A":"A: all-position mean","B":"B: core-context mean","N":"N: non-core contextual mean"}
def save(fig,name,source):
 fig.savefig(ROOT/f"figures/{name}.png",dpi=400,bbox_inches="tight")
 fig.savefig(ROOT/f"figures/{name}.pdf",bbox_inches="tight");plt.close(fig)
 with open(ROOT/"provenance/figure_sources.tsv","a") as f:f.write(name+"\t"+source+"\n")
def main():
 (ROOT/"provenance/figure_sources.tsv").write_text("")
 s=read(ROOT/"tables/A_B_N_augmentation_stability.tsv")
 fig,axs=plt.subplots(1,3,figsize=(12,4))
 for ax,metric,label in zip(axs,["mean_pairwise_cosine","min_pairwise_cosine","mean_centroid_cosine_distance"],["Mean pairwise cosine","Minimum pairwise cosine","Mean centroid cosine distance"]):
  art=ax.boxplot([s[r+"_"+metric] for r in COLORS],patch_artist=True,tick_labels=list(COLORS),whis=(5,95),
    flierprops=dict(markersize=1.3,alpha=.25),medianprops=dict(color="black",linewidth=1.3))
  for patch,color in zip(art["boxes"],COLORS.values()):patch.set_facecolor(color);patch.set_alpha(.65)
  ax.set(ylabel=label,xlabel="Representation")
 fig.suptitle("Augmentation stability: 2,502 BGCs (five states each)")
 fig.tight_layout();save(fig,"A_B_N_augmentation_stability","tables/A_B_N_augmentation_stability.tsv")
 fig,axs=plt.subplots(1,3,figsize=(12,4),sharex=True,sharey=True)
 for ax,r in zip(axs,COLORS):
  ax.scatter(s.core_gene_count,s[r+"_mean_pairwise_cosine"],s=7,alpha=.3,c=COLORS[r],linewidths=0,rasterized=True)
  ax.set(title=LABELS[r],xlabel="Core gene count",ylim=(0,1.02),xlim=(0,120))
 axs[0].set_ylabel("Mean pairwise cosine")
 fig.tight_layout();save(fig,"core_count_vs_stability","tables/A_B_N_augmentation_stability.tsv")
 c=read(ROOT/"tables/core_noncore_cancellation.tsv");b=read(ROOT/"tables/per_BGC_cancellation.tsv")
 fig,axs=plt.subplots(1,3,figsize=(14,4))
 axs[0].hist(c.deviation_cosine,bins=np.linspace(-1,1,51),color=COLORS["N"])
 axs[0].axvline(0,color="black",ls="--",lw=1)
 axs[0].set(xlabel="Cosine(core term, non-core term)",ylabel="Augmentation observations",title="Weighted deviation directions")
 axs[1].scatter(b.core_gene_count,b.mean_cross_term,s=7,alpha=.4,c=COLORS["A"],linewidths=0,rasterized=True)
 axs[1].axhline(0,color="black",ls="--",lw=1)
 axs[1].set(xlabel="Core gene count",ylabel="Per-BGC mean cross term",title="Negative = opposition; positive = reinforcement")
 axs[2].scatter(b.variance_reduction_fraction_from_cross,b.A_minus_B_stability,s=7,alpha=.4,c=COLORS["B"],linewidths=0,rasterized=True)
 axs[2].axvline(0,color="black",ls="--",lw=1);axs[2].axhline(0,color=".6",ls=":",lw=1)
 axs[2].set(xlabel="Signed variance reduction from cross term",ylabel="A minus B mean pairwise cosine",title="Positive x = cancellation")
 fig.tight_layout();save(fig,"core_noncore_deviation_relationship","tables/core_noncore_cancellation.tsv;tables/per_BGC_cancellation.tsv")
 k=read(ROOT/"tables/A_B_N_saccharide_knn_summary.tsv")
 fig,axs=plt.subplots(1,2,figsize=(12,4.5))
 for ax,definition in zip(axs,DEFS):
  for r,marker,ls in zip(COLORS,["o","s","^"],["-","-","--"]):
   g=k[(k.definition==definition)&(k.representation==r)]
   ax.plot(g.k,100*g.precision_at_k,color=COLORS[r],marker=marker,ls=ls,label=LABELS[r],markersize=5)
  ax.axhline(100*g.prevalence.iloc[0],color=".5",ls=":",label="Global prevalence")
  ax.set(title=definition.replace("_"," ").title(),xlabel="k",ylabel="Precision@k (%)",xticks=KS,ylim=(0,100))
  ax.legend(fontsize=8)
 fig.tight_layout();save(fig,"A_B_N_saccharide_knn","tables/A_B_N_saccharide_knn_summary.tsv")
 ov=read(ROOT/"tables/neighbor_overlap_summary.tsv")
 fig,axs=plt.subplots(1,3,figsize=(13,4),sharey=True)
 for ax,scope in zip(axs,["all","any_saccharide","pure_saccharide"]):
  for comparison,color,marker in [("A vs B",COLORS["B"],"s"),("A vs N",COLORS["A"],"o"),("B vs N",COLORS["N"],"^")]:
   g=ov[(ov.group==scope)&(ov.comparison==comparison)]
   ax.plot(g.k,g["mean"],color=color,marker=marker,label=comparison)
  ax.set(title=scope.replace("_"," ").title(),xlabel="k",xticks=KS,ylim=(0,1));ax.legend(fontsize=8)
 axs[0].set_ylabel("Mean neighbor-set overlap")
 fig.tight_layout();save(fig,"A_B_N_neighbor_overlap","tables/neighbor_overlap_summary.tsv")
 js("provenance/plotting.json",dict(colors=COLORS,color_source=str(PROJECT/"analysis/bgcprophet_reference/results/umap_all_bgc_classes_provenance.json"),
  convention="A blue and B orange match corrected Step 3; N uses existing neutral gray. Representation colors do not denote class identities.",
  boxplot="Median, interquartile box, 5th/95th-percentile whiskers; remaining observations shown as outliers",dpi=400))
 print("Five PNG/PDF quantitative figures saved.",flush=True)
if __name__=="__main__":main()
