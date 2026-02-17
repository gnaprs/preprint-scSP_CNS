#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

from scpviz import pAnnData as pAnnData
from scpviz import plotting as scplt
from scpviz import utils as scutils
import scanpy as sc
from anndata import AnnData

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(context='paper', style='ticks')

import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

from scipy.stats import zscore
from scipy.stats import pearsonr

import plotly.graph_objects as go
from plotly.subplots import make_subplots


# # Helper Functions

# In[2]:


def aggregate_by_metadata(df, gnames, value, sort = False):
    grouped = df.groupby(value)[gnames].mean()
    if sort:
        grouped = grouped.sort_values(by=gnames[0], ascending=False)
    return grouped

def plot_umap(xx, yy, cc=None, val=None, fig_width=8, fig_height=8, cmap=None):

    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)

    if cmap is not None :
        plt.scatter(xx, yy, s=0.5, c=val, marker='.', cmap=cmap)
    elif cc is not None :
        plt.scatter(xx, yy, s=0.5, color=cc, marker='.')

    ax.axis('equal')
    ax.set_xlim(-18, 27)
    ax.set_ylim(-18, 27)
    ax.set_xticks([])
    ax.set_yticks([])

    return fig, ax


# In[3]:


def aggregate_by_obs(adata, obs_key, genes=None, layer=None, sort=False):
    """
    Aggregate gene abundances by a categorical column in .obs.
    Computes the mean expression per group.

    Args:
        adata: AnnData object.
        obs_key: Column in adata.obs to group by (categorical or string).
        genes: Optional list of gene names (adata.var_names subset). 
               If None → uses all genes.
        layer: Optional layer name. If None → uses adata.X.
        sort: Whether to sort output by the first gene.

    Returns:
        DataFrame: groups × genes mean abundance.
    """
    # Select matrix
    X = adata.layers[layer] if layer is not None else adata.X

    # Ensure DataFrame view on selected genes
    if genes is None:
        genes = adata.var_names
    gene_idx = adata.var_names.get_indexer_for(genes)
    X_sub = pd.DataFrame(X[:, gene_idx], index=adata.obs_names, columns=genes)

    # Group by metadata
    grouped = X_sub.groupby(adata.obs[obs_key]).mean()

    if sort:
        grouped = grouped.sort_values(by=genes[0], ascending=False)

    return grouped


# In[4]:


from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import numpy as np

def cluster_genes_within_groups(df, groups, metric="correlation", method="average"):
    """
    Cluster genes *within each group block* using hierarchical clustering,
    but safely handle NaNs / invalid genes.
    """

    clustered_order = []

    for group_name, info in groups.items():
        gene_list = [g for g in info["genes"] if g in df.columns]

        if len(gene_list) <= 1:
            clustered_order.extend(gene_list)
            continue

        sub = df[gene_list].T  # genes × samples

        # ---- FILTER OUT PROBLEMATIC GENES ----
        # Drop all-NaN rows
        sub = sub.dropna(axis=0, how="all")

        # Drop zero-variance genes (correlation is undefined)
        sub = sub[sub.var(axis=1, skipna=True) > 1e-12]

        gene_list_clean = list(sub.index)

        if len(gene_list_clean) <= 1:
            # not enough valid genes to cluster
            clustered_order.extend(gene_list_clean)
            continue

        # pairwise distances
        d = pdist(sub.values, metric=metric)

        # Check for NaN or inf in distances
        if np.any(~np.isfinite(d)):
            # fallback: keep original ordering
            clustered_order.extend(gene_list_clean)
            continue

        # hierarchical clustering
        Z = linkage(d, method=method)

        leaf_order = leaves_list(Z)
        ordered = list(sub.index[leaf_order])
        clustered_order.extend(ordered)

    return clustered_order


# In[5]:


def plot_heatmap_with_group_sidebar(
    df,
    groups,
    valid_regions=None,
    fig_width=6,
    fig_height=10,
    cmap=plt.cm.coolwarm,
    vmin=-3,
    vmax=3,
    colorbar_label="Z-score",
    rotation_angle=45,
    width_ratios=[0.4, 6, 0.4],
    return_df=False,
):
    """
    df: aggregated DataFrame with rows = group/category (e.g., regions),
        columns = genes
    groups: dict defining the grouping and colors:
        {
            "GroupName": {"genes": [...], "color": "#RRGGBB"},
            ...
        }
    valid_regions: optional list defining column order for plotting
    """

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch

    # ---- 1. Build ordered gene list and sidebar colors ----
    ordered_genes = cluster_genes_within_groups(df, groups)

    sidebar_colors = []
    for group_name, info in groups.items():
        color = info["color"]
        subgroup = [g for g in ordered_genes if g in info["genes"]]
        sidebar_colors.extend([color] * len(subgroup))

    # ---- 2. Reorder df by genes & transpose (genes as rows) ----
    df = df[ordered_genes].T   # columns = genes → rows = genes

    # ---- 3. Select columns to plot (region order) ----
    if valid_regions is not None:
        df = df[[c for c in valid_regions if c in df.columns]]

    arr = df.to_numpy()

    # ---- 4. Convert group colors to RGB ----
    sidebar_rgb = np.array([mcolors.to_rgb(c) for c in sidebar_colors])
    sidebar_rgb = sidebar_rgb.reshape(len(sidebar_colors), 1, 3)

    # ---- 5. Create figure ----
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(1, 3, width_ratios=width_ratios, figure=fig)

    # Sidebar axis
    ax_side = fig.add_subplot(gs[0, 0])
    ax_side.imshow(sidebar_rgb, aspect="auto")
    ax_side.set_xticks([])
    ax_side.set_yticks([])

    # Heatmap axis
    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(arr, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)

    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels(df.columns, rotation=rotation_angle, ha="right")

    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels(df.index)

    # Colorbar
    ax_cbar = fig.add_subplot(gs[0, 2])
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label(colorbar_label)

    fig.tight_layout()

    if return_df:
        return fig, df
    else:
        return fig

def plot_group_legend(groups, figsize=(2, 1.5), fontsize=12):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    legend_patches = [
        Patch(facecolor=info["color"], label=name)
        for name, info in groups.items()
    ]

    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")

    ax.legend(
        handles=legend_patches,
        loc="center",
        frameon=False,
        fontsize=fontsize,
    )

    fig.tight_layout()
    return fig


# # Allen Brain Atlas data access

# In[6]:


download_base = Path('data/abc_atlas')
abc_cache = AbcProjectCache.from_cache_dir(download_base)

print(abc_cache.list_manifest_file_names)
print(abc_cache.current_manifest)

# load old manifest
# abc_cache.load_manifest('releases/20230630/manifest.json')
# print("old manifest loaded:", abc_cache.current_manifest)

abc_cache.load_latest_manifest()
print("manifest loaded:", abc_cache.current_manifest)


# In[7]:


## CELL METADATA

# cell = abc_cache.get_metadata_dataframe(
#     directory='WMB-10X',
#     file_name='cell_metadata'
# ).set_index('cell_label')

# # cell type classification
# cluster_details = abc_cache.get_metadata_dataframe(
#     directory='WMB-taxonomy',
#     file_name='cluster_to_cluster_annotation_membership_pivoted',
#     keep_default_na=False
# )
# cluster_details.set_index('cluster_alias', inplace=True)

# # cell type colors
# cluster_colors = abc_cache.get_metadata_dataframe(
#     directory='WMB-taxonomy',
#     file_name='cluster_to_cluster_annotation_membership_color'
# )
# cluster_colors.set_index('cluster_alias', inplace=True)

# # cell ROI
# roi = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='region_of_interest_metadata')
# roi.set_index('acronym', inplace=True)
# roi.rename(columns={'order': 'region_of_interest_order',
#                     'color_hex_triplet': 'region_of_interest_color'},
#            inplace=True)

# cell_extended = cell.join(cluster_details, on='cluster_alias')
# cell_extended = cell_extended.join(cluster_colors, on='cluster_alias')
# cell_extended = cell_extended.join(roi[['region_of_interest_order', 'region_of_interest_color']], on='region_of_interest_acronym')
# cell_extended.head(5)


# In[8]:


# Load metadata
cell_extended = abc_cache.get_metadata_dataframe(
    directory='WMB-10X',
    file_name='cell_metadata_with_cluster_annotation'
)
cell_extended.set_index('cell_label', inplace=True)


# In[9]:


cell_extended.columns


# In[10]:


gene = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene').set_index('gene_identifier')
gene.head()


# In[11]:


# rough filter to reduce number of cells we're extracting data from (ROI)
subset = cell_extended[cell_extended["feature_matrix_label"].isin(
    ["WMB-10Xv3-Isocortex-1", "WMB-10Xv3-Isocortex-2", "WMB-10Xv3-MB"]
)]

print(len(subset)) # should match 792308

roi = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='region_of_interest_metadata')
roi.set_index('acronym', inplace=True)
roi.rename(columns={'order': 'region_of_interest_order',
                    'color_hex_triplet': 'region_of_interest_color'},
           inplace=True)

print(subset['region_of_interest_acronym'].unique())
roi.loc[subset['region_of_interest_acronym'].unique()]


# In[12]:


subset.supertype.value_counts()


# In[13]:


# rough filter to reduce number of cells we're extracting data from (supertypes)

keep_supertypes = [
    "0023 L4/5 IT CTX Glut_1",
    "0027 L4/5 IT CTX Glut_5",
    "0024 L4/5 IT CTX Glut_2",
    "0028 L4/5 IT CTX Glut_6",
    "0020 L5 IT CTX Glut_3",
    "0123 L5 NP CTX Glut_2",
    "0124 L5 NP CTX Glut_3",
    "0019 L5 IT CTX Glut_2",
    "0018 L5 IT CTX Glut_1",
    "0026 L4/5 IT CTX Glut_4",
    "0089 L4 RSP-ACA Glut_1",
    "0025 L4/5 IT CTX Glut_3",
    "0883 SNc-VTA-RAmb Foxa1 Dopa_4",
    "0092 L5 ET CTX Glut_3",
    "0021 L5 IT CTX Glut_4",
    "0882 SNc-VTA-RAmb Foxa1 Dopa_3",
    "0095 L5 ET CTX Glut_6",
    "0091 L5 ET CTX Glut_2",
    "0800 SNr-VTA Pax5 Npas1 Gaba_1",
    "0090 L5 ET CTX Glut_1",
    "0880 SNc-VTA-RAmb Foxa1 Dopa_1",
    "0881 SNc-VTA-RAmb Foxa1 Dopa_2",
    "0122 L5 NP CTX Glut_1",
    "0008 L5/6 IT TPE-ENT Glut_2",
    "0125 L5 NP CTX Glut_4",
    "0093 L5 ET CTX Glut_4",
    "0887 SNc-VTA-RAmb Foxa1 Dopa_8",
    "0884 SNc-VTA-RAmb Foxa1 Dopa_5",
    "0885 SNc-VTA-RAmb Foxa1 Dopa_6",
    "0094 L5 ET CTX Glut_5",
    "0022 L5 IT CTX Glut_5",
    "0886 SNc-VTA-RAmb Foxa1 Dopa_7",
    "0126 L5 NP CTX Glut_5",
    "0009 L5/6 IT TPE-ENT Glut_3",
    "0801 SNr-VTA Pax5 Npas1 Gaba_2",
    "0011 L5/6 IT TPE-ENT Glut_5",
    "0007 L5/6 IT TPE-ENT Glut_1",
    "0010 L5/6 IT TPE-ENT Glut_4",
    "0099 L5 PPP Glut_1",
]

subset = subset[subset["supertype"].isin(keep_supertypes)]

print("Remaining cells:", len(subset)) # should match 665549
roi.loc[subset['region_of_interest_acronym'].unique()]


# In[14]:


# Make manuscript label so we can plot easily later

ctx_regions = [
    'RSP', 'ACA', 'PL-ILA-ORB', 'AUD-TEa-PERI-ECT', 'SS-GU-VISC',
    'MO-FRP', 'AI', 'VIS-PTLp', 'VIS', 'MOp'
]
snpc_regions = ['MB']

cell_extended['manuscript_region'] = 'Other'  # default if needed

cell_extended.loc[
    cell_extended['region_of_interest_acronym'].isin(ctx_regions),
    'manuscript_region'
] = 'Ctx'

cell_extended.loc[
    cell_extended['region_of_interest_acronym'].isin(snpc_regions),
    'manuscript_region'
] = 'Snpc'

cell_extended['manuscript_region'].value_counts()


# #### (optional, skip) BRIEF CHECK: double check that subset is what we want

# In[15]:


subset.columns


# In[16]:


subset['neurotransmitter'].value_counts()
subset['class'].value_counts()
# subset['subclass'].value_counts()
# subset['supertype'].value_counts()


# ## get DEPs from proteomics data

# In[17]:


## GET GENES FROM PROTEOMICS DATA - DIRECTLFQ
volcano_df = pd.read_csv('MANUSCRIPT_Fig2-CNS_DE-directlfq.csv', index_col='Genes')

all_genes = volcano_df.index.to_list()
up_list = volcano_df[volcano_df['significance']=='upregulated'].index.to_list()
down_list = volcano_df[volcano_df['significance']=='downregulated'].index.to_list()
sig_list = up_list+down_list

print("Up-regulated genes:", len(up_list))
print("Down-regulated genes:", len(down_list))
print("All genes:", len(all_genes))
all_genes.append('Otx2')


# In[18]:


## GET GENES FROM PROTEOMICS DATA
volcano_df_impute = pd.read_csv('MANUSCRIPT_Fig2-CNS_DE-impute.csv', index_col='Genes')

all_genes_impute = volcano_df_impute.index.to_list()
up_list_impute = volcano_df_impute[volcano_df_impute['significance']=='upregulated'].index.to_list()
down_list_impute = volcano_df_impute[volcano_df_impute['significance']=='downregulated'].index.to_list()
sig_list_impute = up_list_impute+down_list_impute

print("Up-regulated genes:", len(up_list_impute))
print("Down-regulated genes:", len(down_list_impute))
print("All genes:", len(all_genes_impute))
all_genes_impute.append('Otx2')


# In[19]:


## GET GENES FROM PROTEOMICS DATA - DIRECTLFQ
volcano_df_anxa = pd.read_csv('MANUSCRIPT_Fig3-Anxa_DE-directlfq-relaxed.csv', index_col='Genes')

all_genes_anxa = volcano_df_anxa.index.to_list()
up_list_anxa = volcano_df_anxa[volcano_df_anxa['significance']=='upregulated'].index.to_list()
down_list_anxa = volcano_df_anxa[volcano_df_anxa['significance']=='downregulated'].index.to_list()
sig_list_anxa = up_list_anxa+down_list_anxa

print("Up-regulated genes:", len(up_list_anxa))
print("Down-regulated genes:", len(down_list_anxa))
print("All genes:", len(all_genes_anxa))

## the true all list would be all of the detected proteins...

all_genes_anxa.append('Otx2')


# In[20]:


len(set(all_genes).union(set(all_genes_impute)))
len(set(all_genes).union(set(all_genes_anxa)))

all_genes_triple = list(set(all_genes).union(set(all_genes_anxa)))
len(all_genes_triple)


# ## (optional, skip if parquets available) get expression data
# 
# If joined_all/joined_sig.parquet is already downlaoded, skip this step and load in parquet directly

# In[ ]:


# from abc_atlas_access.abc_atlas_cache.anndata_utils import get_gene_data

# gene_names = all_genes_triple

# # this will take a decent amount of time (and will download large files)
# gene_data = get_gene_data(
#     abc_atlas_cache=abc_cache,
#     all_cells=subset,
#     all_genes=gene,
#     selected_genes=gene_names
# )

# gene_data


# In[ ]:


# import matplotlib.pyplot as plt

# fig, ax = plot_umap(subset['x'], subset['y'], cc=subset['region_of_interest_color'])
# res = ax.set_title("Dissection Region Of Interest")
# plt.show()


# In[ ]:


# # WARNING: WILL TAKE A LONG TIME (assuming you have gene_data from above)
# # all data
# cols = set(gene_data.columns)          # O(1) lookup
# valid_genes = [g for g in all_genes if g in cols]
# missing_genes = [g for g in all_genes if g not in cols]

# print(f"{len(missing_genes)} missing genes ignored.")

# filtered_all = gene_data[valid_genes]
# joined_all = cell_extended.merge(filtered_all, left_index=True, right_index=True, how="right")

# cols = joined_all.columns
# dupes = cols[cols.duplicated()].unique()

# print("Duplicate columns:", dupes)
# joined_all = joined_all.loc[:, ~cols.duplicated()]

# # just so we can load more easily next time, save joined_all
# joined_all.to_parquet("joined_all.parquet")


# # Clustering and Annotation

# ## clean data into anndata format

# In[21]:


# read in transcriptomics and proteomics data and format nicely for comparison

df_transcriptomics = pd.read_parquet("joined_all_triple.parquet", engine="pyarrow")
df_proteomics_region = pd.read_csv("proteomics_formatted.csv")


# In[22]:


df_transcriptomics


# In[23]:


# clean df_transcriptomics to keep cells that we want:
start = df_transcriptomics.columns.get_loc("manuscript_region") + 1
gene_cols = df_transcriptomics.columns[start:]
obs_cols = ["supertype", "Region","region_of_interest_acronym","class","donor_label","cluster"]
df_tx = df_transcriptomics[["supertype", "manuscript_region","region_of_interest_acronym","class","donor_label","cluster"] + gene_cols.tolist()].copy() # Subset only needed columns: supertype, manuscript_region, and all genes

# Drop unwanted stuff 
df_tx = df_tx[df_tx['region_of_interest_acronym'].isin(['MOp','MB'])]
df_tx = df_tx[~df_tx["supertype"].str.contains(r"L5 NP CTX Glut|SNr", na=False)]

# Normalize region labels to lowercase, clean naming
df_tx["manuscript_region"] = df_tx["manuscript_region"].replace({"Ctx": "ctx", "Snpc": "snpc"})
df_tx = df_tx.rename(columns={"manuscript_region": "Region"})

# Clean 0s to NaNs
gene_cols_tx = df_tx.columns.difference(obs_cols)
df_tx[gene_cols_tx] = df_tx[gene_cols_tx].apply(pd.to_numeric, errors="coerce")

df_transcriptomics_clean = df_tx
df_transcriptomics_clean


# In[24]:


df_transcriptomics_clean['Region'].value_counts()


# In[25]:


df_transcriptomics_clean['cluster'].value_counts()


# In[26]:


df_transcriptomics_clean['supertype'].value_counts()


# In[27]:


# Treat zeros as missing to standardize with transcriptomics
df_prot_region = df_proteomics_region.copy()
df_prot_region = df_prot_region.rename(columns={"Grouping": "Region"})
gene_cols = df_prot_region.columns.difference(obs_cols)
df_prot_region[gene_cols] = df_prot_region[gene_cols].apply(pd.to_numeric, errors="coerce")

df_proteomics_region_clean = df_prot_region
df_proteomics_region_clean


# ## region: cortex/snpc 

# ### clustering

# In[28]:


tx_region = df_transcriptomics_clean.copy()
prot_region = df_proteomics_region_clean.copy()

# Identify gene columns (i.e. no obs cols)
gene_cols_region = df_transcriptomics_clean.columns.difference(obs_cols)
prot_cols_region = df_proteomics_region_clean.columns.difference(obs_cols)


# In[29]:


# Subset transcriptomics and proteomics by region
tx_region = df_transcriptomics_clean.copy()
prot_region = df_proteomics_region_clean.copy()

# Identify gene columns (i.e. no obs cols)
gene_cols_region = df_transcriptomics_clean.columns.difference(obs_cols)
prot_cols_region = df_proteomics_region_clean.columns.difference(obs_cols)

# 3. Common genes between RNA and protein
# Region
genes_tx_region = gene_cols_region[tx_region[gene_cols_region].notna().any(axis=0)]
genes_prot_region = prot_cols_region[prot_region[prot_cols_region].notna().any(axis=0)]
common_genes_region = genes_tx_region.intersection(genes_prot_region)

print("Genes: ", len(genes_tx_region))
print("Prot: ", len(genes_prot_region))
print("Region:", len(common_genes_region), "genes")


# In[30]:


obs = df_transcriptomics_clean[obs_cols].copy()
obs["supertype"] = obs["supertype"].str[5:-2]

# 2. Expression matrix = all non-Grouping columns
expr = df_transcriptomics_clean.drop(columns=obs_cols)
expr = expr[common_genes_region].copy()
print(expr.shape)

# 3. Create AnnData
adata_tx = AnnData(
    X=expr.values,
    obs=obs,
    var=pd.DataFrame(index=expr.columns)
)

adata_tx_base = adata_tx.copy()
adata_tx


# In[31]:


# adata_tx.X = np.nan_to_num(adata_tx.X, nan=0.0)

# Scale each gene
sc.pp.scale(adata_tx)
sc.tl.pca(adata_tx, n_comps=30, svd_solver='arpack')
sc.pp.neighbors(adata_tx)

sc.tl.leiden(
    adata_tx, 
    flavor="leidenalg",  
    n_iterations=2,
    resolution=0.1,
)


# In[34]:


sc.tl.umap(adata_tx)
sc.pl.umap(adata_tx, color=['Syngap1', 'Homer1', 'Camk2a', 'Sncg', 'Ddc','Th'])
sc.pl.umap(adata_tx, color=["leiden", "supertype"])
sc.pl.umap(adata_tx, color=["Region","region_of_interest_acronym","class"])


# In[35]:


# create a dictionary to map cluster to annotation label
marker_genes_dict = {
    "Cortex": ['Homer1', 'Camk2a','Slc17a7'],
    "SNpC": ["Th", "Ddc","Sncg"],
}

cluster2annotation = {
    "0": "L4/5 IT CTX Glut",
    "1": "L4/5 IT CTX Glut",
    "2": "L5 IT CTX Glut",
    "3": "SNpC",
    "4": "L5 ET CTX Glut",
}
adata_tx.obs["cell type"] = adata_tx.obs["leiden"].map(cluster2annotation).astype("category")

sc.pl.dotplot(adata_tx, marker_genes_dict, "cell type", dendrogram=True)


# In[36]:


sc.tl.dendrogram(adata_tx, groupby='Region')

sc.tl.rank_genes_groups(adata_tx, groupby="Region", method="wilcoxon")
dp = sc.pl.rank_genes_groups_dotplot(adata_tx, n_genes=10, swap_axes=False, return_fig=True,)

# dp.savefig("rank_genes.svg", bbox_inches="tight")


# In[41]:


# create a dictionary to map cluster to annotation label
marker_genes_dict = {
    "Cortex": ['Nrgn', 'Celf2', 'Camk2a','Slc17a7', 'Homer1', 'Olfm1', 'Phactr1','Khdrbs3','Ptprd','Cacnb4','Gria3','Icam5'],
    "SNpC": ["Th", "Ddc", 'Ncam1', 'Gap43', 'Vat1', 'Scg2', "Chl1", 'Dpp6', 'Gaa', 'Rps8', 'Vat1l','Slc6a3',"Sncg",'Calb2','Aldh1a1'],
}

dp = sc.pl.dotplot(adata_tx, marker_genes_dict, "Region", dendrogram=True, return_fig=True)

dp.savefig("rank_genes.svg", bbox_inches="tight")


# In[33]:


adata_tx.obs.Region.value_counts()

adata_tx.obs["Region"] = adata_tx.obs["Region"].replace({
    "ctx": "Cortex",
    "snpc": "SNpc"
})


# In[35]:


cns_color={
    "Cortex": "#D19DCB",
    "SNpc": "#85BE9E"
}

sc.settings.set_figure_params(fontsize=14)

fig=sc.pl.umap(adata_tx, color=["Region"], palette=cns_color,    show=False,
    return_fig=True,)

# plt.show()
fig.savefig("../MANUSCRIPT/FigS3_UMAP.svg", bbox_inches="tight")


# In[37]:


sc.settings.set_figure_params(fontsize=14)

fig = sc.pl.umap(
    adata_tx,
    color=['Slc17a7', 'Nrgn', 'Homer1', 'Camk2a'],
    legend_fontsize='medium',
    color_map='inferno',
    show=False,
    return_fig=True,
)
fig.savefig('../MANUSCRIPT/FigS3_umap_cortex_markers.svg', bbox_inches='tight')

fig = sc.pl.umap(
    adata_tx,
    color=['Ddc', 'Ncam1', 'Th', 'Sncg'],
    legend_fontsize='medium',
    color_map='inferno',
    show=False,
    return_fig=True,
)
fig.savefig('../MANUSCRIPT/FigS3_umap_snpc_markers.svg', bbox_inches='tight')


# In[45]:


prot_list = ['Ryr3',
 'Syndig1',
 'Atp8b5',
 'Plxnb2',
 'Myh2',
 'Trip12',
 'Col15a1',
 'Cfdp1',
 'Ctsl',
 'Cck',
 'Chmp6',
 'Rrm2',
 'Chrm1',
 'Myh8',
 'Gstm3',
 'Cbl',
 'Cryaa',
 'Ptprk',
 'Grm8',
 'Fosl2',
 'Npy',
 'Aaas',
 'Kcnj3',
 'Tgfbi',
 'Kcnb1',
 'Mlip',
 'Myh4',
 'Icam5',
 'Col12a1',
 'Pzp',
 'Cep290',
 'Lpcat4',
 'Ttbk1',
 'Tiam2',
 'Abhd17b',
 'Serpina12',
 'Nup214',
 'Slitrk1',
 'Thoc5',
 'Metap1',
 'Synpo',
 'Nrn1',
 'Adgrb2',
 'Ntng2',
 'Satb2',
 'Sorcs3',
 'Dgkg',
 'Tmem132a',
 'Eps8l2',
 'Cacng3',
 'Ybx3',
 'Rgs17',
 'Mbd3',
 'Cdkn2b',
 'Crybb2',
 'Myl11',
 'Arhgef18',
 'Lrrc7',
 'Glrx5',
 'Htatsf1',
 'Hat1',
 'Septin10',
 'Aqr',
 'Rgs12',
 'Gtf2h5',
 'Exosc4',
 'Polr2h',
 'Pwp1',
 'Rnf141',
 'Acox3',
 'Zc3h18',
 'Arhgef40',
 'Ipcef1',
 'Rasl2-9',
 'Sipa1l2',
 'Mrps5',
 'Cdip1',
 'Cldn12',
 'Chm',
 'Cpne9',
 'Asphd1',
 'Them6',
 'Iffo2',
 'Oxnad1',
 'Kiaa1671']

common = sorted(set(prot_list) & set(adata_tx.var_names))
common


# In[46]:


fig = sc.pl.umap(
    adata_tx,
    color=['Icam5', 'Synpo'],
    legend_fontsize='medium',
    color_map='inferno',
    show=False,
    return_fig=True,
)


# ### abundance ROI graph

# In[38]:


## filter to only genes that were DEPs (else it's a very LONG heatmap)

# directlfq/imputation comparison
# RUN IF PARQUET ALREADY EXPORTED BEFORE
valid_genes = [g for g in sig_list if g in adata_tx.var_names]
valid_genes_impute = [g for g in sig_list_impute if g in adata_tx.var_names]


# In[33]:


# for directlfq
up_list_valid = [g for g in up_list if g in adata_tx.var_names]
down_list_valid = [g for g in down_list if g in adata_tx.var_names]

groups = {
    "up cortex":    {"genes": up_list_valid,    "color": "#EF767A"},
    "up snpc":  {"genes": down_list_valid,  "color": "#4F9FD1"},
}

# and for impute
up_list_valid_impute = [g for g in up_list_impute if g in adata_tx.var_names]
down_list_valid_impute = [g for g in down_list_impute if g in adata_tx.var_names]

groups_impute = {
    "up cortex":    {"genes": up_list_valid_impute,    "color": "#EF767A"},
    "up snpc":  {"genes": down_list_valid_impute,  "color": "#4F9FD1"},
}

# other aesthetics 
cmap = plt.cm.bwr


# In[35]:


agg_ROI = aggregate_by_obs(
    adata_tx,
    obs_key='region_of_interest_acronym',
    genes=valid_genes
)

roi_order = roi.index
valid_order = roi_order.intersection(agg_ROI.index)
agg_ROI = agg_ROI.loc[valid_order]
agg_ROI = agg_ROI.replace("nan", np.nan).apply(pd.to_numeric, errors='coerce')

agg_ROI_impute = aggregate_by_obs(
    adata_tx,
    obs_key='region_of_interest_acronym',
    genes=valid_genes_impute
)

roi_order = roi.index
valid_order_impute = roi_order.intersection(agg_ROI_impute.index)
agg_ROI_impute = agg_ROI_impute.loc[valid_order_impute]
agg_ROI_impute = agg_ROI_impute.replace("nan", np.nan).apply(pd.to_numeric, errors='coerce')


# In[36]:


# roi plot
agg_ROI_z = agg_ROI.apply(lambda col: zscore(col, nan_policy='omit'), axis=0)
agg = agg_ROI_z
valid_regions = agg_ROI_z.index.to_list()
fig_size=(45,4) # height, width

fig, df_directlfq = plot_heatmap_with_group_sidebar(
    agg,
    groups,
    valid_regions=valid_regions,
    vmin=-3,
    vmax=3,
    colorbar_label="Z-score",
    fig_height=fig_size[0], 
    fig_width=fig_size[1], 
    cmap=cmap,
    width_ratios=[0.2, 6, 0.2],
    return_df=True
)
plt.title("ROI")
# plt.show()
plt.savefig('DEP_ROI_directlfq.svg')


# In[37]:


# roi plot - impute
agg_ROI_z_impute = agg_ROI_impute.apply(lambda col: zscore(col, nan_policy='omit'), axis=0)
agg_impute = agg_ROI_z_impute
valid_regions_impute = agg_ROI_z_impute.index.to_list()
fig_size=(45,4) # height, width

fig, df_impute = plot_heatmap_with_group_sidebar(
    agg_impute,
    groups_impute,
    valid_regions=valid_regions_impute,
    vmin=-3,
    vmax=3,
    colorbar_label="Z-score",
    fig_height=fig_size[0], 
    fig_width=fig_size[1], 
    cmap=cmap,
    width_ratios=[0.2, 6, 0.2],
    return_df=True
)
plt.title("ROI_impute")
# plt.show()
plt.savefig('DEP_ROI_impute.svg')


# ### DE

# In[ ]:


if not hasattr(np, "Inf"):
    np.Inf = np.inf

values=[{'region':'Cortex'},{'region':'SNpC'}]

fig, ax = plt.subplots(figsize = (3.5,3.5))
ax, df_snpc = scplt.plot_volcano_adata(ax, adata_tx, values=values, return_df=True, log2fc=0.5)

df_snpc.to_csv("transcriptomics-snpc_DE_volcano.csv")


# In[49]:


fig, ax = plt.subplots(figsize = (3.5,3.5))
ax, df_mark = scplt.plot_volcano_adata(ax, adata_snpc, values=values, no_marks=True, return_df=True,log2fc=0.5)
scplt.mark_volcano(ax, df_mark, label=['Anxa1','Aldh1a7','Aldh1a1'], label_color='red')
scplt.mark_volcano(ax, df_mark, label=['Prph'], label_color='blue')


# # Omics comparisons

# ## region: cortex/snpc

# In[16]:


adata_tx = adata_tx_base.copy()

adata_tx.obs.region_of_interest_acronym.value_counts()

# only comparing MOp (cortex) and MB (SNpC)


# In[17]:


def proteomics_direction(gene, up_list, down_list):
    if gene in up_list:
        return "proteomics_up"
    elif gene in down_list:
        return "proteomics_down"
    else:
        return "proteomics_NS"

def combined_consistency(row):
    t = row.transcriptomics_direction
    p = row.proteomics_direction

    # --- Strong consistent calls ---
    if t == "transcriptomics_up" and p == "proteomics_up":
        return "consistent_up"
    if t == "transcriptomics_down" and p == "proteomics_down":
        return "consistent_down"

    # --- Strong opposite calls ---
    if t == "transcriptomics_up" and p == "proteomics_down":
        return "opposite"
    if t == "transcriptomics_down" and p == "proteomics_up":
        return "opposite"

    # --- Explicit no-signal states ---
    if t == "transcriptomics_NS" and p != "proteomics_NS":
        return "transcriptomics_NS"
    if p == "proteomics_NS" and t != "transcriptomics_NS":
        return "proteomics_NS"

    # Both NS
    if t == "transcriptomics_NS" and p == "proteomics_NS":
        return "both_NS"

    # --- Residual category: defined directions but mismatched ---
    return "inconsistent"


# ### directlfq

# In[26]:


agg_ROI = aggregate_by_obs(
    adata_tx,
    obs_key='region_of_interest_acronym',
    genes=adata_tx.var_names
)

roi_order = roi.index
valid_order = roi_order.intersection(agg_ROI.index)

agg_ROI = agg_ROI.loc[valid_order]
agg_ROI = agg_ROI.replace("nan", np.nan).apply(pd.to_numeric, errors='coerce')


# In[27]:


agg_linear = 2 ** agg_ROI      # undo log2

cortex_regions = ["MOp"]
mb_region = "MB"

mb_linear = agg_linear.loc[mb_region]

# region-wise fold-change vs MB
fc_linear = agg_linear.loc[cortex_regions].div(mb_linear)
fc_log2 = np.log2(fc_linear)


# In[34]:


THRESH = 1.0 

# Zero out small changes before sign calling
fc_sign = np.where(np.abs(fc_log2) >= THRESH, np.sign(fc_log2), 0)

consistency = pd.Series(fc_sign.sum(axis=0), index=fc_log2.columns)

direction_label = consistency.map({
    1:  "transcriptomics_up",
    -1: "transcriptomics_down",
}).fillna("transcriptomics_NS")

fold_change_df = fc_log2.T.copy()
fold_change_df.columns = [f"transcriptomics_log2fc" for r in cortex_regions]
fold_change_df["proteomics_log2fc"] = volcano_df["log2fc"].reindex(fold_change_df.index)

fold_change_df["transcriptomics_direction"] = direction_label
fold_change_df["proteomics_direction"] = [
    proteomics_direction(g, up_list_valid, down_list_valid) for g in fold_change_df.index
]
fold_change_df["combined_consistency"] = fold_change_df.apply(combined_consistency, axis=1)

# add to fold_change_df
fold_change_df.to_csv('DEG_consistency_directlfq_MOPonly.csv')
fold_change_df


# In[53]:


x = fold_change_df["proteomics_log2fc"]
y = fold_change_df["transcriptomics_log2fc"]

cats = fold_change_df["combined_consistency"]

# keep only genes with directional proteomics call
prot_dir = fold_change_df["proteomics_direction"]
mask = (
    ~(x.isna() | y.isna()) &
    prot_dir.isin(["proteomics_up", "proteomics_down"])
)

xv = x[mask]
yv = y[mask]

print(f"Using {mask.sum()} genes with proteomics_up/down for correlation")

# Pearson correlation on directional-only subset
r, p = pearsonr(xv, yv)
print(f"Pearson r = {r:.3f}, p = {p:.3e}")

# Fit line on the same subset
m, b = np.polyfit(xv, yv, 1)
x_line = np.linspace(xv.min(), xv.max(), 100)
y_line = m * x_line + b


# In[54]:


# CORE
color_map_core = {
    "consistent_up":   "#86BCD8",
    "consistent_down": "#22316B",
    "opposite":        "#B4353B",
    "inconsistent":    "#EFEFEF",
}
GREY = "#D3D3D3"

alpha_map = {
    "consistent_up":   0.9,
    "consistent_down": 0.9,
    "opposite":        0.9,
    "inconsistent":    0.15,
    # everything else (NS categories)
    "default":         0.15,   # faded grey points
}

base_colors = cats.map(lambda c: color_map_core.get(c, GREY))
alphas = cats.map(lambda c: alpha_map.get(c, alpha_map["default"])).values
rgba = np.array([mcolors.to_rgba(c) for c in base_colors])
rgba[:, 3] = alphas

edge_rgba = np.zeros_like(rgba)
edge_rgba[:, :3] = 0      # black RGB
edge_rgba[:, 3] = alphas

fig, ax = plt.subplots(figsize=(3, 3))

ax.scatter(x, y, c=rgba, s=10)

ax.axhline(0, color="black", linestyle="--", linewidth=1)
ax.axvline(0, color="black", linestyle="--", linewidth=1)
ax.plot(x_line, y_line, color="black", linewidth=1)

handles = [
    plt.Line2D([0],[0], marker="o", markersize=9, color="w",
               markerfacecolor=color_map_core[k], label=k)
    for k in color_map_core
]
ax.legend(handles=handles, title="Category", frameon=False)
scplt.shift_legend(ax)

# ---- Linear fit line + r annotation ----
ax.plot(x_line, y_line, color="black", linewidth=1)

# place the r text just above the middle of the line, in data coords
xmid = 0.5 * (xv.min() + xv.max())
ymid = m * xmid + b

# small vertical offset so it sits "above" the line
ymin, ymax = ax.get_ylim()
offset = 0.02 * (ymax - ymin)

ax.text(
    xmid,
    ymid + offset,
    f"r = {r:.2f}",
    ha="center",
    va="bottom",
    fontsize=8,
    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="black", alpha=0.8),
)

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# ax.set_title("Proteomics vs Transcriptomics Fold Change (directLFQ)")
ax.set_xlabel("Proteomics $\log_2$FC", fontsize=13, labelpad=4)
ax.set_ylabel("Transcriptomics $\log_2$FC", fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

ax.set_ylim(-10,10)
ax.set_xlim(-10,10)

ax.xaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))

plt.savefig('../MANUSCRIPT/MANUSCRIPT_Fig3_OmicsComparison-directLFQ.svg', bbox_inches="tight",)


# In[124]:


# FULL
color_map_full = {
    "consistent_up":        "#B2182B",  # dark red
    "consistent_down":      "#2166AC",  # dark blue
    "opposite":             "#762A83",  # dark purple

    "transcriptomics_NS":   "#FFD580",
    "proteomics_NS":        "#B0E57C",
    "both_NS":              "#DDDDDD",

    "inconsistent":         "darkgray",
}

GREY = "#D3D3D3"

# Alpha per category
alpha_map_full = {
    "consistent_up":        0.3,
    "consistent_down":      0.3,
    "opposite":             0.3,
    "transcriptomics_NS":   0.90,
    "proteomics_NS":        0.90,
    "both_NS":              0.90,
    "inconsistent":         0.15,
    "default":              0.15,
}

cats = fold_change_df["combined_consistency"]

# Base face colors per point
base_colors = cats.map(lambda c: color_map_full.get(c, GREY)).values

# Alpha per point
alphas = cats.map(lambda c: alpha_map_full.get(c, alpha_map_full["default"])).values

# Convert hex → RGBA
rgba = np.array([mcolors.to_rgba(c) for c in base_colors])
rgba[:, 3] = alphas

# Edge RGBA (black with matching alpha)
edge_rgba = np.zeros_like(rgba)
edge_rgba[:, :3] = 0       # black RGB
edge_rgba[:, 3] = alphas   # alpha per point

fig, ax = plt.subplots(figsize=(4, 4))

# Scatter with full RGBA control
ax.scatter(
    x,
    y,
    facecolors=rgba,
    edgecolors=edge_rgba,
    linewidths=0.3,
    s=12,
)

# Axes overlays
ax.axhline(0, color="black", linestyle="--", linewidth=1)
ax.axvline(0, color="black", linestyle="--", linewidth=1)
ax.plot(x_line, y_line, color="black", linewidth=1)

ax.set_xlabel("Proteomics log2FC")
ax.set_ylabel("Transcriptomics log2FC")
ax.set_title("directlfq: All Categories")

# Legend
legend_order = list(color_map_full.keys())
present = cats.dropna().unique()
legend_filtered = [k for k in legend_order if k in present]

handles = [
    plt.Line2D(
        [0],[0],
        marker="o", markersize=9, color="w",
        markerfacecolor=color_map_full[k],
        markeredgecolor="k",
        label=k
    )
    for k in legend_filtered
]

ax.legend(handles=handles, title="Category", frameon=False)
scplt.shift_legend(ax)

# ---- Linear fit line + r annotation ----
ax.plot(x_line, y_line, color="black", linewidth=1)

# place the r text just above the middle of the line, in data coords
xmid = 0.5 * (xv.min() + xv.max())
ymid = m * xmid + b

# small vertical offset so it sits "above" the line
ymin, ymax = ax.get_ylim()
offset = 0.02 * (ymax - ymin)

ax.text(
    xmid,
    ymid + offset,
    f"r = {r:.2f}",
    ha="center",
    va="bottom",
    fontsize=8,
    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="black", alpha=0.8),
)

ax.set_xlim(-8,8)


# ### impute

# In[45]:


THRESH = 1.0 

# Zero out small changes before sign calling
fc_sign = np.where(np.abs(fc_log2) >= THRESH, np.sign(fc_log2), 0)

consistency = pd.Series(fc_sign.sum(axis=0), index=fc_log2.columns)

direction_label = consistency.map({
    1:  "transcriptomics_up",
    -1: "transcriptomics_down",
}).fillna("transcriptomics_NS")

fold_change_df_impute = fc_log2.T.copy()
fold_change_df_impute.columns = [f"transcriptomics_log2fc" for r in cortex_regions]
fold_change_df_impute["proteomics_log2fc"] = volcano_df_impute["log2fc"].reindex(fold_change_df_impute.index)

fold_change_df_impute["transcriptomics_direction"] = direction_label
fold_change_df_impute["proteomics_direction"] = [
    proteomics_direction(g, up_list_valid_impute, down_list_valid_impute) for g in fold_change_df_impute.index
]
fold_change_df_impute["combined_consistency"] = fold_change_df_impute.apply(combined_consistency, axis=1)

# add to fold_change_df
fold_change_df_impute.to_csv('DEG_consistency_impute_MOPonly.csv')
fold_change_df_impute


# In[46]:


x = fold_change_df_impute["proteomics_log2fc"]
y = fold_change_df_impute["transcriptomics_log2fc"]

cats = fold_change_df_impute["combined_consistency"]

# keep only genes with directional proteomics call
prot_dir = fold_change_df_impute["proteomics_direction"]
mask = (
    ~(x.isna() | y.isna()) &
    prot_dir.isin(["proteomics_up", "proteomics_down"])
)

xv = x[mask]
yv = y[mask]

print(f"Using {mask.sum()} genes with proteomics_up/down for correlation")

# Pearson correlation on directional-only subset
r, p = pearsonr(xv, yv)
print(f"Pearson r = {r:.3f}, p = {p:.3e}")

# Fit line on the same subset
m, b = np.polyfit(xv, yv, 1)
x_line = np.linspace(xv.min(), xv.max(), 100)
y_line = m * x_line + b


# In[52]:


# CORE
color_map_core = {
    "consistent_up":   "#86BCD8",
    "consistent_down": "#22316B",
    "opposite":        "#B4353B",
    "inconsistent":    "#EFEFEF",
}
GREY = "#D3D3D3"

alpha_map = {
    "consistent_up":   0.9,
    "consistent_down": 0.9,
    "opposite":        0.9,
    "inconsistent":    0.15,
    # everything else (NS categories)
    "default":         0.15,   # faded grey points
}

base_colors = cats.map(lambda c: color_map_core.get(c, GREY))
alphas = cats.map(lambda c: alpha_map.get(c, alpha_map["default"])).values
rgba = np.array([mcolors.to_rgba(c) for c in base_colors])
rgba[:, 3] = alphas

edge_rgba = np.zeros_like(rgba)
edge_rgba[:, :3] = 0      # black RGB
edge_rgba[:, 3] = alphas

fig, ax = plt.subplots(figsize=(3, 3))

ax.scatter(x, y, c=rgba, s=10)

ax.axhline(0, color="black", linestyle="--", linewidth=1)
ax.axvline(0, color="black", linestyle="--", linewidth=1)
ax.plot(x_line, y_line, color="black", linewidth=1)

handles = [
    plt.Line2D([0],[0], marker="o", markersize=9, color="w",
               markerfacecolor=color_map_core[k], label=k)
    for k in color_map_core
]
ax.legend(handles=handles, title="Category", frameon=False)
scplt.shift_legend(ax)

# ---- Linear fit line + r annotation ----
ax.plot(x_line, y_line, color="black", linewidth=1)

# place the r text just above the middle of the line, in data coords
xmid = 0.5 * (xv.min() + xv.max())
ymid = m * xmid + b

# small vertical offset so it sits "above" the line
ymin, ymax = ax.get_ylim()
offset = 0.02 * (ymax - ymin)

ax.text(
    xmid,
    ymid + offset,
    f"r = {r:.2f}",
    ha="center",
    va="bottom",
    fontsize=8,
    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="black", alpha=0.8),
)

ax.set_xlabel("Proteomics $\log_2$FC", fontsize=13, labelpad=4)
ax.set_ylabel("Transcriptomics $\log_2$FC", fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

ax.set_ylim(-10,10)
ax.set_xlim(-10,10)

ax.xaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))

ax.set_xlim(-8,8)

plt.savefig('../MANUSCRIPT/MANUSCRIPT_Fig3_OmicsComparison-imputeNorm.svg', bbox_inches="tight",)


# ### comparison of methods

# In[ ]:


categories = [
    "consistent_up",
    "consistent_down",
    "opposite",
    "NS",
]

collapse_map = {
    "transcriptomics_NS": "NS",
    "proteomics_NS":      "NS",
    "both_NS":            "NS",
    "inconsistent":       "NS",
    "consistent_up":      "consistent_up",
    "consistent_down":    "consistent_down",
    "opposite":           "opposite",
}
collapsed_impute = fold_change_df_impute["combined_consistency"].map(collapse_map)
collapsed_directlfq = fold_change_df["combined_consistency"].map(collapse_map)

counts_impute = collapsed_impute.value_counts().reindex(categories, fill_value=0)
counts_directlfq = collapsed_directlfq.value_counts().reindex(categories, fill_value=0)

# Build combined table
plot_df = pd.DataFrame({
    "category": categories,
    "imputed": counts_impute.values,
    "directlfq": counts_directlfq.values,
})
plot_df


# Convert to long format for convenience
plot_long = plot_df.melt(id_vars="category", var_name="dataset", value_name="count")

color_map = {
    "consistent_up":   "#86BCD8",
    "consistent_down": "#22316B",
    "opposite":        "#B4353B",
    "NS":              "#EFEFEF",
}

# Prepare figure
fig, ax = plt.subplots(figsize=(6, 4))

datasets = ["directlfq", "imputed"]
x_positions = [0, 1]

categories = ["consistent_up", "consistent_down", "opposite","NS"]

# Stacked bars for each dataset
for i, ds in enumerate(datasets):
    bottom = 0
    for cat in categories:  
        value = plot_df.loc[plot_df["category"] == cat, ds].iloc[0]
        ax.bar(
            x_positions[i],
            value,
            bottom=bottom,
            label=cat if i == 0 else None,   # only label first stack
            color=color_map[cat],
        )
        bottom += value

# Formatting
ax.set_xticks(x_positions)
ax.set_xticklabels(["directLFQ", "Imputed"], fontsize=12)
ax.set_ylabel("Protein count", fontsize=12)
ax.set_title("Consistency category comparison: directLFQ vs Imputed (All Proteins)", fontsize=13)

# Legend
ax.legend(
    title="Combined Consistency",
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

plt.tight_layout()
plt.show()


# Comparing DE **proteins**
# 
# Out of those: 
# - how many transcriptomics agree (consistent_up/down)
# - how many disagree (opposite)
# - how many transcriptomics are not significant (NS)
# 
# to answer the biologically meaningful question:
# *"When proteomics says a protein is significantly up/down, how does transcriptomics behave?"*

# In[72]:


# 1. Filter to proteomics-significant proteins (up/down only)
sig_direct = fold_change_df[
    fold_change_df["proteomics_direction"].isin(["proteomics_up", "proteomics_down"])
]

sig_impute = fold_change_df_impute[
    fold_change_df_impute["proteomics_direction"].isin(["proteomics_up", "proteomics_down"])
]

# 2. Map to collapsed categories
collapse_map = {
    "transcriptomics_NS": "NS",
    "proteomics_NS":      "NS",
    "both_NS":            "NS",
    "inconsistent":       "NS",
    "consistent_up":      "consistent_up",
    "consistent_down":    "consistent_down",
    "opposite":           "opposite",
}

categories = ["consistent_up", "consistent_down", "opposite", "NS"]

collapsed_direct = sig_direct["combined_consistency"].map(collapse_map)
collapsed_impute = sig_impute["combined_consistency"].map(collapse_map)

# 3. Count categories ONLY within the significant-proteomics subset
counts_directlfq = collapsed_direct.value_counts().reindex(categories, fill_value=0)
counts_impute    = collapsed_impute.value_counts().reindex(categories, fill_value=0)

# 4. Build plotting DataFrame
plot_df = pd.DataFrame({
    "category": categories,
    "directlfq": counts_directlfq.values,
    "imputed": counts_impute.values,
})

# 5. Plot
color_map = {
    "consistent_up":   "#86BCD8",
    "consistent_down": "#22316B",
    "opposite":        "#B4353B",
    "NS":              "#EFEFEF",
}

fig, ax = plt.subplots(figsize=(4, 6.2))

datasets = ["directlfq", "imputed"]
x_positions = [0, 1]

categories = ["consistent_up", "consistent_down", "opposite","NS"]

# Stacked bars
for i, ds in enumerate(datasets):
    bottom = 0
    for cat in categories:
        value = plot_df.loc[plot_df["category"] == cat, ds].iloc[0]
        ax.bar(
            x_positions[i],
            value,
            bottom=bottom,
            label=cat if i == 0 else None,
            color=color_map[cat],
        )
        bottom += value

# Axis labels
ax.set_xticks(x_positions)
ax.set_xticklabels(["directLFQ", "Imputed\n+ Median"], fontsize=12)
ax.set_ylabel("Proteomics DEP (number)", fontsize=12)
# ax.set_title("Agreement with Transcriptomics for Proteomics-DE Proteins", fontsize=13)

ax.legend(
    title="Consistency",
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)


plt.tight_layout()
# plt.show()
plt.savefig("../MANUSCRIPT/MANUSCRIPT_Fig3_Omics_Counts.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


collapse_map2 = {
    "consistent_up":   "matching",
    "consistent_down": "matching",
    "opposite":        "opposite",
    # everything else → inconsistent
    # "inconsistent":         "inconsistent",
    # "transcriptomics_NS":   "inconsistent",
    # "proteomics_NS":        "inconsistent",
    # "both_NS":              "inconsistent",
}

collapsed_direct = fold_change_df["combined_consistency"].map(collapse_map2)
collapsed_impute = fold_change_df_impute["combined_consistency"].map(collapse_map2)

categories2 = ["matching", "opposite"]

df2 = pd.DataFrame({
    "category": categories2,
    "directlfq": collapsed_direct.value_counts().reindex(categories2, fill_value=0).values,
    "imputed": collapsed_impute.value_counts().reindex(categories2, fill_value=0).values,
})
df2


df2_pct = df2.copy()
df2_pct["directlfq"] = df2["directlfq"] / df2["directlfq"].sum() * 100
df2_pct["imputed"]   = df2["imputed"]   / df2["imputed"].sum()   * 100

color_map2 = {
    "matching": "#42C478",   # green
    "opposite": "#C64D4A",   # red
}

fig, ax = plt.subplots(figsize=(6, 4))

datasets = ["directlfq", "imputed"]
x_positions = [0, 1]
categories2 = ["matching", "opposite"]

for i, ds in enumerate(datasets):
    bottom = 0  # now in percent
    for cat in categories2:
        value = df2_pct.loc[df2_pct["category"] == cat, ds].iloc[0]

        ax.bar(
            x_positions[i],
            value,
            bottom=bottom,
            color=color_map2[cat],
            label=cat if i == 0 else None,
        )

        if value > 0:
            ax.text(
                x_positions[i],
                bottom + value / 2,
                f"{value:.1f}%",
                ha='center',
                va='center',
                fontsize=10,
                color='k'
            )

        bottom += value

# Formatting
ax.set_xticks(x_positions)
ax.set_xticklabels(["directLFQ", "Imputed"], fontsize=12)
ax.set_ylabel("Percentage (%)", fontsize=12)
ax.set_ylim(0, 100)
ax.set_title("Matching vs Opposite — 100% Stacked Comparison", fontsize=13)

# Legend
ax.legend(
    title="Category",
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

plt.tight_layout()
plt.show()


# In[190]:


# Filter to significant proteomics DEPs
sig_direct = fold_change_df[
    fold_change_df["proteomics_direction"].isin(["proteomics_up", "proteomics_down"])
]

sig_impute = fold_change_df_impute[
    fold_change_df_impute["proteomics_direction"].isin(["proteomics_up", "proteomics_down"])
]

# Collapse mapping
collapse_map2 = {
    "consistent_up":   "matching",
    "consistent_down": "matching",
    "opposite":        "opposite",

    # Everything else → NS (not matching and not opposite)
    "inconsistent":       "NS",
    "transcriptomics_NS": "NS",
    "proteomics_NS":      "NS",
    "both_NS":            "NS",
}

# Apply collapsing
collapsed_direct = sig_direct["combined_consistency"].map(collapse_map2)
collapsed_impute = sig_impute["combined_consistency"].map(collapse_map2)

categories2 = ["matching", "opposite", "NS"]

# Count occurrences
counts_direct = collapsed_direct.value_counts().reindex(categories2, fill_value=0)
counts_impute = collapsed_impute.value_counts().reindex(categories2, fill_value=0)

df2 = pd.DataFrame({
    "category": categories2,
    "directlfq": counts_direct.values,
    "imputed":   counts_impute.values,
})

df2_pct = df2.copy()
df2_pct["directlfq"] = df2["directlfq"] / df2["directlfq"].sum() * 100
df2_pct["imputed"]   = df2["imputed"]   / df2["imputed"].sum()   * 100

color_map2 = {
    "matching": "#42C478",   # green
    "opposite": "#C64D4A",   # red
    "NS":       "#D3D3D3",   # light gray
}

fig, ax = plt.subplots(figsize=(6, 4))

datasets = ["directlfQ", "imputed"]
x_positions = [0, 1]

for i, ds in enumerate(["directlfq", "imputed"]):
    bottom = 0
    for cat in categories2:
        value = df2_pct.loc[df2_pct["category"] == cat, ds].iloc[0]

        ax.bar(
            x_positions[i],
            value,
            bottom=bottom,
            color=color_map2[cat],
            label=cat if i == 0 else None
        )

        # Add percentage label inside bar
        if value > 0:
            ax.text(
                x_positions[i],
                bottom + value / 2,
                f"{value:.1f}%",
                ha='center',
                va='center',
                fontsize=10,
                color='black'
            )

        bottom += value

# Formatting
ax.set_xticks(x_positions)
ax.set_xticklabels(["directLFQ", "Imputed"], fontsize=12)
ax.set_ylabel("Percentage (%)", fontsize=12)
ax.set_ylim(0, 100)
ax.set_title("Matching / Opposite / NS for DE Proteins — 100% Stacked", fontsize=13)

ax.legend(
    title="Category",
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

plt.tight_layout()
plt.show()


# In[73]:


sig_direct = fold_change_df[
    fold_change_df["proteomics_direction"].isin(["proteomics_up", "proteomics_down"])
]
sig_impute = fold_change_df_impute[
    fold_change_df_impute["proteomics_direction"].isin(["proteomics_up", "proteomics_down"])
]

collapse_map = {
    "transcriptomics_NS": "NS",
    "proteomics_NS":      "NS",
    "both_NS":            "NS",
    "inconsistent":       "NS",
    "consistent_up":      "consistent_up",
    "consistent_down":    "consistent_down",
    "opposite":           "opposite",
}
categories = ["consistent_up", "consistent_down", "opposite", "NS"]

collapsed_direct = sig_direct["combined_consistency"].map(collapse_map)
collapsed_impute = sig_impute["combined_consistency"].map(collapse_map)

counts_directlfq = collapsed_direct.value_counts().reindex(categories, fill_value=0)
counts_impute    = collapsed_impute.value_counts().reindex(categories, fill_value=0)

plot_df = pd.DataFrame({
    "category": categories,
    "directlfq": counts_directlfq.values,
    "imputed": counts_impute.values,
})

# --- convert to 100% stacked (column-wise percentages) ---
pct_df = plot_df.set_index("category")[["directlfq", "imputed"]]
pct_df = pct_df.div(pct_df.sum(axis=0), axis=1).fillna(0)  # each column sums to 1

# --- plot ---
color_map = {
    "consistent_up":   "#86BCD8",
    "consistent_down": "#22316B",
    "opposite":        "#B4353B",
    "NS":              "#EFEFEF",
}

fig, ax = plt.subplots(figsize=(4, 6.2))

datasets = ["directlfq", "imputed"]
x_positions = [0, 1]

for i, ds in enumerate(datasets):
    bottom = 0
    for cat in categories:
        value = pct_df.loc[cat, ds]
        ax.bar(
            x_positions[i],
            value,
            bottom=bottom,
            label=cat if i == 0 else None,
            color=color_map[cat],
        )
        bottom += value

ax.set_xticks(x_positions)
ax.set_xticklabels(["directLFQ", "Imputed \n+ Median"], fontsize=12)

from matplotlib.ticker import PercentFormatter

ax.set_ylim(0, 1)
ax.yaxis.set_major_formatter(PercentFormatter(1))
ax.set_ylabel("Proteomics DEP (%)", fontsize=12)
# ax.set_title("Agreement with Transcriptomics for Proteomics-DE Proteins", fontsize=13)
ax.legend(title="Consistency", bbox_to_anchor=(1.05, 1), loc="upper left")


plt.tight_layout()
plt.savefig("../MANUSCRIPT/MANUSCRIPT_Fig3_Omics_Stacked.svg", bbox_inches="tight", pad_inches=0.1)
# plt.show()


# ### sankey

# #### directlfq

# In[182]:


mask = fold_change_df["proteomics_direction"].isin(["proteomics_up", "proteomics_down"])
df_sankey = fold_change_df[mask].copy()

agreement_collapsed = df_sankey["combined_consistency"].map({
    "consistent_up":   "matching",
    "consistent_down": "matching",
    "opposite":        "opposite",
    "inconsistent":         "NS_inconsistent",
    "transcriptomics_NS":   "NS_inconsistent",
    "proteomics_NS":        "NS_inconsistent",
    "both_NS":              "NS_inconsistent",
})
df_sankey["agreement"] = agreement_collapsed

flow = (
    df_sankey.groupby([
        "proteomics_direction",
        "agreement",
        "transcriptomics_direction"
    ]).size().reset_index(name="count")
)

nodes = [
    # Layer 1
    "proteomics_up",
    "proteomics_down",

    # Layer 2
    "matching",
    "opposite",
    "NS_inconsistent",

    # Layer 3
    "transcriptomics_up",
    "transcriptomics_down",
    "transcriptomics_NS",
]

node_index = {name: i for i, name in enumerate(nodes)}

node_colors = [
    "#B2182B",  # proteomics_up
    "#2166AC",  # proteomics_down

    "#42C478",  # matching
    "#762A83",  # opposite
    "#D3D3D3",  # NS_inconsistent

    "#B2182B",  # transcriptomics_up
    "#2166AC",  # transcriptomics_down
    "#D3D3D3",  # transcriptomics_NS
]

LINK_GREY = "rgba(150,150,150,0.35)"

sources = []
targets = []
values = []
link_colors = []

for _, row in flow.iterrows():
    src = node_index[row["proteomics_direction"]]
    mid = node_index[row["agreement"]]
    tgt = node_index[row["transcriptomics_direction"]]
    count = row["count"]

    # proteomics → agreement
    sources.append(src)
    targets.append(mid)
    values.append(count)
    link_colors.append(LINK_GREY)

    # agreement → transcriptomics
    sources.append(mid)
    targets.append(tgt)
    values.append(count)
    link_colors.append(LINK_GREY)

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=25,
        thickness=15,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=node_colors,
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=link_colors,
    ),
)])

fig.update_layout(
    title_text="Directional Consistency Sankey (DEPs only) — directLFQ",
    font=dict(size=11),
)


fig.show()


# In[164]:


# Proteomics + Transcriptomics directional colors
COLOR_UP   = "#B2182B"   # dark red
COLOR_DOWN = "#2166AC"   # dark blue
COLOR_NS   = "#D3D3D3"   # light grey

# Agreement colors
COLOR_MATCHING = "#42C478"   # green
COLOR_OPPOSITE = "#762A83"   # purple
COLOR_NSI      = "#D3D3D3"   # grey for NS_inconsistent

# Link color (uniform)
LINK_GREY = "rgba(150,150,150,0.35)"

agreement_collapsed_direct = fold_change_df["combined_consistency"].map({
    "consistent_up":   "matching",
    "consistent_down": "matching",
    "opposite":        "opposite",
    "inconsistent":         "NS_inconsistent",
    "transcriptomics_NS":   "NS_inconsistent",
    "proteomics_NS":        "NS_inconsistent",
    "both_NS":              "NS_inconsistent",
})

df_sankey = fold_change_df.copy()
df_sankey["agreement"] = agreement_collapsed_direct

flow = (
    df_sankey.groupby([
        "proteomics_direction",
        "agreement",
        "transcriptomics_direction"
    ]).size().reset_index(name="count")
)

nodes = [
    # Layer 1: Proteomics
    "proteomics_up",
    "proteomics_down",
    "proteomics_NS",

    # Layer 2: Agreement
    "matching",
    "opposite",
    "NS_inconsistent",

    # Layer 3: Transcriptomics
    "transcriptomics_up",
    "transcriptomics_down",
    "transcriptomics_NS",
]

node_index = {name: i for i, name in enumerate(nodes)}


node_colors = [
    COLOR_UP,        # proteomics_up
    COLOR_DOWN,      # proteomics_down
    COLOR_NS,        # proteomics_NS

    COLOR_MATCHING,  # matching
    COLOR_OPPOSITE,  # opposite
    COLOR_NSI,       # NS_inconsistent

    COLOR_UP,        # transcriptomics_up
    COLOR_DOWN,      # transcriptomics_down
    COLOR_NS,        # transcriptomics_NS
]

sources = []
targets = []
values  = []
link_colors = []

for _, row in flow.iterrows():

    src = node_index[row["proteomics_direction"]]
    mid = node_index[row["agreement"]]
    tgt = node_index[row["transcriptomics_direction"]]
    count = row["count"]

    # proteomics → agreement
    sources.append(src)
    targets.append(mid)
    values.append(count)
    link_colors.append(LINK_GREY)

    # agreement → transcriptomics
    sources.append(mid)
    targets.append(tgt)
    values.append(count)
    link_colors.append(LINK_GREY)


fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=25,
        thickness=15,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=node_colors,
        hovertemplate="%{label}<extra></extra>",
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=link_colors,
        hovertemplate="Count: %{value}<extra></extra>",
    ),
)])

fig.update_layout(
    title_text="Directional Consistency Sankey (all proteins) — directLFQ",
    font=dict(size=14),
)

fig.show()



# #### impute

# In[171]:


mask = fold_change_df_impute["proteomics_direction"].isin(["proteomics_up", "proteomics_down"])
df_sankey = fold_change_df_impute[mask].copy()

agreement_collapsed = df_sankey["combined_consistency"].map({
    "consistent_up":   "matching",
    "consistent_down": "matching",
    "opposite":        "opposite",
    "inconsistent":         "NS_inconsistent",
    "transcriptomics_NS":   "NS_inconsistent",
    "proteomics_NS":        "NS_inconsistent",
    "both_NS":              "NS_inconsistent",
})
df_sankey["agreement"] = agreement_collapsed

flow = (
    df_sankey.groupby([
        "proteomics_direction",
        "agreement",
        "transcriptomics_direction"
    ]).size().reset_index(name="count")
)

nodes = [
    # Layer 1
    "proteomics_up",
    "proteomics_down",

    # Layer 2
    "matching",
    "opposite",
    "NS_inconsistent",

    # Layer 3
    "transcriptomics_up",
    "transcriptomics_down",
    "transcriptomics_NS",
]

node_index = {name: i for i, name in enumerate(nodes)}

node_colors = [
    "#B2182B",  # proteomics_up
    "#2166AC",  # proteomics_down

    "#42C478",  # matching
    "#762A83",  # opposite
    "#D3D3D3",  # NS_inconsistent

    "#B2182B",  # transcriptomics_up
    "#2166AC",  # transcriptomics_down
    "#D3D3D3",  # transcriptomics_NS
]

LINK_GREY = "rgba(150,150,150,0.35)"

sources = []
targets = []
values = []
link_colors = []

for _, row in flow.iterrows():
    src = node_index[row["proteomics_direction"]]
    mid = node_index[row["agreement"]]
    tgt = node_index[row["transcriptomics_direction"]]
    count = row["count"]

    # proteomics → agreement
    sources.append(src)
    targets.append(mid)
    values.append(count)
    link_colors.append(LINK_GREY)

    # agreement → transcriptomics
    sources.append(mid)
    targets.append(tgt)
    values.append(count)
    link_colors.append(LINK_GREY)

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=25,
        thickness=15,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=node_colors,
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=link_colors,
    ),
)])

fig.update_layout(
    title_text="Directional Consistency Sankey (DEPs only) — Impute Min/Normalize Median",
    font=dict(size=11),
)

fig.show()


# In[172]:


# Proteomics + Transcriptomics directional colors
COLOR_UP   = "#B2182B"   # dark red
COLOR_DOWN = "#2166AC"   # dark blue
COLOR_NS   = "#D3D3D3"   # light grey

# Agreement colors
COLOR_MATCHING = "#42C478"   # green
COLOR_OPPOSITE = "#762A83"   # purple
COLOR_NSI      = "#D3D3D3"   # grey for NS_inconsistent

# Link color (uniform)
LINK_GREY = "rgba(150,150,150,0.35)"

agreement_collapsed_direct = fold_change_df_impute["combined_consistency"].map({
    "consistent_up":   "matching",
    "consistent_down": "matching",
    "opposite":        "opposite",
    "inconsistent":         "NS_inconsistent",
    "transcriptomics_NS":   "NS_inconsistent",
    "proteomics_NS":        "NS_inconsistent",
    "both_NS":              "NS_inconsistent",
})

df_sankey = fold_change_df_impute.copy()
df_sankey["agreement"] = agreement_collapsed_direct

flow = (
    df_sankey.groupby([
        "proteomics_direction",
        "agreement",
        "transcriptomics_direction"
    ]).size().reset_index(name="count")
)

nodes = [
    # Layer 1: Proteomics
    "proteomics_up",
    "proteomics_down",
    "proteomics_NS",

    # Layer 2: Agreement
    "matching",
    "opposite",
    "NS_inconsistent",

    # Layer 3: Transcriptomics
    "transcriptomics_up",
    "transcriptomics_down",
    "transcriptomics_NS",
]

node_index = {name: i for i, name in enumerate(nodes)}


node_colors = [
    COLOR_UP,        # proteomics_up
    COLOR_DOWN,      # proteomics_down
    COLOR_NS,        # proteomics_NS

    COLOR_MATCHING,  # matching
    COLOR_OPPOSITE,  # opposite
    COLOR_NSI,       # NS_inconsistent

    COLOR_UP,        # transcriptomics_up
    COLOR_DOWN,      # transcriptomics_down
    COLOR_NS,        # transcriptomics_NS
]

sources = []
targets = []
values  = []
link_colors = []

for _, row in flow.iterrows():

    src = node_index[row["proteomics_direction"]]
    mid = node_index[row["agreement"]]
    tgt = node_index[row["transcriptomics_direction"]]
    count = row["count"]

    # proteomics → agreement
    sources.append(src)
    targets.append(mid)
    values.append(count)
    link_colors.append(LINK_GREY)

    # agreement → transcriptomics
    sources.append(mid)
    targets.append(tgt)
    values.append(count)
    link_colors.append(LINK_GREY)


fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=25,
        thickness=15,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=node_colors,
        hovertemplate="%{label}<extra></extra>",
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=link_colors,
        hovertemplate="Count: %{value}<extra></extra>",
    ),
)])

fig.update_layout(
    title_text="Directional Consistency Sankey (all proteins) — Impute Min/Normalize Median",
    font=dict(size=11),
)

fig.show()


