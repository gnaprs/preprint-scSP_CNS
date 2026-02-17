#!/usr/bin/env python
# coding: utf-8

# # package imports

# In[1]:


from scpviz import pAnnData as pAnnData
from scpviz import plotting as scplt
from scpviz import utils as scutils
import scanpy as sc

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(context='paper', style='ticks')

import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D


# # Combined data files
# 
# The same report file is used in Fig 2, Fig 4 and Fig 5 - `pdata_all_region_filtered`

# In[3]:


import os

# obs_columns = ['date', 'acquisition', 'sample_id','size','confirmation','thickness','type','organism','region','well_position']
pdata_all_region = pAnnData.import_data(source_type='diann',report_file = 'data/2505_rna_prot_full/report.tsv')

# data clean-up | uses file_annotation.csv to help with annotations

file_annotation = pd.read_csv('abc_analysis/file_annotation.csv')

file_annotation["parsed_filename"] = (
    file_annotation["RAW FILE"]
    .apply(os.path.basename)
    .str.replace(".raw", "", regex=False)
)

summary = pdata_all_region.summary

summary_files = set(summary.index)
annotation_files = set(file_annotation["parsed_filename"])
missing_in_summary = annotation_files - summary_files
missing_in_annotation = summary_files - annotation_files
pdata_all_region_filtered = pdata_all_region.filter_sample(exclude_file_list=list(missing_in_annotation))

fa_sub = file_annotation[["parsed_filename", "File Name", "Grouping", "Sub grouping", "Region", "Batch"]]
fa_sub = fa_sub.set_index("parsed_filename")
fa_aligned = fa_sub.reindex(pdata_all_region_filtered.summary.index)
cols_to_add = ["File Name", "Grouping", "Sub grouping", "Region", "Batch"]
pdata_all_region_filtered.summary[cols_to_add] = fa_aligned[cols_to_add]

pdata_all_region_filtered.update_summary()

pd.DataFrame(pdata_all_region_filtered.summary)


# In[4]:


pd.DataFrame(pdata_all_region.summary)


# # Fig 1 - size

# In[2]:


pdata = pAnnData.import_data(source_type = 'diann', report_file = 'data/2411_size/report.tsv', obs_columns = ['date','gradient','size','batch','replicate'], prot_value = 'PG.Quantity')

# change size = 10k2 to 10k
pdata.summary['size'] = pdata.summary['size'].replace('10k2','10k')
pdata.update_summary()


# In[3]:


pd.DataFrame(pdata.summary)


# In[4]:


pdata_filter = pdata.filter_sample(min_prot=1000)
pdata_filter = pdata_filter.filter_prot_significant()
pdata_filter = pdata_filter.filter_sample(condition = "(size == 'sc') or (protein_quant > 0.35)")
pdata_filter = pdata_filter.filter_sample(condition = "not(size == '2k' and protein_count > 3000)") # 2000 outlier

pdata_filter = pdata_filter.filter_sample(exclude_file_list=['20240425_Aur60minDIA_10k2_T2_01'])

pd.DataFrame(pdata_filter.summary)


# In[5]:


pdata_filter.summary['size'].value_counts()


# In[ ]:


# # import files.txt
# files = pd.read_csv('data/2411_size/files.txt', sep = '\t', header = None)[0].tolist()
# files.remove('20240425_Aur60minDIA_10k2_T2_01')
# files.remove('20241105_Aur60minDIA_2k_TEAB35_H8')

# pdata_filter = pdata.filter_sample(file_list = files, return_copy = True)
# pd.DataFrame(pdata_filter.summary)


# ### size plots

# In[6]:


# PROTEIN COUNT
amnt_order = ['sc', '2k','5k', '10k', '20k']

# create cmap for replicate from 1 to 5, using blues, and invert the order
# cmap = sns.color_palette("Blues", 4)
cmap = ['#F5793A', '#CFCFCF', '#AFAFAF','#7F7F7F', '#4F4F4F']

sns.barplot(x='size', y='protein_count', data=pdata_filter.summary, errorbar='sd', capsize=.2, saturation=1, palette = cmap, order=amnt_order, hue_order=amnt_order)
sns.stripplot(x='size', y='protein_count', data=pdata_filter.summary, color='black', alpha=0.4, order=amnt_order, s=4.5, jitter=0.15)

# change 'sc' in x-axis to '1000 (single-cell)'
plt.xticks([0, 1, 2, 3, 4], ['800\n(sc)', '2000', '5000', '10000', '20000'])

plt.ylabel('Protein Count')
plt.xlabel('Area ($\mu m^2$)')
plt.ylim(0, 4500)
plt.gcf().set_size_inches(2.63,3)

# Create a second x-axis for number of cells
ax2 = plt.twiny()
ax2.set_xticks([0.45, 1.26, 2.13, 2.98, 3.83, 4.25])  # Shift the x-ticks by 0.5
ax2.set_xticklabels(['1', '2.5', '6.25', '12.5', '25', ''])
ax2.set_xlabel('Number of cells')
ax2.set_xlim(plt.xlim())

# save figure
plt.savefig('MANUSCRIPT/MANUSCRIPT_Fig1_protein_count_bar.svg', bbox_inches="tight", pad_inches=0.1)
plt.show()

# PROTEIN COUNT
amnt_order = ['sc', '2k','5k', '10k', '20k']

# create cmap for replicate from 1 to 5, using blues, and invert the order
# cmap = sns.color_palette("Blues", 4)
cmap = ['#F5793A', '#CFCFCF', '#AFAFAF','#7F7F7F', '#4F4F4F']

sns.violinplot(x='size', y='protein_count', data=pdata_filter.summary, saturation=1, palette = cmap, order=amnt_order, hue_order=amnt_order, edgecolor='black', linewidth=0.8,alpha=0.8, inner="points")

# change 'sc' in x-axis to '1000 (single-cell)'
plt.xticks([0, 1, 2, 3, 4], ['800\n(sc)', '2000', '5000', '10000', '20000'])

plt.ylabel('Protein Count')
plt.xlabel('Area ($\mu m^2$)')
plt.ylim(0, 4500)
plt.gcf().set_size_inches(2.63,3)

# Create a second x-axis for number of cells
ax2 = plt.twiny()
ax2.set_xticks([0.45, 1.26, 2.13, 2.98, 3.83, 4.25])  # Shift the x-ticks by 0.5
ax2.set_xticklabels(['1', '2.5', '6.25', '12.5', '25', ''])
ax2.set_xlabel('Number of cells')
ax2.set_xlim(plt.xlim())

# save figure
plt.show()
# plt.savefig('MANUSCRIPT/MANUSCRIPT_Fig1_protein_count_violin.svg', bbox_inches="tight", pad_inches=0.1)


# In[12]:


# PEPTIDE COUNT
amnt_order = ['sc', '2k','5k', '10k', '20k']

# create cmap for replicate from 1 to 5, using blues, and invert the order
# cmap = sns.color_palette("Blues", 4)
cmap = ['#F5793A', '#CFCFCF', '#AFAFAF','#7F7F7F', '#4F4F4F']

sns.barplot(x='size', y='peptide_count', data=pdata_filter.summary, errorbar='sd', capsize=.1, saturation=1, palette = cmap, order=amnt_order, hue_order=amnt_order)
sns.stripplot(x='size', y='peptide_count', data=pdata_filter.summary, color='black', alpha=0.4, order=amnt_order, s=4.5, jitter=0.15)

# change 'sc' in x-axis to '1000 (single-cell)'
plt.xticks([0, 1, 2, 3, 4], ['800\n(sc)', '2000', '5000', '10000', '20000'])

plt.ylabel('Peptide Count')
plt.xlabel('Area ($\mu m^2$)')
# plt.ylim(0, 4500)
plt.gcf().set_size_inches(2.4,3)

# Create a second x-axis for number of cells
ax2 = plt.twiny()
ax2.set_xticks([0.45, 1.26, 2.13, 2.98, 3.83, 4.25])  # Shift the x-ticks by 0.5
ax2.set_xticklabels(['1', '2.5', '6.25', '12.5', '25', ''])
ax2.set_xlabel('Number of cells')
ax2.set_xlim(plt.xlim())

# save figure
plt.savefig('MANUSCRIPT/MANUSCRIPT_Fig1_peptide_count_bar.svg', bbox_inches="tight", pad_inches=0.1)
plt.show()

# PEPTIDE COUNT

sns.violinplot(x='size', y='peptide_count', data=pdata_filter.summary, saturation=1, palette = cmap, order=amnt_order, hue_order=amnt_order, edgecolor='black', inner="points", linewidth=0.8,alpha=0.8, )

# change 'sc' in x-axis to '1000 (single-cell)'
plt.xticks([0, 1, 2, 3, 4], ['800\n(sc)', '2000', '5000', '10000', '20000'])

plt.ylabel('Peptide Count')
plt.xlabel('Area ($\mu m^2$)')
# plt.ylim(0, 4500)
plt.gcf().set_size_inches(2.63,3)

# Create a second x-axis for number of cells
ax2 = plt.twiny()
ax2.set_xticks([0.45, 1.26, 2.13, 2.98, 3.83, 4.25])  # Shift the x-ticks by 0.5
ax2.set_xticklabels(['1', '2.5', '6.25', '12.5', '25', ''])
ax2.set_xlabel('Number of cells')
ax2.set_xlim(plt.xlim())

# save figure
# plt.show()
plt.savefig('MANUSCRIPT/MANUSCRIPT_Fig1_peptide_count_violin.svg', bbox_inches="tight", pad_inches=0.1)


# In[ ]:


# CV
order = ['sc', '2k', '5k','10k', '20k']
colors = ['#F5793A', '#CFCFCF', '#AFAFAF','#7F7F7F', '#4F4F4F']

classes='size'

fig, ax = plt.subplots(figsize=(2.795, 3))
cv_df = scplt.plot_cv(ax, pdata_filter, classes = classes, return_df=True)

cv_df = cv_df.reset_index()
sns.violinplot(data=cv_df, y='Class', x='CV', orient='h', order=order, palette=colors, linewidth=1, inner='quartile', saturation = 1, ax=ax)
# ax.set_title(classes)
# legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
# plt.legend(handles=legend_patches, bbox_to_anchor=(.75, 1), loc=2, borderaxespad=0., frameon=False)

from matplotlib.ticker import PercentFormatter
ax.xaxis.set_major_formatter(PercentFormatter(1))
plt.yticks([0, 1, 2, 3, 4], ['800\n(sc)', '2000', '5000', '10000', '20000'])
ax.set_ylabel('Area ($\mu m^2$)')
ax.set_xlabel('Protein Abundance CV')


# plt.show()
plt.savefig('MANUSCRIPT/MANUSCRIPT_Fig1_cv.svg', bbox_inches="tight", pad_inches=0.1)


# In[10]:


import upsetplot

size_upset = scutils.get_upset_contents(pdata_filter, classes = 'size')
size_upset=size_upset.reorder_levels(["sc", "2k", "5k", "10k", "20k"])
upplot = upsetplot.UpSet(size_upset, subset_size="count", show_counts=True, facecolor = 'black', sort_categories_by='-input')
upplot.style_subsets(present=["sc"], absent=['2k','5k','10k','20k'],edgecolor='black', facecolor="darkorange", linewidth=2, label="sc only")
upplot.style_subsets(absent=["sc"], present=['2k','5k','10k','20k'],edgecolor = 'black', facecolor="#7F7F7F", linewidth=2, label="in all but sc")
uplot = upplot.plot()
uplot["intersections"].set_ylabel("Subset size")
uplot["totals"].set_xlabel("Protein count")

legend = uplot["intersections"].get_legend()
legend.get_frame().set_facecolor('white') 

# set size
plt.gcf().set_size_inches(10.3,3)
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-UpsetPlot.svg", bbox_inches="tight")
# plt.show()


# In[13]:


prot_sc_df = scutils.get_upset_query(size_upset, present=["sc"], absent=['2k','5k','10k','20k'])
prot_not_sc_df = scutils.get_upset_query(size_upset, absent=["sc"], present=['2k','5k','10k','20k'])
prot_all_df = scutils.get_upset_query(size_upset, present=['sc','2k','5k','10k','20k'], absent=[])


# In[ ]:


# prep data for area

df = pdata_filter.get_abundance(namelist=prot_all_df.accession.to_list())

# Map size → numeric area
area_map = {
    "sc": 800,
    "2k": 2000,
    "5k": 5000,
    "10k": 10000,
    "20k": 20000,
}
df["area"] = df["size"].map(area_map)

# Aggregate + log transforms
agg = (
    df.groupby(["accession", "area", "gene"])["abundance"]
      .mean()
      .reset_index()
)

agg["log_ab"] = np.log10(agg["abundance"])
agg["log_area"] = np.log10(agg["area"])

# Compute slopes (\beta_1) 
slopes = agg.groupby("accession").apply(
    lambda x: np.polyfit(x["log_area"], x["log_ab"], 1)[0]
)

# Mean abundance (per protein)
mean_ab = (
    df.groupby(["accession", "gene"])["abundance"]
      .mean()
      .reset_index()
)

# Remove keratins, keep to proteins of reasonable slopes
keratins = mean_ab[mean_ab["gene"].str.startswith("Krt", na=False)]["accession"].unique()
mean_ab_nonker = mean_ab[~mean_ab["accession"].isin(keratins)].copy()
mean_ab_nonker["slope"] = mean_ab_nonker["accession"].map(slopes)
reasonable = mean_ab_nonker.query("0.8 < slope < 1.1").copy()

sorted_nonker = reasonable.sort_values("abundance")

# Low abundance example (Plxnb1 since it's more biologically relevant here)
low_prot_acc  = mean_ab[mean_ab["gene"] == "Plxnb1"]["accession"].iloc[0]
low_prot_gene = "Plxnb1"

# Medium abundance
mid = sorted_nonker.iloc[len(sorted_nonker)//2 + 200]
mid_prot_acc, mid_prot_gene = mid["accession"], mid["gene"]

# High abundance
high = sorted_nonker.iloc[-1]
high_prot_acc, high_prot_gene = high["accession"], high["gene"]

print("Low:   ", low_prot_acc, low_prot_gene, slopes[low_prot_acc])
print("Medium:", mid_prot_acc, mid_prot_gene, slopes[mid_prot_acc])
print("High:  ", high_prot_acc, high_prot_gene, slopes[high_prot_acc])

keratin_slopes = slopes[slopes.index.isin(keratins)]
keratin_mean = keratin_slopes.mean()

nonkeratin_slopes = slopes[~slopes.index.isin(keratins)]
nonkeratin_mean = nonkeratin_slopes.mean()
textstr = f"Mean slope (non-keratin): {nonkeratin_mean:.2f}"

highlight = {
    "Low-Abundance":    low_prot_acc,
    "Medium-Abundance": mid_prot_acc,
    "High-Abundance":   high_prot_acc,
}


# In[ ]:


from matplotlib.transforms import offset_copy
import matplotlib.lines as mlines

palette = {
    "Low-Abundance":    "#6baed6",
    "Medium-Abundance": "#3182bd",
    "High-Abundance":   "#08519c",
}

fig, ax = plt.subplots(figsize=(3, 3))

# Background grey (non-keratin)
for acc, sub in agg.groupby("accession"):
    if acc in keratins:
        continue
    ax.plot(sub["area"], sub["log_ab"],
            color="lightgrey", linewidth=0.5, alpha=0.07)

# Highlight: low, medium, high abundance proteins
highlight_lines = []  # for inline labels

for key, acc in highlight.items():
    gene = mean_ab.loc[mean_ab["accession"] == acc, "gene"].values[0]
    beta = slopes[acc]

    sub = agg[agg.accession == acc]
    (line,) = ax.plot(
        sub["area"], sub["log_ab"],
        linewidth=1.5, alpha=0.9,
        marker="o", markersize=3,
        color=palette[key],
        label=key,   # legend entry: Low/Medium/High-Abundance
    )
    highlight_lines.append((line, gene, beta))

# Inline labels above each line with gene name + β
for line, gene, beta in highlight_lines:
    x = line.get_xdata()
    y = line.get_ydata()
    x_end, y_end = x[-1], y[-1]

    text_transform = offset_copy(ax.transData, fig=fig, x=3, y=3, units="points")

    ax.text(
        x_end, y_end,
        f"{gene} (β={beta:.2f})",
        transform=text_transform,
        color="black",
        fontsize=7,
        ha="right", va="bottom",
    )

# Rounded Mean β box (your style preserved)
props = dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="grey")
ax.text(
    0.05, 0.95,
    f"Mean β: {nonkeratin_mean:.2f}",
    transform=ax.transAxes,
    fontsize=8,
    va="top",
    bbox=props,
)

# Compact legend inside plot (only the 3 abundance classes)
legend = ax.legend(
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(0.03, 0.9),
    fontsize=7,
)

ax.set_xscale("log")
ax.set_xlabel("Sampled Area (µm²)")
ax.set_ylabel(r"$\log_{10}$ Abundance")
ax.set_ylim(top=10)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig1-Abundance_Area.svg", bbox_inches="tight")


# In[ ]:


# KERATIN FIGURE: potentially supplementary?

palette = {
    "Low-Abundance":    "#6baed6",  # light blue
    "Medium-Abundance": "#3182bd",  # medium blue
    "High-Abundance":   "#08519c",  # dark blue
}

# PLOT

fig, ax = plt.subplots(figsize=(2.64,3))

# Background grey (non-keratin)
for acc, sub in agg.groupby("accession"):
    if acc in keratins:
        continue
    ax.plot(sub["area"], sub["log_ab"],
            color="lightgrey", linewidth=0.5, alpha=0.07)

keratin_color = "#4a2c2a"   # muted red
keratin_alpha = 0.25
keratin_lw = 0.9

# Keratins (orange)
for acc in keratins:
    sub = agg[agg["accession"] == acc]
    ax.plot(sub["area"], sub["log_ab"],
            color=keratin_color, linewidth=keratin_lw, alpha=keratin_alpha)

# Highlight: low, medium, high abundance proteins
# for key, acc in highlight.items():
#     gene = mean_ab.loc[mean_ab["accession"] == acc, "gene"].values[0]
#     beta = slopes[acc]
#     label = f"{key}: {acc}/{gene} (β={beta:.2f})"

#     sub = agg[agg.accession == acc]
#     ax.plot(
#         sub["area"], sub["log_ab"],
#         linewidth=1.5, alpha=0.9,
#         marker="o", markersize=3,
#         color=palette[key],
#         label=label
#     )

# Add slope text box
props = dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="grey")
ax.text(
    0.05, 0.92, 
    f"Mean β (Keratins): {keratin_mean:.2f}",
    transform=ax.transAxes,
    fontsize=8,
    verticalalignment="bottom",
    bbox=props
)

ax.set_xscale("log")
ax.set_xlabel("Sampled Area (µm²)")
ax.set_ylabel("$\log_{10}$ Abundance")
ax.set_ylim(top=10)
# ax.set_title("Scaling Behaviour of Protein Abundance")
ax.legend(frameon=False)

scplt.shift_legend(ax)
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-Abundance_Area_Keratin.svg", bbox_inches="tight")


# In[14]:


order = ['sc','2k', '5k','10k', '20k']

# colors = sns.color_palette("Blues", 4)
cmap = ['#F5793A', '#CFCFCF', '#AFAFAF','#7F7F7F', '#4F4F4F']
colors = ['#F5793A', '#CFCFCF', '#AFAFAF','#7F7F7F', '#4F4F4F']
cmaps = ['Oranges','Greys','Greys','Greys','Greys']
# cmaps = []
# for i, color in enumerate(colors):
#     cmap = LinearSegmentedColormap.from_list(f'custom_blues_{i}', ['white', color])
#     cmaps.append(cmap)
ylims = (2,10.5)


fig, ax = plt.subplots(figsize=(6.88,3), nrows = 1, ncols = 2, sharey='col', sharex='col', gridspec_kw={'width_ratios': [14.61, 9.03], 'hspace': 0.08, 'wspace': 0})
labels = order
line_width = 1

ax = plt.subplot(1, 2, 1)
scplt.plot_rankquant(ax, pdata_filter, classes = 'size', order = order, cmap = cmaps, color=colors, s=10, calpha=1, alpha=0.005)
plt.ylabel('MS Intensity')
ax.set_ylim(10**ylims[0],10**ylims[1])
legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
plt.legend(handles=legend_patches, bbox_to_anchor=(.75, 1), loc=2, borderaxespad=0., frameon=False)
scplt.mark_rankquant(ax, pdata_filter, mark_df=prot_sc_df, class_values=['sc'], show_label=False, color = 'black', label_type='gene')
scplt.mark_rankquant(ax, pdata_filter, mark_df=prot_not_sc_df, class_values=['20k'], s = 5,show_label=False, color = 'white', label_type='gene', alpha=0.5)

ax = plt.subplot(1, 2, 2)
scplt.plot_raincloud(ax, pdata_filter, classes='size', order=order, color=colors, linewidth=line_width, debug = False)
scplt.mark_raincloud(ax, pdata_filter, mark_df=prot_sc_df, class_values=['sc'], color = 'black')
scplt.mark_raincloud(ax, pdata_filter, mark_df=prot_not_sc_df, class_values=['20k'], color = 'white', alpha=0.3, lowest_index = 4)
ax.set_yticks([],[])
# if multiple labels
# ax.set_xticks(np.arange(len(order))+1, [case[1]+'\n'+case[2] for case in cases])

order_label = ['sc','2000', '5000','10000', '20000']
ax.set_xticks(np.arange(len(order))+1, order_label)
ax.set_ylim(ylims[0],ylims[1])
ax.set_xlabel("Area (µm²)")

plt.subplots_adjust(wspace=0, hspace=0)
sns.despine()

# save as svg
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig1-RankQuant.png", dpi=450, bbox_inches="tight")
# plt.show()


# #### other exploration

# In[ ]:


fig, ax = plt.subplots(figsize=(4, 3))

beta_df = mean_ab.copy()
beta_df["beta"] = beta_df["accession"].map(slopes)

is_keratin = beta_df["accession"].isin(keratins)
nk = beta_df[~is_keratin]      # non-keratin
k  = beta_df[ is_keratin]      # keratin

ax.scatter(
    np.log10(nk["abundance"]),
    nk["beta"],
    s=6,
    alpha=0.5,
    color="steelblue",
    label="Non-keratin",
)

ax.scatter(
    np.log10(k["abundance"]),
    k["beta"],
    s=10,
    alpha=0.8,
    facecolor="orange",
    edgecolor="black",
    label="Keratin",
)

ax.set_xlabel("log10 Mean Abundance")
ax.set_ylabel("β (scaling exponent)")
ax.set_title("Scaling Exponent vs. Protein Abundance")
ax.grid(alpha=0.2)
ax.legend(frameon=False, fontsize=8, loc="lower right")

plt.show()


# In[ ]:


all_slopes = np.concatenate([nonkeratin_slopes, keratin_slopes])
bins = np.linspace(all_slopes.min(), all_slopes.max(), 30)

fig, ax = plt.subplots(figsize=(3.36, 3))

# Non-keratin: filled
ax.hist(
    nonkeratin_slopes,
    bins=bins,
    alpha=0.6,
    color="steelblue",
    label="Non-keratin",
)

# Keratin: outline-only, same bins
ax.hist(
    keratin_slopes,
    bins=bins,
    histtype="step",
    linewidth=1.8,
    color="orange",
    label="Keratin",
)

# # Optional: rug plot for keratins at the bottom
# for x in keratin_slopes:
#     ax.axvline(x, ymin=0, ymax=0.03, color="orange", linewidth=0.8, alpha=0.9)

# Medians
nk_med = np.median(nonkeratin_slopes)
k_med  = np.median(keratin_slopes)

ax.axvline(
    nk_med,
    color="steelblue",
    linestyle="--",
    linewidth=1.5,
    label=f"Non-keratin median = {nk_med:.2f}",
)
ax.axvline(
    k_med,
    color="orange",
    linestyle="--",
    linewidth=1.5,
    label=f"Keratin median = {k_med:.2f}",
)

ax.set_xlabel("β (scaling exponent)")
ax.set_ylabel("Count")
ax.set_title("Scaling Exponent Distributions")
ax.grid(alpha=0.2)

# Remove duplicate legend entries
handles, labels = ax.get_legend_handles_labels()
unique = dict(zip(labels, handles))
ax.legend(unique.values(), unique.keys(), frameon=False, fontsize=8, loc="upper left")

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-Histogram_Keratin.svg", bbox_inches="tight")


# ### marker - bulk

# In[ ]:


import re

# (plasma membrane separate; organelle membranes go with organelle)
selected_components = [
    "cytoplasm",
    "plasma membrane",
    "extracellular region",
    "mitochondrion",
    "nucleus",
    "endoplasmic reticulum",
    "golgi apparatus",
    "cytoskeleton",
]

label_map = {
    "cytoplasm": "Cytoplasm",
    "plasma membrane": "Plasma membrane",
    "extracellular region": "Extracellular region",
    "mitochondrion": "Mitochondria",
    "nucleus": "Nucleus",
    "endoplasmic reticulum": "ER",
    "golgi apparatus": "Golgi apparatus",
    "cytoskeleton": "Cytoskeleton",
    "others": "Others",
}

plot_order = selected_components + ["others"]

def _clean_cc_term(s: str) -> str:
    if pd.isna(s):
        return ""
    s = str(s).strip()
    s = re.sub(r"\[.*?\]", "", s).strip()
    return s

def map_cc_term(term: str) -> str:
    t = term.lower()

    # cytoskeleton (keep separate)
    if "cytoskeleton" in t:
        return "cytoskeleton"

    # nucleus
    if re.search(r"\bnucleus\b|nucleoplasm|nuclear", t):
        return "nucleus"

    # cytoplasm (cytosol/cytoplasmic included)
    if re.search(r"cytosol|cytoplasm|cytoplasmic", t):
        return "cytoplasm"

    # extracellular
    if "extracellular" in t:
        return "extracellular region"

    # ER (includes ER membrane etc.)
    if "endoplasmic reticulum" in t:
        return "endoplasmic reticulum"

    # Golgi
    if "golgi" in t:
        return "golgi apparatus"

    # mitochondrion (includes mitochondrial membrane etc.)
    if "mitochondri" in t:
        return "mitochondrion"

    # plasma membrane (keep separate from generic "membrane")
    if "plasma membrane" in t:
        return "plasma membrane"

    # everything else
    return "others"

def build_cc_df(df, cc_col="gene_ontology_cellular_component", protein_col=None, source_name=None):
    prot = df.index.to_series() if protein_col is None else df[protein_col]
    out = df[[cc_col]].copy()
    out["protein"] = prot.values
    out = out.rename(columns={cc_col: "cc_raw"})
    if source_name is not None:
        out["source"] = source_name

    out["cc_raw"] = out["cc_raw"].astype(str).str.split(";")
    out = out.explode("cc_raw", ignore_index=True)

    out["cc_clean"] = out["cc_raw"].map(_clean_cc_term)
    out = out[out["cc_clean"].ne("") & out["cc_clean"].ne("nan")].copy()

    out["cc_mapped"] = out["cc_clean"].map(map_cc_term)
    out["in_selected"] = out["cc_mapped"].isin(selected_components)
    return out

def counts_from_cc_df(cc_df):
    return cc_df["cc_mapped"].value_counts()


# #### IF (ab vs pbs)

# In[66]:


# Sample, n/a, 100ng, n/a, n/a, n/a, n/a, PBS3, n/a
obs_columns = ['sample', 'load', 'pair_sub']

pdata_ab = pAnnData.import_proteomeDiscoverer(prot_file = 'data/MouseFF_100ng_AbPBS_43minDDA_20230616_Proteins.txt', pep_file='data/MouseFF_100ng_AbPBS_43minDDA_20230616_PeptideGroups.txt', obs_columns=obs_columns)


# In[67]:


pdata_ab.summary['pair'] = pdata_ab.summary['pair_sub'].str[-1]
pdata_ab.summary['type'] = pdata_ab.summary['pair_sub'].str[:-1]

# for column "type", change Ab to IF and PBS to control
pdata_ab.summary['type'] = pdata_ab.summary['type'].replace({'Ab': 'IF', 'PBS': 'Control'})

# move pair and type to the front
cols = pdata_ab.summary.columns.tolist()
cols = cols[-2:] + cols[:-2]
pdata_ab.summary = pdata_ab.summary[cols]

pdata_ab.update_summary()

pd.DataFrame(pdata_ab.summary)


# In[68]:


averaged_data = pdata_ab.summary.groupby(['pair', 'type'])['protein_count'].mean().reset_index()
pivoted_data = averaged_data.pivot(index='pair', columns='type', values='protein_count')

print(pivoted_data)

set1 = pivoted_data['IF'].values
set2 = pivoted_data['Control'].values

data = pd.melt(pivoted_data)


# In[69]:


ab_summary = pdata_ab.summary

# make color palette for Ab and PBS
colors = sns.color_palette("Paired")[2:4]
treatment_order = ['Control', 'IF']
color_dict = dict(zip(treatment_order, colors))


# Plot
fig, ax = plt.subplots(figsize=(1.74, 2.13))
sns.swarmplot(data=data, x="type", y="value", ax=ax, color='k', order = treatment_order)
sns.barplot(data=data, x="type", y="value", ax=ax, errorbar = 'ci', alpha=1, palette=color_dict, order = treatment_order)

# Now connect the dots
# Find idx0 and idx1 by inspecting the elements return from ax.get_children()
# ... or find a way to automate it
idx0 = 0
idx1 = 1
locs1 = ax.get_children()[idx0].get_offsets()
locs2 = ax.get_children()[idx1].get_offsets()

# before plotting, we need to sort so that the data points
# correspond to each other as they did in "set1" and "set2"
sort_idxs1 = np.argsort(set1)
sort_idxs2 = np.argsort(set2)

# revert "ascending sort" through sort_idxs2.argsort(),
# and then sort into order corresponding with set1
locs2_sorted = locs2[sort_idxs2.argsort()][sort_idxs1]

for i in range(locs1.shape[0]):
    x = [locs1[i, 0], locs2_sorted[i, 0]]
    y = [locs1[i, 1], locs2_sorted[i, 1]]
    ax.plot(x, y, color="black", alpha=0.3)

# calculate p-value for significance of "control" vs "stained" protein_count
from scipy.stats import ttest_rel
control = pivoted_data['Control']
stained = pivoted_data['IF']

# scplt.plot_significance(ax, 2450, 20, pval=ttest_rel(control, stained).pvalue, fontsize=8)

plt.ylabel('Protein count')
plt.xlabel('Paired Sample')

plt.ylim(2000, 2550)
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-IF_Bulk_ProtCount.svg", bbox_inches="tight")


# In[70]:


fig, ax = plt.subplots(figsize=(3, 3))

scplt.plot_venn(ax, pdata_ab, classes = 'type', set_colors = sns.color_palette("Paired")[2:4], label_order = ['Control', 'IF'])
plt.rcParams.update({'font.size': 12})
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-IF_Bulk_Venn.svg", bbox_inches="tight")

plt.rcParams.update({"font.size": 14})

fig, ax = plt.subplots(figsize=(3, 3))
ax, contents = scplt.plot_venn(
    ax, pdata_ab,
    classes="type",
    set_colors=sns.color_palette("Paired")[2:4],
    label_order=["Control", "IF"],
    weighted=True,
    return_contents=True,
    fixed_subset_sizes=(1, 1, 3),
)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Control", "IF"}]:
    x, y = t.get_position()
    if t.get_text() == "Control":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")

plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-IF_Bulk_Venn_Weighted.svg", bbox_inches="tight")


# In[72]:


ab_upset = scutils.get_upset_contents(pdata_ab, classes = ['type'])

pdata_ab_only = pdata_ab.filter_sample(condition='type == "IF"', return_copy=True)
pdata_pbs_only = pdata_ab.filter_sample(condition='type == "Control"', return_copy=True)

ab_only = scutils.get_upset_query(ab_upset, present='IF', absent='Control')
pbs_only = scutils.get_upset_query(ab_upset, present='Control', absent='IF')

# searching 2000+ proteins on Uniprot, this will take awhile...
ab_all = scutils.get_upset_query(ab_upset, present='IF', absent=[])
pbs_all = scutils.get_upset_query(ab_upset, present='Control', absent=[])

# export ab_only and pbs_only to csv
# ab_only.to_csv('ab_only.csv')
# pbs_only.to_csv('pbs_only.csv')


# In[73]:


from matplotlib.colors import LinearSegmentedColormap

light_green_cmap = mcolors.LinearSegmentedColormap.from_list("light_green_cmap", ["#E8F5E9", "#388E3C"])
dark_green_cmap = mcolors.LinearSegmentedColormap.from_list("dark_green_cmap", ["#1B5E20", "#004D00"])

colors = sns.color_palette("Paired")[2:4]
order = ['Control', 'IF']
cmaps = [light_green_cmap,dark_green_cmap]
ylims = (3.5,10.5)

fig, ax = plt.subplots(figsize=(5.62,2.13), nrows = 1, ncols = 2, sharey='col', sharex='col', gridspec_kw={'width_ratios': [11.34, 6.85], 'hspace': 0.08, 'wspace': 0})
labels = order
line_width = 1

ax = plt.subplot(1, 2, 1)
scplt.plot_rankquant(ax, pdata_ab, classes = ['type'], order = order, cmap = cmaps, color=colors, s=10, calpha=1, alpha=0.005)
plt.ylabel('MS Intensity')
ax.set_ylim(10**ylims[0],10**ylims[1])
legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
plt.legend(handles=legend_patches, bbox_to_anchor=(.6, 1), loc=2, borderaxespad=0., frameon=False)
scplt.mark_rankquant(ax, pdata_ab, mark_df=pbs_only, class_values=['Control'], show_label=False, color = 'black', alpha=0.3, label_type='gene')
scplt.mark_rankquant(ax, pdata_ab, mark_df=ab_only, class_values=['IF'], show_label=False, color = 'white', alpha=0.3, label_type='gene')

ax = plt.subplot(1, 2, 2)
scplt.plot_raincloud(ax, pdata_ab, classes=['type'], order=order, color=colors, linewidth=line_width)
scplt.mark_raincloud(ax, pdata_ab, mark_df=pbs_only, class_values=['Control'], color = 'black', alpha=0.3)
scplt.mark_raincloud(ax, pdata_ab, mark_df=ab_only, class_values=['IF'], color = 'white', alpha=0.3, lowest_index = 1)
ax.set_yticks([],[])
ax.set_xticks(np.arange(len(order))+1, [case for case in order])
ax.set_ylim(ylims[0],ylims[1])

plt.subplots_adjust(wspace=0, hspace=0)
sns.despine()

# save as svg
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-IF_Bulk_RankQuant.png", dpi=450, bbox_inches="tight")
# plt.show()


# In[74]:


# check df
ab_cc_df  = build_cc_df(ab_all,  protein_col='id', source_name="ab")
pbs_cc_df = build_cc_df(pbs_all, protein_col='id', source_name="pbs")

# counts 
ab_counts  = counts_from_cc_df(ab_cc_df).reindex(plot_order, fill_value=0)
pbs_counts = counts_from_cc_df(pbs_cc_df).reindex(plot_order, fill_value=0)

# convert to %
ab_pct  = ab_counts  / ab_counts.sum()  * 100
pbs_pct = pbs_counts / pbs_counts.sum() * 100

# order top→bottom like your plot
labels = plot_order[::-1]
display_labels = [label_map[l] for l in labels]
ab_pct  = ab_pct.reindex(labels)
pbs_pct = pbs_pct.reindex(labels)

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

# bar chart
fig, ax = plt.subplots(figsize=(2.7, 2.7))

y = np.arange(len(labels))
h = 0.35

colors = sns.color_palette("Paired")[2:4]

ax.barh(y - h/2, ab_pct.values, height=h, label="IF",
        color=colors[1], edgecolor="black", linewidth=0.6)
ax.barh(y + h/2, pbs_pct.values,  height=h, label="Control",
        color=colors[0], edgecolor="black", linewidth=0.6)

ax.set_yticks(y)
ax.set_yticklabels(display_labels)
ax.set_xlabel("% proteins")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", width=1.0, length=8)

ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-IF_CC.svg", bbox_inches="tight")


# In[79]:


ab_cc_df


# #### HCR (hcr vs pbs)

# In[2]:


# Sample, n/a, 100ng, n/a, n/a, n/a, n/a, PBS3, n/a
obs_columns = ['sample', 'treatment', 'type']

pdata_hcr = pAnnData.import_proteomeDiscoverer(prot_file = 'data/Marion_20241021_OTE_Aur60min_FFmouse_HCR_Proteins.txt', pep_file='data/Marion_20241021_OTE_Aur60min_FFmouse_HCR_PeptideGroups.txt', obs_columns=obs_columns)


# In[3]:


# rename "stained" to "HCR" for the treatment column
pdata_hcr.summary['treatment'] = pdata_hcr.summary['treatment'].replace('stained', 'HCR')
pdata_hcr.summary['treatment'] = pdata_hcr.summary['treatment'].replace('control', 'Control')
pdata_hcr.update_summary()


# In[5]:


hcr_summary = pdata_hcr.summary

treatment_order = ['Control', 'HCR']

# make color palette for Ab and PBS
colors = sns.color_palette("Paired")[4:6]
color_dict = dict(zip(treatment_order, colors))

# Plot
fig, ax = plt.subplots(figsize=(1.74, 2.13))
sns.swarmplot(data=hcr_summary, x="treatment", y="protein_count", ax=ax, color='k')
sns.barplot(data=hcr_summary, x="treatment", y="protein_count", ax=ax, errorbar = 'ci', alpha=1, palette=color_dict)

# calculate p-value for significance of "control" vs "stained" protein_count
from scipy.stats import ttest_ind
control = hcr_summary[hcr_summary['treatment'] == 'Control']['protein_count']
stained = hcr_summary[hcr_summary['treatment'] == 'HCR']['protein_count']

# scplt.plot_significance(ax, 2630, 30, pval=ttest_ind(control, stained).pvalue, fontsize=8)


# Add significance
plt.ylabel('Protein count')
plt.xlabel('')

plt.ylim(2000, 2800)
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-HCR_Bulk_ProtCount.svg", bbox_inches="tight")


# In[24]:


fig, ax = plt.subplots(figsize=(3, 3))

scplt.plot_venn(ax, pdata_hcr, classes = 'treatment', set_colors = sns.color_palette("Paired")[4:6])
plt.rcParams.update({'font.size': 12})
# plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-HCR_Bulk_Venn.svg", bbox_inches="tight")
plt.show()


# In[6]:


plt.rcParams.update({"font.size": 14})

fig, ax = plt.subplots(figsize=(3, 3))
ax, contents = scplt.plot_venn(
    ax, pdata_hcr,
    classes="treatment",
    set_colors=sns.color_palette("Paired")[4:6],
    weighted=True,
    fixed_subset_sizes=(1, 1, 3),
    return_contents=True,
)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Control", "IF"}]:
    x, y = t.get_position()
    if t.get_text() == "Control":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")

plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-HCR_Bulk_Venn_Weighted.svg", bbox_inches="tight")
# plt.show()


# In[7]:


hcr_upset = scutils.get_upset_contents(pdata_hcr, classes = ['treatment'])

pdata_hcr_only = pdata_hcr.filter_sample(condition='treatment == "HCR"', return_copy=True)
pdata_pbs_only = pdata_hcr.filter_sample(condition='treatment == "Control"', return_copy=True)

hcr_only = scutils.get_upset_query(hcr_upset, present='HCR', absent='Control')
pbs_only = scutils.get_upset_query(hcr_upset, present='Control', absent='HCR')

# searching 2000+ proteins on Uniprot, this will take awhile...
hcr_all = scutils.get_upset_query(hcr_upset, present='HCR', absent=None)
pbs_all = scutils.get_upset_query(hcr_upset, present='Control', absent=None)


# In[8]:


from matplotlib.colors import LinearSegmentedColormap

light_red_cmap = mcolors.LinearSegmentedColormap.from_list("light_red_cmap", ["#FFEBEE", "#D32F2F"])
dark_red_cmap = mcolors.LinearSegmentedColormap.from_list("dark_red_cmap", ["#B71C1C", "#4A0000"])

colors = sns.color_palette("Paired")[4:6]
order = ['Control', 'HCR']
cmaps = [light_red_cmap,dark_red_cmap]
ylims = (3,10.5)

fig, ax = plt.subplots(figsize=(5.62,2.13), nrows = 1, ncols = 2, sharey='col', sharex='col', gridspec_kw={'width_ratios': [11.34, 6.85], 'hspace': 0.08, 'wspace': 0})
labels = order
line_width = 1

ax = plt.subplot(1, 2, 1)
scplt.plot_rankquant(ax, pdata_hcr, classes = ['treatment'], order = order, cmap = cmaps, color=colors, s=10, calpha=1, alpha=0.005)
plt.ylabel('MS Intensity')
ax.set_ylim(10**ylims[0],10**ylims[1])
legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
plt.legend(handles=legend_patches, bbox_to_anchor=(.6, 1), loc=2, borderaxespad=0., frameon=False)
scplt.mark_rankquant(ax, pdata_hcr, mark_df=pbs_only, class_values=['Control'], show_label=False, color = 'black', alpha=0.3, label_type='gene')
scplt.mark_rankquant(ax, pdata_hcr, mark_df=hcr_only, class_values=['HCR'], show_label=False, color = 'white', alpha=0.3, label_type='gene')

ax = plt.subplot(1, 2, 2)
scplt.plot_raincloud(ax, pdata_hcr, classes=['treatment'], order=order, color=colors, linewidth=line_width)
scplt.mark_raincloud(ax, pdata_hcr, mark_df=pbs_only, class_values=['Control'], color = 'black', alpha=0.3)
scplt.mark_raincloud(ax, pdata_hcr, mark_df=hcr_only, class_values=['HCR'], color = 'white', alpha=0.3, lowest_index = 1)
ax.set_yticks([],[])
ax.set_xticks(np.arange(len(order))+1, [case for case in order])
ax.set_ylim(ylims[0],ylims[1])

plt.subplots_adjust(wspace=0, hspace=0)
sns.despine()

# save as svg
# plt.savefig('ab_rankquant.svg', dpi=300, bbox_inches='tight')
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-HCR_Bulk_RankQuant.png", dpi=450, bbox_inches="tight")
# plt.show()


# In[64]:


# check df
hcr_cc_df  = build_cc_df(hcr_all,  protein_col='id', source_name="hcr")
pbs_cc_df = build_cc_df(pbs_all, protein_col='id', source_name="pbs")

# counts 
hcr_counts  = counts_from_cc_df(hcr_cc_df).reindex(plot_order, fill_value=0)
pbs_counts = counts_from_cc_df(pbs_cc_df).reindex(plot_order, fill_value=0)

# convert to %
hcr_pct  = hcr_counts  / hcr_counts.sum()  * 100
pbs_pct = pbs_counts / pbs_counts.sum() * 100

# order top→bottom like your plot
labels = plot_order[::-1]
display_labels = [label_map[l] for l in labels]
hcr_pct  = hcr_pct.reindex(labels)
pbs_pct = pbs_pct.reindex(labels)

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

# bar chart
fig, ax = plt.subplots(figsize=(2.7, 2.7))

y = np.arange(len(labels))
h = 0.35

colors = sns.color_palette("Paired")[4:6]

ax.barh(y - h/2, hcr_pct.values, height=h, label="HCR",
        color=colors[1], edgecolor="black", linewidth=0.6)
ax.barh(y + h/2, pbs_pct.values,  height=h, label="Control",
        color=colors[0], edgecolor="black", linewidth=0.6)

ax.set_yticks(y)
ax.set_yticklabels(display_labels)
ax.set_xlabel("% proteins")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", width=1.0, length=8)

ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-HCR_CC.svg", bbox_inches="tight")


# ### marker - sc

# In[2]:


import re

# (plasma membrane separate; organelle membranes go with organelle)
selected_components = [
    "cytoplasm",
    "plasma membrane",
    "extracellular region",
    "mitochondrion",
    "nucleus",
    "endoplasmic reticulum",
    "golgi apparatus",
    "cytoskeleton",
]

label_map = {
    "cytoplasm": "Cytoplasm",
    "plasma membrane": "Plasma membrane",
    "extracellular region": "Extracellular region",
    "mitochondrion": "Mitochondria",
    "nucleus": "Nucleus",
    "endoplasmic reticulum": "ER",
    "golgi apparatus": "Golgi apparatus",
    "cytoskeleton": "Cytoskeleton",
    "others": "Others",
}

plot_order = selected_components + ["others"]

def _clean_cc_term(s: str) -> str:
    if pd.isna(s):
        return ""
    s = str(s).strip()
    s = re.sub(r"\[.*?\]", "", s).strip()
    return s

def map_cc_term(term: str) -> str:
    t = term.lower()

    # cytoskeleton (keep separate)
    if "cytoskeleton" in t:
        return "cytoskeleton"

    # nucleus
    if re.search(r"\bnucleus\b|nucleoplasm|nuclear", t):
        return "nucleus"

    # cytoplasm (cytosol/cytoplasmic included)
    if re.search(r"cytosol|cytoplasm|cytoplasmic", t):
        return "cytoplasm"

    # extracellular
    if "extracellular" in t:
        return "extracellular region"

    # ER (includes ER membrane etc.)
    if "endoplasmic reticulum" in t:
        return "endoplasmic reticulum"

    # Golgi
    if "golgi" in t:
        return "golgi apparatus"

    # mitochondrion (includes mitochondrial membrane etc.)
    if "mitochondri" in t:
        return "mitochondrion"

    # plasma membrane (keep separate from generic "membrane")
    if "plasma membrane" in t:
        return "plasma membrane"

    # everything else
    return "others"

def build_cc_df(df, cc_col="gene_ontology_cellular_component", protein_col=None, source_name=None):
    prot = df.index.to_series() if protein_col is None else df[protein_col]
    out = df[[cc_col]].copy()
    out["protein"] = prot.values
    out = out.rename(columns={cc_col: "cc_raw"})
    if source_name is not None:
        out["source"] = source_name

    out["cc_raw"] = out["cc_raw"].astype(str).str.split(";")
    out = out.explode("cc_raw", ignore_index=True)

    out["cc_clean"] = out["cc_raw"].map(_clean_cc_term)
    out = out[out["cc_clean"].ne("") & out["cc_clean"].ne("nan")].copy()

    out["cc_mapped"] = out["cc_clean"].map(map_cc_term)
    out["in_selected"] = out["cc_mapped"].isin(selected_components)
    return out

def counts_from_cc_df(cc_df):
    return cc_df["cc_mapped"].value_counts()


# In[3]:


pdata_fig1_supplementary = pAnnData.import_data(source_type='diann', report_file = 'data/2512_Ab-HCR/report.parquet')


# In[4]:


pdata_fig1_supplementary.summary = scutils.parse_filename_index(pdata_fig1_supplementary.summary, obs_columns=['date','acquisition','size','buffer','well_position'], condition='parsingType == "5-tokens"')
pdata_fig1_supplementary.summary = scutils.parse_filename_index(pdata_fig1_supplementary.summary, obs_columns=['date','acquisition','sample_id','size','confirmation','thickness','type','organism','region','well_position'], condition='parsingType == "10-tokens"')

mapping = {"800": "sc",}
pdata_fig1_supplementary.summary["size"] = pdata_fig1_supplementary.summary["size"].replace(mapping)
pdata_fig1_supplementary.update_summary()

pdata_fig1_supplementary_sc = pdata_fig1_supplementary.filter_sample(condition="size=='sc'")

# for dataset with single cell stained vs non-stained
# ab stained is april, pbs stained is november
summary = pdata_fig1_supplementary_sc.summary
summary.loc[summary["date"].astype(str).str.startswith("202404"), "type"] = "IF"
summary.loc[summary["date"].astype(str).str.startswith("202411"), "type"] = "Control"

pdata_fig1_supplementary_sc.summary = summary
pdata_fig1_supplementary_sc.update_summary()

pd.DataFrame(pdata_fig1_supplementary_sc.summary)


# In[5]:


pdata_fig1_supplementary_sc.summary.type.value_counts()


# #### IF (ab vs pbs)

# In[6]:


pdata_ab_sc = pdata_fig1_supplementary_sc.filter_sample(condition = 'type in ("IF", "Control")')
pdata_ab_sc = pdata_ab_sc.filter_sample(min_prot = 1000)
pdata_ab_sc = pdata_ab_sc.filter_prot_significant()

pdata_ab_sc.summary.type.value_counts()


# In[7]:


ab_summary = pdata_ab_sc.summary

type_order = ['Control', 'IF']

# make color palette for Ab and PBS
colors = sns.color_palette("Paired")[2:4]
color_dict = dict(zip(type_order, colors))

# Plot
fig, ax = plt.subplots(figsize=(3, 3))
sns.violinplot(data=ab_summary, x="type", y="protein_count", ax=ax, alpha=1, inner='points', palette=color_dict)

# calculate p-value for significance of "control" vs "stained" protein_count
from scipy.stats import ttest_ind
control = ab_summary[ab_summary['type'] == 'Control']['protein_count']
stained = ab_summary[ab_summary['type'] == 'IF']['protein_count']

# scplt.plot_significance(ax, 2800, 30, pval=ttest_ind(control, stained, alternative='greater').pvalue, fontsize=11)


# Add significance
plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.ylim(0, 3200)
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_proteinCount.svg", bbox_inches="tight")


# In[22]:


ab_summary = pdata_ab_sc.summary

type_order = ['Control', 'IF']

# make color palette for Ab and PBS
colors = sns.color_palette("Paired")[2:4]
color_dict = dict(zip(type_order, colors))

# Plot
fig, ax = plt.subplots(figsize=(3, 3))
sns.barplot(data=ab_summary, x="type", y="protein_count", ax=ax, alpha=1, palette=color_dict, errorbar='sd', capsize=0.2, order=type_order)
sns.stripplot(x='type', y='protein_count', data=ab_summary, color='black', alpha=0.4, s=4.5, jitter=0.15, order=type_order)

# Add significance
plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.ylim(0, 3200)
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_proteinCount_bar.svg", bbox_inches="tight")


# In[61]:


fig, ax = plt.subplots(figsize=(3, 3))

scplt.plot_venn(ax, pdata_ab_sc, classes = 'type', set_colors = sns.color_palette("Paired")[2:4], label_order = ['Control', 'IF'])
plt.rcParams.update({'font.size': 12})
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_Venn.svg", bbox_inches="tight")

plt.rcParams.update({"font.size": 14})

fig, ax = plt.subplots(figsize=(3, 3))
ax, contents = scplt.plot_venn(
    ax, pdata_ab_sc,
    classes="type",
    set_colors=sns.color_palette("Paired")[2:4],
    label_order=["Control", "IF"],
    weighted=True,
    return_contents=True,
    fixed_subset_sizes=(1, 1, 3),
)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Control", "IF"}]:
    x, y = t.get_position()
    if t.get_text() == "Control":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")


plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-IF_sc_Venn_Weighted.svg", bbox_inches="tight")


# In[75]:


ab_upset = scutils.get_upset_contents(pdata_ab_sc, classes = ['type'])

pdata_ab_only = pdata_ab_sc.filter_sample(query_mode=True, values = 'type == "IF"', return_copy=True)
pdata_pbs_only = pdata_ab_sc.filter_sample(query_mode=True, values ='type == "Control"', return_copy=True)

ab_only = scutils.get_upset_query(ab_upset, present='IF', absent='Control')
pbs_only = scutils.get_upset_query(ab_upset, present='Control', absent='IF')

ab_all = scutils.get_upset_query(ab_upset, present='IF', absent=[])
pbs_all = scutils.get_upset_query(ab_upset, present='Control', absent=[])


# In[62]:


light_green_cmap = mcolors.LinearSegmentedColormap.from_list("light_green_cmap", ["#E8F5E9", "#388E3C"])
dark_green_cmap = mcolors.LinearSegmentedColormap.from_list("dark_green_cmap", ["#1B5E20", "#004D00"])

colors = sns.color_palette("Paired")[2:4]
order = ['Control', 'IF']
cmaps = [light_green_cmap,dark_green_cmap]
ylims = (2.5,9.5)

fig, ax = plt.subplots(figsize=(8,3), nrows = 1, ncols = 2, sharey='col', sharex='col', gridspec_kw={'width_ratios': [5, 3], 'hspace': 0.08, 'wspace': 0})
labels = order
line_width = 1

ax = plt.subplot(1, 2, 1)
scplt.plot_rankquant(ax, pdata_ab_sc, classes = ['type'], order = order, cmap = cmaps, color=colors, s=10, calpha=1, alpha=0.005)
plt.ylabel('MS Intensity')
ax.set_ylim(10**ylims[0],10**ylims[1])
legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
plt.legend(handles=legend_patches, bbox_to_anchor=(.75, 1), loc=2, borderaxespad=0., frameon=False)
# scplt.mark_rankquant(ax, pdata_ab_sc, mark_df=pbs_only, class_values=['Control'], show_label=False, color = 'black', alpha=0.3, label_type='gene')
# scplt.mark_rankquant(ax, pdata_ab_sc, mark_df=ab_only, class_values=['IF'], show_label=False, color = 'white', alpha=0.3, label_type='gene')

ax = plt.subplot(1, 2, 2)
scplt.plot_raincloud(ax, pdata_ab_sc, classes=['type'], order=order, color=colors, linewidth=line_width)
# scplt.mark_raincloud(ax, pdata_ab_sc, mark_df=pbs_only, class_values=['Control'], color = 'black', alpha=0.3)
# scplt.mark_raincloud(ax, pdata_ab_sc, mark_df=ab_only, class_values=['IF'], color = 'white', alpha=0.3, lowest_index = 1)
ax.set_yticks([],[])
ax.set_xticks(np.arange(len(order))+1, [case for case in order])
ax.set_ylim(ylims[0],ylims[1])

plt.subplots_adjust(wspace=0, hspace=0)
sns.despine()

# save as svg
plt.savefig('MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_rankquant.png', dpi=450, bbox_inches='tight')
# plt.show()


# In[79]:


# check df
ab_cc_df  = build_cc_df(ab_all,  protein_col='id', source_name="ab")
pbs_cc_df = build_cc_df(pbs_all, protein_col='id', source_name="pbs")

# counts 
ab_counts  = counts_from_cc_df(ab_cc_df).reindex(plot_order, fill_value=0)
pbs_counts = counts_from_cc_df(pbs_cc_df).reindex(plot_order, fill_value=0)

# convert to %
ab_pct  = ab_counts  / ab_counts.sum()  * 100
pbs_pct = pbs_counts / pbs_counts.sum() * 100

# order top→bottom like your plot
labels = plot_order[::-1]
display_labels = [label_map[l] for l in labels]
ab_pct  = ab_pct.reindex(labels)
pbs_pct = pbs_pct.reindex(labels)

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

# bar chart
fig, ax = plt.subplots(figsize=(2.7, 2.7))

y = np.arange(len(labels))
h = 0.35

colors = sns.color_palette("Paired")[2:4]

ax.barh(y - h/2, ab_pct.values, height=h, label="IF",
        color=colors[1], edgecolor="black", linewidth=0.6)
ax.barh(y + h/2, pbs_pct.values,  height=h, label="Control",
        color=colors[0], edgecolor="black", linewidth=0.6)

ax.set_yticks(y)
ax.set_yticklabels(display_labels)
ax.set_xlabel("% proteins")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", width=1.0, length=8)

ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_CC.svg", bbox_inches="tight")


# In[ ]:


# plot protein level of Ab vs PBS (log2)
pdata_p1 = pdata_ab_sc.copy()
pdata_p1.prot.layers['log2'] = np.log2(pdata_p1.prot.X.toarray() + 1)

# Get the positional indices for Ab and PBS samples
ab_pos_indices = [pdata_p1.prot.obs.index.get_loc(idx) for idx in pdata_p1.prot.obs.index[pdata_p1.prot.obs['type'] == 'IF']]
pbs_pos_indices = [pdata_p1.prot.obs.index.get_loc(idx) for idx in pdata_p1.prot.obs.index[pdata_p1.prot.obs['type'] == 'Control']]

log2_data = pdata_p1.prot.layers['log2']

ab_avg = np.mean(log2_data[ab_pos_indices, :], axis=0)
pbs_avg = np.mean(log2_data[pbs_pos_indices, :], axis=0)


valid_indices = ~np.isnan(ab_avg) & ~np.isnan(pbs_avg)
ab_avg = ab_avg[valid_indices]
pbs_avg = pbs_avg[valid_indices]

from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(ab_avg, pbs_avg)
line = slope * ab_avg + intercept

# Determine the range for the axes
min_val = 15  # Starting point as per your request
max_ab = np.max(ab_avg)
max_pbs = np.max(pbs_avg)
max_val = max(max_ab, max_pbs)
tick_step = 5

# Create a scatter plot
plt.figure(figsize=(3, 3))
plt.scatter(ab_avg, pbs_avg, c=ab_avg + pbs_avg, cmap='Greens', alpha=0.5)
plt.plot(ab_avg, line, color='black', label=f'$R^2$ = {r_value**2:.2f}')
plt.xlabel('Log2 Protein Abundance (IF)')
plt.ylabel('Log2 Protein Abundance (Control)')
# plt.title('Scatter Plot of Average Protein Abundance in Ab vs PBS Samples')
plt.legend()
plt.xlim(15-2, max_val + 0.5)
plt.ylim(15-2, max_val + 0.5)
ticks = np.arange(min_val, max_val + tick_step, tick_step)
plt.xticks(ticks)
plt.yticks(ticks)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_corr.svg", bbox_inches="tight")



# In[ ]:


from matplotlib.colors import LogNorm
from scipy.stats import linregress

# --- Prep data (log10) ---
pdata_p1 = pdata_ab_sc.copy()
pdata_p1.prot.layers["log10"] = np.log10(pdata_p1.prot.X.toarray() + 1)

if_idx = pdata_p1.prot.obs.index[pdata_p1.prot.obs["type"] == "IF"]
ct_idx = pdata_p1.prot.obs.index[pdata_p1.prot.obs["type"] == "Control"]

if_pos = [pdata_p1.prot.obs.index.get_loc(i) for i in if_idx]
ct_pos = [pdata_p1.prot.obs.index.get_loc(i) for i in ct_idx]

log10_data = pdata_p1.prot.layers["log10"]
x = np.mean(log10_data[if_pos, :], axis=0)       # IF
y = np.mean(log10_data[ct_pos, :], axis=0)       # Control

valid = ~np.isnan(x) & ~np.isnan(y)
x = x[valid]
y = y[valid]

res = linregress(x, y)
slope, intercept, r2 = res.slope, res.intercept, res.rvalue**2

lo = np.nanpercentile(np.r_[x, y], 0.5)
hi = np.nanpercentile(np.r_[x, y], 99.5)
pad = 0.03 * (hi - lo)
lo, hi = lo - pad, hi + pad

# density via 2D hist
bins = 120
H, xedges, yedges = np.histogram2d(x, y, bins=bins, range=[[lo, hi], [lo, hi]])
xi = np.clip(np.searchsorted(xedges, x, side="right") - 1, 0, bins - 1)
yi = np.clip(np.searchsorted(yedges, y, side="right") - 1, 0, bins - 1)
dens = H[xi, yi]

fig, ax = plt.subplots(figsize=(1.7, 2.52))
order = np.argsort(dens)  # low → high
x_s    = x[order]
y_s    = y[order]
dens_s = dens[order]

sc = ax.scatter(
    x_s, y_s,
    c=dens_s,
    s=9,
    cmap="YlGn",
    alpha=0.9,
    linewidths=0,
)

cb = plt.colorbar(sc, ax=ax, pad=0.02)
cb.set_label("Local density (proteins)", fontsize=8)
cb.ax.tick_params(labelsize=7)

# ax.plot([lo, hi], [lo, hi], "--", color="black", linewidth=1)
xs = np.array([lo, hi])
ax.plot(xs, slope * xs + intercept, color="black", linewidth=1)

ax.text(
    0.05, 0.95,
    f"$r^2$ = {r2:.2f}\nSlope = {slope:.2f}",
    transform=ax.transAxes,
    va="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.7, linewidth=0.8,edgecolor="grey"),
    fontsize=8,
)

ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel("$\log_{10}$ Protein Abundance (IF)")
ax.set_ylabel("$\log_{10}$ Protein Abundance (Control)")

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_corr_log10.svg", bbox_inches="tight")


# In[8]:


pdata_ab_sc_processed = pdata_ab_sc.copy()
pdata_ab_sc_processed = pdata_ab_sc_processed.filter_prot_found(min_ratio=0.4, group = ['type'], match_any=True)
pdata_ab_sc_processed = pdata_ab_sc_processed.filter_prot(valid_genes=True, unique_profiles=True)

pdata_ab_sc_norm = pdata_ab_sc_processed.copy()
# pdata_ab_sc_norm.impute(method='min', min_scale=0.2)
# pdata_ab_sc_norm.normalize(method='median')

# pdata_ab_sc_norm2 = pdata_ab_sc_processed.copy()
# pdata_ab_sc_norm.normalize(method='directlfq')


# In[48]:


fig, ax = plt.subplots(figsize=(2.7, 2.7))
ax, pca = scplt.plot_pca(ax, pdata_ab_sc_norm, classes='type', cmap=colors, return_fit=True, add_ellipses=False)
scplt.shift_legend(ax)
ax.legend_.remove()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_pca.svg", bbox_inches="tight")

fig, ax = plt.subplots(figsize=(2.7, 2.7))
ax, umap = scplt.plot_umap(ax, pdata_ab_sc_norm, classes='type', cmap=colors, umap_params={'min_dist': 0.1}, force=True, return_fit=True)
scplt.shift_legend(ax)
# remove legend
ax.legend_.remove()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_umap.svg", bbox_inches="tight")


# In[9]:


classes='type'
values=['IF','Control']

colors = sns.color_palette("Paired")[2:4]
color_dict = dict(zip(['downregulated', 'upregulated'], colors))

fig, ax = plt.subplots(1,1, figsize = (3,3))
scplt.plot_volcano(ax, pdata_ab_sc_norm, classes=classes, values=values, log2fc=1, pval=0.05, color=color_dict,
                   group_annot_kwargs={"pos": {"group1_xy": (0.98, 1.08), "group2_xy": (0.02, 1.08)}},)
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_IF_volcano.svg", bbox_inches="tight")
df_ab_norm = pdata_ab_sc_norm.stats["[{'type': 'IF'}] vs [{'type': 'Control'}]"]
df_ab_norm.to_csv('MANUSCRIPT/MANUSCRIPT_FigS1_DE-Ab.csv')


# In[ ]:


pdata_ab_sc_norm.list_enrichments()


# In[ ]:


pdata_ab_sc_norm.enrichment_functional(from_de=True, de_key="IF vs Control")


# In[ ]:


pdata_ab_sc_norm.plot_enrichment_svg("IF vs Control", direction="up")


# In[ ]:


pdata_ab_sc_norm.plot_enrichment_svg("IF vs Control", direction="down", category='RCTM')

# no enrichments


# #### HCR (hcr vs pbs)

# In[10]:


# obs_columns
pdata_hcr_sc = pdata_fig1_supplementary_sc.filter_sample(condition = 'type in ("GADHCR", "Control")')
pdata_hcr_sc = pdata_hcr_sc.filter_sample(min_prot = 1000)
pdata_hcr_sc = pdata_hcr_sc.filter_prot_significant()

mapping = {
    "GADHCR": "HCR",
}
pdata_hcr_sc.summary["type"] = pdata_hcr_sc.summary["type"].replace(mapping)
pdata_hcr_sc.update_summary()

pdata_hcr_sc.summary.type.value_counts()


# In[ ]:


hcr_summary = pdata_hcr_sc.summary

treatment_order = ['HCR', 'Control']

# make color palette for Ab and PBS
colors = sns.color_palette("Paired")[4:6]
color_dict = dict(zip(treatment_order, colors[::-1]))

# Plot
fig, ax = plt.subplots(figsize=(3, 3))
sns.violinplot(x='type', y='protein_count', ax=ax,
            data=hcr_summary, alpha=1, palette=color_dict, inner='points')

# calculate p-value for significance of "control" vs "stained" protein_count
from scipy.stats import ttest_ind
control = hcr_summary[hcr_summary['type'] == 'Control']['protein_count']
stained = hcr_summary[hcr_summary['type'] == 'HCR']['protein_count']

# scplt.plot_significance(ax, 3300, 100, pval=ttest_ind(control, stained).pvalue, fontsize=11)

# Add significance
plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 4000)
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_proteinCount.svg", bbox_inches="tight")


# In[25]:


hcr_summary = pdata_hcr_sc.summary

treatment_order = ['Control', 'HCR']

# make color palette for Ab and PBS
colors = sns.color_palette("Paired")[4:6]
color_dict = dict(zip(treatment_order, colors))

# Plot
fig, ax = plt.subplots(figsize=(3, 3))
sns.barplot(data=hcr_summary, x="type", y="protein_count", ax=ax, alpha=1, palette=color_dict, errorbar='sd', capsize=0.2, order=treatment_order)
sns.stripplot(x='type', y='protein_count', data=hcr_summary, color='black', alpha=0.4, s=4.5, jitter=0.15, order=treatment_order)

# Add significance
plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.ylim(0, 3200)
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_proteinCount_bar.svg", bbox_inches="tight")


# In[64]:


fig, ax = plt.subplots(figsize=(3, 3))

scplt.plot_venn(ax, pdata_hcr_sc, classes = 'type', set_colors = sns.color_palette("Paired")[4:6][::-1])
plt.rcParams.update({'font.size': 12})
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_Venn.svg", bbox_inches="tight")

plt.rcParams.update({"font.size": 14})

fig, ax = plt.subplots(figsize=(3, 3))
ax, contents = scplt.plot_venn(
    ax, pdata_hcr_sc,
    classes="type",
    set_colors=sns.color_palette("Paired")[4:6],
    label_order=["Control", "HCR"],
    weighted=True,
    return_contents=True,
    fixed_subset_sizes=(1, 1, 3),
)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Control", "HCR"}]:
    x, y = t.get_position()
    if t.get_text() == "Control":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")


plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-HCR_sc_Venn_Weighted.svg", bbox_inches="tight")


# In[80]:


ab_upset = scutils.get_upset_contents(pdata_hcr_sc, classes = ['type'])

pdata_HCR_only = pdata_hcr_sc.filter_sample(query_mode=True, values='type == "HCR"', return_copy=True)
pdata_control_only = pdata_hcr_sc.filter_sample(query_mode=True, values='type == "Control"', return_copy=True)


HCR_only = scutils.get_upset_query(ab_upset, present='HCR', absent='Control')
Control_only = scutils.get_upset_query(ab_upset, present='Control', absent='HCR')

hcr_all = scutils.get_upset_query(ab_upset, present='HCR', absent=[])
Control_all = scutils.get_upset_query(ab_upset, present='Control', absent=[])

# export ab_only and pbs_only to csv
# ab_only.to_csv('sc_HCR+_only.csv')
# pbs_only.to_csv('sc_HCR-_only.csv')


# In[81]:


# check df
hcr_cc_df  = build_cc_df(hcr_all,  protein_col='id', source_name="hcr")
pbs_cc_df = build_cc_df(Control_all, protein_col='id', source_name="pbs")

# counts 
hcr_counts  = counts_from_cc_df(hcr_cc_df).reindex(plot_order, fill_value=0)
pbs_counts = counts_from_cc_df(pbs_cc_df).reindex(plot_order, fill_value=0)

# convert to %
hcr_pct  = hcr_counts  / hcr_counts.sum()  * 100
pbs_pct = pbs_counts / pbs_counts.sum() * 100

# order top→bottom like your plot
labels = plot_order[::-1]
display_labels = [label_map[l] for l in labels]
hcr_pct  = hcr_pct.reindex(labels)
pbs_pct = pbs_pct.reindex(labels)

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

# bar chart
fig, ax = plt.subplots(figsize=(2.7, 2.7))

y = np.arange(len(labels))
h = 0.35

colors = sns.color_palette("Paired")[4:6]

ax.barh(y - h/2, hcr_pct.values, height=h, label="HCR",
        color=colors[1], edgecolor="black", linewidth=0.6)
ax.barh(y + h/2, pbs_pct.values,  height=h, label="Control",
        color=colors[0], edgecolor="black", linewidth=0.6)

ax.set_yticks(y)
ax.set_yticklabels(display_labels)
ax.set_xlabel("% proteins")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", width=1.0, length=8)

ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_CC.svg", bbox_inches="tight")


# In[ ]:


from matplotlib.colors import LinearSegmentedColormap

light_red_cmap = mcolors.LinearSegmentedColormap.from_list("light_red_cmap", ["#FFEBEE", "#D32F2F"])
dark_red_cmap = mcolors.LinearSegmentedColormap.from_list("dark_red_cmap", ["#B71C1C", "#4A0000"])

colors = sns.color_palette("Paired")[4:6][::-1]
order = ['HCR', 'Control']
cmaps = [light_red_cmap,dark_red_cmap]
ylims = (2,9.5)

fig, ax = plt.subplots(figsize=(8,3), nrows = 1, ncols = 2, sharey='col', sharex='col', gridspec_kw={'width_ratios': [5, 3], 'hspace': 0.08, 'wspace': 0})
labels = order
line_width = 1

ax = plt.subplot(1, 2, 1)
scplt.plot_rankquant(ax, pdata_hcr_sc, classes = ['type'], order = order, cmap = cmaps, color=colors, s=10, calpha=1, alpha=0.005)
plt.ylabel('Abundance')
ax.set_ylim(10**ylims[0],10**ylims[1])
legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
plt.legend(handles=legend_patches, bbox_to_anchor=(.75, 1), loc=2, borderaxespad=0., frameon=False)
# scplt.mark_rankquant(ax, pdata_hcr_sc, mark_df=pbs_only, class_values=['Control'], show_label=False, color = 'black', alpha=0.3, label_type='gene')
# scplt.mark_rankquant(ax, pdata_hcr_sc, mark_df=ab_only, class_values=['HCR'], show_label=False, color = 'white', alpha=0.3, label_type='gene')

ax = plt.subplot(1, 2, 2)
scplt.plot_raincloud(ax, pdata_hcr_sc, classes=['type'], order=order, color=colors, linewidth=line_width)
# scplt.mark_raincloud(ax, pdata_hcr_sc, mark_df=pbs_only, class_values=['Control'], color = 'black', alpha=0.3)
# scplt.mark_raincloud(ax, pdata_hcr_sc, mark_df=ab_only, class_values=['HCR'], color = 'white', alpha=0.3, lowest_index = 1)
ax.set_yticks([],[])
ax.set_xticks(np.arange(len(order))+1, [case for case in order])
ax.set_ylim(ylims[0],ylims[1])

plt.subplots_adjust(wspace=0, hspace=0)
sns.despine()

# save as svg
plt.savefig('MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_rankquant.png', dpi=450, bbox_inches='tight')
# plt.show()


# In[ ]:


# plot protein level of HCR vs Control (log2)
pdata_p1 = pdata_hcr_sc.copy()
pdata_p1.prot.layers['log2'] = np.log2(pdata_p1.prot.X.toarray() + 1)

# Get the positional indices for Ab and PBS samples
hcr_pos_indices = [pdata_p1.prot.obs.index.get_loc(idx) for idx in pdata_p1.prot.obs.index[pdata_p1.prot.obs['type'] == 'HCR']]
pbs_pos_indices = [pdata_p1.prot.obs.index.get_loc(idx) for idx in pdata_p1.prot.obs.index[pdata_p1.prot.obs['type'] == 'Control']]

log2_data = pdata_p1.prot.layers['log2']


hcr_avg = np.mean(log2_data[hcr_pos_indices, :], axis=0)
pbs_avg = np.mean(log2_data[pbs_pos_indices, :], axis=0)


valid_indices = ~np.isnan(hcr_avg) & ~np.isnan(pbs_avg)
hcr_avg = hcr_avg[valid_indices]
pbs_avg = pbs_avg[valid_indices]

from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(hcr_avg, pbs_avg)
line = slope * hcr_avg + intercept

# Determine the range for the axes
min_val = 15  # Starting point as per your request
max_hcr = np.max(hcr_avg)
max_pbs = np.max(pbs_avg)
max_val = max(max_hcr, max_pbs)
tick_step = 5

# Create a scatter plot
plt.figure(figsize=(3, 3))
plt.scatter(hcr_avg, pbs_avg, c=hcr_avg + pbs_avg, cmap='Reds', alpha=0.5)
plt.plot(hcr_avg, line, color='black', label=f'$R^2$ = {r_value**2:.2f}')
plt.xlabel('Log2 Protein Abundance (HCR)')
plt.ylabel('Log2 Protein Abundance (Control)')
# plt.title('Scatter Plot of Average Protein Abundance in Ab vs PBS Samples')
plt.legend()
plt.xlim(15-4, max_val + 0.5)
plt.ylim(15-4, max_val + 0.5)
ticks = np.arange(min_val, max_val + tick_step, tick_step)
plt.xticks(ticks)
plt.yticks(ticks)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_corr.svg", bbox_inches="tight")



# In[ ]:


# --- Prep data (log10) ---
pdata_p1 = pdata_hcr_sc.copy()
pdata_p1.prot.layers["log10"] = np.log10(pdata_p1.prot.X.toarray() + 1)

hcr_idx = pdata_p1.prot.obs.index[pdata_p1.prot.obs["type"] == "HCR"]
ct_idx  = pdata_p1.prot.obs.index[pdata_p1.prot.obs["type"] == "Control"]

hcr_pos = [pdata_p1.prot.obs.index.get_loc(i) for i in hcr_idx]
ct_pos  = [pdata_p1.prot.obs.index.get_loc(i) for i in ct_idx]

log10_data = pdata_p1.prot.layers["log10"]
x = np.mean(log10_data[hcr_pos, :], axis=0)  # HCR
y = np.mean(log10_data[ct_pos,  :], axis=0)  # Control

valid = ~np.isnan(x) & ~np.isnan(y)
x = x[valid]
y = y[valid]

# Regression
res = linregress(x, y)
slope, intercept, r2 = res.slope, res.intercept, res.rvalue**2

# Limits
lo = np.nanpercentile(np.r_[x, y], 0.5)
hi = np.nanpercentile(np.r_[x, y], 99.5)
pad = 0.03 * (hi - lo)
lo, hi = lo - pad, hi + pad

# --- Density via 2D histogram ---
bins = 120
H, xedges, yedges = np.histogram2d(x, y, bins=bins, range=[[lo, hi], [lo, hi]])
xi = np.clip(np.searchsorted(xedges, x, side="right") - 1, 0, bins - 1)
yi = np.clip(np.searchsorted(yedges, y, side="right") - 1, 0, bins - 1)
dens = H[xi, yi]

# Sort so dense points are drawn last
order = np.argsort(dens)
x_s, y_s, dens_s = x[order], y[order], dens[order]

fig, ax = plt.subplots(figsize=(1.7, 2.52))

sc = ax.scatter(
    x_s, y_s,
    c=dens_s,
    s=9,
    cmap="Reds",
    alpha=0.9,
    linewidths=0,
)

cb = plt.colorbar(sc, ax=ax, pad=0.02)
cb.set_label("Local density (proteins)", fontsize=8)
cb.ax.tick_params(labelsize=7)

# Regression line
xs = np.array([lo, hi])
ax.plot(xs, slope * xs + intercept, color="black", linewidth=1)

ax.text(
    0.05, 0.95,
    f"$r^2$ = {r2:.2f}\nSlope = {slope:.2f}",
    transform=ax.transAxes,
    va="top",
    bbox=dict(
        boxstyle="round",
        facecolor="white",
        alpha=0.7,
        linewidth=0.8,
        edgecolor="grey",
    ),
    fontsize=8,
)

ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel("log10 Protein Abundance (HCR)")
ax.set_ylabel("log10 Protein Abundance (Control)")

plt.savefig(
    "MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_corr_log10.svg",
    bbox_inches="tight",
)


# In[11]:


pdata_hcr_sc_processed = pdata_hcr_sc.copy()
pdata_hcr_sc_processed = pdata_hcr_sc_processed.filter_prot_found(min_ratio=0.4, group = ['type'], match_any=True)
pdata_hcr_sc_processed = pdata_hcr_sc_processed.filter_prot(valid_genes=True, unique_profiles=True)

pdata_hcr_sc_norm = pdata_hcr_sc_processed.copy()
# pdata_ab_sc_norm.impute(method='min', min_scale=0.2)
# pdata_ab_sc_norm.normalize(method='median')

# pdata_ab_sc_norm2 = pdata_ab_sc_processed.copy()
# pdata_ab_sc_norm.normalize(method='directlfq')


# In[52]:


fig, ax = plt.subplots(figsize=(2.7, 2.7))
ax, pca = scplt.plot_pca(ax, pdata_hcr_sc_norm, classes='type', cmap=colors[::-1], return_fit=True, add_ellipses=False)
scplt.shift_legend(ax)
ax.legend_.remove()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_pca.svg", bbox_inches="tight")

fig, ax = plt.subplots(figsize=(2.7, 2.7))
ax, umap = scplt.plot_umap(ax, pdata_hcr_sc_norm, classes='type', cmap=colors[::-1], umap_params={'min_dist': 0.1}, force=True, return_fit=True)
scplt.shift_legend(ax)
# remove legend
ax.legend_.remove()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_umap.svg", bbox_inches="tight")


# In[12]:


classes='type'
values=['HCR','Control']

colors = sns.color_palette("Paired")[4:6]
color_dict = dict(zip(['downregulated', 'upregulated'], colors))

fig, ax = plt.subplots(1,1, figsize = (3.5,3.5))
scplt.plot_volcano(ax, pdata_hcr_sc_norm, classes=classes, values=values, log2fc=1, pval=0.05, color=color_dict)
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-sc_HCR_volcano.svg", bbox_inches="tight")
df_hcr_norm = pdata_hcr_sc_norm.stats["[{'type': 'HCR'}] vs [{'type': 'Control'}]"]
df_hcr_norm.to_csv('MANUSCRIPT/MANUSCRIPT_FigS1_DE-HCR.csv')


# In[ ]:


pdata_hcr_sc_norm.list_enrichments()


# In[ ]:


pdata_hcr_sc_norm.enrichment_functional(from_de=True, de_key="HCR vs Control")


# In[ ]:


pdata_hcr_sc_norm.plot_enrichment_svg("HCR vs Control", direction="up", category="Component")


# In[ ]:


pdata_hcr_sc_norm.plot_enrichment_svg("HCR vs Control", direction="down")


# ## fixed vs frozen

# In[15]:


import re

# (plasma membrane separate; organelle membranes go with organelle)
selected_components = [
    "cytoplasm",
    "plasma membrane",
    "extracellular region",
    "mitochondrion",
    "nucleus",
    "endoplasmic reticulum",
    "golgi apparatus",
    "cytoskeleton",
]

label_map = {
    "cytoplasm": "Cytoplasm",
    "plasma membrane": "Plasma membrane",
    "extracellular region": "Extracellular region",
    "mitochondrion": "Mitochondria",
    "nucleus": "Nucleus",
    "endoplasmic reticulum": "ER",
    "golgi apparatus": "Golgi apparatus",
    "cytoskeleton": "Cytoskeleton",
    "others": "Others",
}

plot_order = selected_components + ["others"]

def _clean_cc_term(s: str) -> str:
    if pd.isna(s):
        return ""
    s = str(s).strip()
    s = re.sub(r"\[.*?\]", "", s).strip()
    return s

def map_cc_term(term: str) -> str:
    t = term.lower()

    # cytoskeleton (keep separate)
    if "cytoskeleton" in t:
        return "cytoskeleton"

    # nucleus
    if re.search(r"\bnucleus\b|nucleoplasm|nuclear", t):
        return "nucleus"

    # cytoplasm (cytosol/cytoplasmic included)
    if re.search(r"cytosol|cytoplasm|cytoplasmic", t):
        return "cytoplasm"

    # extracellular
    if "extracellular" in t:
        return "extracellular region"

    # ER (includes ER membrane etc.)
    if "endoplasmic reticulum" in t:
        return "endoplasmic reticulum"

    # Golgi
    if "golgi" in t:
        return "golgi apparatus"

    # mitochondrion (includes mitochondrial membrane etc.)
    if "mitochondri" in t:
        return "mitochondrion"

    # plasma membrane (keep separate from generic "membrane")
    if "plasma membrane" in t:
        return "plasma membrane"

    # everything else
    return "others"

def build_cc_df(df, cc_col="gene_ontology_cellular_component", protein_col=None, source_name=None):
    prot = df.index.to_series() if protein_col is None else df[protein_col]
    out = df[[cc_col]].copy()
    out["protein"] = prot.values
    out = out.rename(columns={cc_col: "cc_raw"})
    if source_name is not None:
        out["source"] = source_name

    out["cc_raw"] = out["cc_raw"].astype(str).str.split(";")
    out = out.explode("cc_raw", ignore_index=True)

    out["cc_clean"] = out["cc_raw"].map(_clean_cc_term)
    out = out[out["cc_clean"].ne("") & out["cc_clean"].ne("nan")].copy()

    out["cc_mapped"] = out["cc_clean"].map(map_cc_term)
    out["in_selected"] = out["cc_mapped"].isin(selected_components)
    return out

def counts_from_cc_df(cc_df):
    return cc_df["cc_mapped"].value_counts()


# In[16]:


obs_columns = ['sample', 'amount', 'type']

pdata_ff = pAnnData.import_proteomeDiscoverer(prot_file = 'data/MouseFF_100ng_FixedFrozen_43minDDA_20230623_Proteins.txt', pep_file='data/MouseFF_100ng_FixedFrozen_43minDDA_20230623_PeptideGroups.txt', obs_columns=obs_columns)


# In[17]:


pdata_ff.summary['type'] = pdata_ff.summary['type'].replace({'fixed': 'Fixed', 'frozen': 'Frozen'})
pdata_ff.update_summary()


# In[18]:


ff_summary = pdata_ff.summary

treatment_order = ['Frozen', 'Fixed']

# make color palette for Ab and PBS
colors = sns.color_palette("Paired")[0:2]
color_dict = dict(zip(treatment_order, colors))

# Plot
fig, ax = plt.subplots(figsize=(1.74, 2.13))
sns.swarmplot(data=ff_summary, x="type", y="protein_count", ax=ax, color='k', order=treatment_order)
sns.barplot(data=ff_summary, x="type", y="protein_count", ax=ax, errorbar = 'ci', alpha=1, palette=color_dict, order=treatment_order)

# calculate p-value for significance of "control" vs "stained" protein_count
from scipy.stats import ttest_ind
control = ff_summary[ff_summary['type'] == 'Frozen']['protein_count']
stained = ff_summary[ff_summary['type'] == 'Fixed']['protein_count']

# scplt.plot_significance(ax, 2200, 50, pval=ttest_ind(control, stained).pvalue, fontsize=11)

# Add significance
plt.ylabel('Protein Count')
plt.xlabel('')
# plt.ylabel('Protein Count',fontsize=13, labelpad=4)
# plt.xlabel('', fontsize=13, labelpad=4)
# ax.tick_params(axis='x', labelsize=11)
# ax.tick_params(axis='y', labelsize=11)

plt.ylim(1000, 2400)
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-FF_proteinCount.svg", bbox_inches="tight")
# plt.show()


# In[7]:


treatment_order = ['Frozen', 'Fixed']

fig, ax = plt.subplots(figsize=(3, 3))
scplt.plot_venn(ax, pdata_ff, classes = 'type', set_colors = sns.color_palette("Paired")[0:2], label_order=treatment_order)
plt.rcParams.update({'font.size': 12})
# plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-FF_Venn.svg", bbox_inches="tight")
plt.show()


# In[36]:


import matplotlib.patheffects as pe

treatment_order = ['Frozen', 'Fixed']

fig, ax = plt.subplots(figsize=(3, 3))
ax = scplt.plot_venn(ax, pdata_ff, classes = 'type', set_colors = sns.color_palette("Paired")[0:2], label_order=treatment_order, weighted=True, fixed_subset_sizes=(1, 1, 3),)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Frozen", "Fixed"}]:
    x, y = t.get_position()
    if t.get_text() == "Frozen":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")

plt.rcParams.update({'font.size': 12})
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-FF_Venn_Weighted.svg", bbox_inches="tight")
# plt.show()


# In[9]:


# found in all 3 replicates
pdata_ff_filtered = pdata_ff.copy()
pdata_ff_filtered = pdata_ff_filtered.filter_prot_found('type', min_ratio=1, match_any=True)

fig, ax = plt.subplots(figsize=(3, 3))
scplt.plot_venn(ax, pdata_ff_filtered, classes = 'type', set_colors = sns.color_palette("Paired")[0:2], label_order=treatment_order)
plt.rcParams.update({'font.size': 12})
plt.show()


# In[37]:


ff_upset = scutils.get_upset_contents(pdata_ff, classes = ['type'])

pdata_fixed_only = pdata_ff.filter_sample(condition = 'type == "Fixed"')
pdata_frozen_only = pdata_ff.filter_sample(condition = 'type == "Frozen"')

fixed_only = scutils.get_upset_query(ff_upset, present='Fixed', absent='Frozen')
frozen_only = scutils.get_upset_query(ff_upset, present='Frozen', absent='Fixed')

fixed_all = scutils.get_upset_query(ff_upset, present='Fixed', absent=[])
frozen_all = scutils.get_upset_query(ff_upset, present='Frozen', absent=[])

# export ab_only and pbs_only to csv
# fixed_only.to_csv('fixed_only.csv')
# frozen_only.to_csv('frozen_only.csv')


# In[85]:


# check df
fixed_cc_df  = build_cc_df(fixed_all,  protein_col=None, source_name="fixed")
frozen_cc_df = build_cc_df(frozen_all, protein_col=None, source_name="frozen")

# counts 
fixed_counts  = counts_from_cc_df(fixed_cc_df).reindex(plot_order, fill_value=0)
frozen_counts = counts_from_cc_df(frozen_cc_df).reindex(plot_order, fill_value=0)

# convert to %
fixed_pct  = fixed_counts  / fixed_counts.sum()  * 100
frozen_pct = frozen_counts / frozen_counts.sum() * 100

# order top→bottom like your plot
labels = plot_order[::-1]
display_labels = [label_map[l] for l in labels]
fixed_pct  = fixed_pct.reindex(labels)
frozen_pct = frozen_pct.reindex(labels)

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})


# bar chart
fig, ax = plt.subplots(figsize=(2.7, 2.7))

y = np.arange(len(labels))
h = 0.35

ax.barh(y - h/2, frozen_pct.values, height=h, label="Frozen ",
        color=(0.80, 0.88, 0.98), edgecolor="black", linewidth=0.6)
ax.barh(y + h/2, fixed_pct.values,  height=h, label="Fixed ",
        color=(0.00, 0.45, 0.70), edgecolor="black", linewidth=0.6)

ax.set_yticks(y)
ax.set_yticklabels(display_labels)
ax.set_xlabel("% proteins")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", width=1.0, length=8)

ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-FF_CC.svg", bbox_inches="tight")


# In[39]:


from matplotlib.colors import LinearSegmentedColormap

light_blue_cmap = mcolors.LinearSegmentedColormap.from_list("light_blue_cmap", ["#E0F7FA", "#0277BD"])
dark_blue_cmap = mcolors.LinearSegmentedColormap.from_list("dark_blue_cmap", ["#01579B", "#002171"])

colors = sns.color_palette("Paired")[0:2]
order = ['Frozen', 'Fixed']
cmaps = [light_blue_cmap,dark_blue_cmap]
ylims = (3,10.5)

fig, ax = plt.subplots(figsize=(5.62,2.13), nrows = 1, ncols = 2, sharey='col', sharex='col', gridspec_kw={'width_ratios': [11.34, 6.85], 'hspace': 0.08, 'wspace': 0})
labels = order
line_width = 1

ax = plt.subplot(1, 2, 1)
scplt.plot_rankquant(ax, pdata_ff, classes = ['type'], order = order, cmap = cmaps, color=colors, s=10, calpha=1, alpha=0.005)
plt.ylabel('MS Intensity')
ax.set_ylim(10**ylims[0],10**ylims[1])
legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
plt.legend(handles=legend_patches, bbox_to_anchor=(.6, 1), loc=2, borderaxespad=0., frameon=False)
scplt.mark_rankquant(ax, pdata_ff, mark_df=frozen_only, class_values=['Frozen'], show_label=False, color = 'black', alpha=0.3, label_type='gene')
scplt.mark_rankquant(ax, pdata_ff, mark_df=fixed_only, class_values=['Fixed'], show_label=False, color = 'white', alpha=0.3, label_type='gene')

ax = plt.subplot(1, 2, 2)
scplt.plot_raincloud(ax, pdata_ff, classes=['type'], order=order, color=colors, linewidth=line_width)
scplt.mark_raincloud(ax, pdata_ff, mark_df=frozen_only, class_values=['Frozen'], color = 'black', alpha=0.3)
scplt.mark_raincloud(ax, pdata_ff, mark_df=fixed_only, class_values=['Fixed'], color = 'white', alpha=0.3, lowest_index = 1)
ax.set_yticks([],[])
ax.set_xticks(np.arange(len(order))+1, [case for case in order])
ax.set_ylim(ylims[0],ylims[1])

plt.subplots_adjust(wspace=0, hspace=0)
sns.despine()

# save as svg
plt.savefig('MANUSCRIPT/MANUSCRIPT_FigS1-FF_rankquant.png', dpi=450, bbox_inches='tight')
# plt.show()


# In[ ]:


# Get the positional indices for Ab and PBS samples
from math import log2

pdata_p1 = pdata_ff.copy()
pdata_p1.prot.layers['log2'] = np.log2(pdata_p1.prot.X.toarray() + 1)

ab_pos_indices = [pdata_p1.prot.obs.index.get_loc(idx) for idx in pdata_p1.prot.obs.index[pdata_p1.prot.obs['type'] == 'Fixed']]
pbs_pos_indices = [pdata_p1.prot.obs.index.get_loc(idx) for idx in pdata_p1.prot.obs.index[pdata_p1.prot.obs['type'] == 'Frozen']]

log2_data = pdata_p1.prot.layers['log2']

ab_avg = np.mean(log2_data[ab_pos_indices, :], axis=0)
pbs_avg = np.mean(log2_data[pbs_pos_indices, :], axis=0)

valid_indices = ~np.isnan(ab_avg) & ~np.isnan(pbs_avg)
ab_avg = ab_avg[valid_indices]
pbs_avg = pbs_avg[valid_indices]

from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(ab_avg, pbs_avg)
line = slope * ab_avg + intercept

# Determine the range for the axes
min_val = 10  # Starting point as per your request
max_ab = np.max(ab_avg)
max_pbs = np.max(pbs_avg)
max_val = max(max_ab, max_pbs)
tick_step = 5

# Create a scatter plot
plt.figure(figsize=(3, 3))
plt.scatter(ab_avg, pbs_avg, c=ab_avg + pbs_avg, cmap='Blues', alpha=0.5)
plt.plot(ab_avg, line, color='black', label=f'$R^2$ = {r_value**2:.2f}')
plt.xlabel('Log2 Protein Abundance (Fixed)')
plt.ylabel('Log2 Protein Abundance (Frozen)')
# plt.title('Scatter Plot of Average Protein Abundance in Ab vs PBS Samples')
plt.legend()
plt.xlim(15, max_val + 0.5)
plt.ylim(15, max_val + 0.5)
ticks = np.arange(min_val, max_val + tick_step, tick_step)
plt.xticks(ticks)
plt.yticks(ticks)

plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-FF_corr.svg", bbox_inches="tight")


# In[ ]:


# --- Prep data (log10) ---
pdata_p1 = pdata_ff.copy()
pdata_p1.prot.layers["log10"] = np.log10(pdata_p1.prot.X.toarray() + 1)

fixed_idx = pdata_p1.prot.obs.index[pdata_p1.prot.obs["type"] == "Fixed"]
frozen_idx  = pdata_p1.prot.obs.index[pdata_p1.prot.obs["type"] == "Frozen"]

fixed_pos = [pdata_p1.prot.obs.index.get_loc(i) for i in fixed_idx]
frozen_pos  = [pdata_p1.prot.obs.index.get_loc(i) for i in frozen_idx]

log10_data = pdata_p1.prot.layers["log10"]
x = np.mean(log10_data[fixed_pos, :], axis=0)  # Fixed
y = np.mean(log10_data[frozen_pos,  :], axis=0)  # Frozen

valid = ~np.isnan(x) & ~np.isnan(y)
x = x[valid]
y = y[valid]

# Regression
res = linregress(x, y)
slope, intercept, r2 = res.slope, res.intercept, res.rvalue**2

# Limits
lo = np.nanpercentile(np.r_[x, y], 0.5)
hi = np.nanpercentile(np.r_[x, y], 99.5)
pad = 0.03 * (hi - lo)
lo, hi = lo - pad, hi + pad

# --- Density via 2D histogram ---
bins = 120
H, xedges, yedges = np.histogram2d(x, y, bins=bins, range=[[lo, hi], [lo, hi]])
xi = np.clip(np.searchsorted(xedges, x, side="right") - 1, 0, bins - 1)
yi = np.clip(np.searchsorted(yedges, y, side="right") - 1, 0, bins - 1)
dens = H[xi, yi]

# Sort so dense points are drawn last
order = np.argsort(dens)
x_s, y_s, dens_s = x[order], y[order], dens[order]

fig, ax = plt.subplots(figsize=(2.6, 3))

sc = ax.scatter(
    x_s, y_s,
    c=dens_s,
    s=9,
    cmap="Blues",
    alpha=0.9,
    linewidths=0,
)

cb = plt.colorbar(sc, ax=ax, pad=0.02)
cb.set_label("Local density (proteins)", fontsize=8)
cb.ax.tick_params(labelsize=7)

# Regression line
xs = np.array([lo, hi])
ax.plot(xs, slope * xs + intercept, color="black", linewidth=1)

ax.text(
    0.05, 0.95,
    f"$r^2$ = {r2:.2f}\nSlope = {slope:.2f}",
    transform=ax.transAxes,
    va="top",
    bbox=dict(
        boxstyle="round",
        facecolor="white",
        alpha=0.7,
        linewidth=0.8,
        edgecolor="grey",
    ),
    fontsize=8,
)

ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel("log10 Protein Abundance (HCR)")
ax.set_ylabel("log10 Protein Abundance (Control)")

plt.savefig(
    "MANUSCRIPT/MANUSCRIPT_FigS1-FF_corr.svg",
    bbox_inches="tight",
)


# ### Mouse proteome

# In[5]:


df_mouse_proteome = pd.read_csv("MANUSCRIPT/uniprotkb_proteome_UP000000589_AND_revi_2026_01_19.tsv", delimiter='\t')
df_mouse_proteome


# In[10]:


# check df
fixed_cc_df  = build_cc_df(fixed_all,  protein_col=None, source_name="fixed")
frozen_cc_df = build_cc_df(frozen_all, protein_col=None, source_name="frozen")
mouse_cc_df  = build_cc_df(df_mouse_proteome, protein_col="Entry Name",
                           cc_col="Gene Ontology (cellular component)", source_name="mouse")

# counts
fixed_counts  = counts_from_cc_df(fixed_cc_df).reindex(plot_order, fill_value=0)
frozen_counts = counts_from_cc_df(frozen_cc_df).reindex(plot_order, fill_value=0)
mouse_counts  = counts_from_cc_df(mouse_cc_df).reindex(plot_order, fill_value=0)

# convert to %
fixed_pct  = fixed_counts  / fixed_counts.sum()  * 100
frozen_pct = frozen_counts / frozen_counts.sum() * 100
mouse_pct  = mouse_counts  / mouse_counts.sum()  * 100

# order top→bottom like your plot
labels = plot_order[::-1]
display_labels = [label_map[l] for l in labels]
fixed_pct  = fixed_pct.reindex(labels)
frozen_pct = frozen_pct.reindex(labels)
mouse_pct  = mouse_pct.reindex(labels)

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

# bar chart
fig, ax = plt.subplots(figsize=(3.2, 2.7))  # slightly wider for 3 bars

y = np.arange(len(labels))
h = 0.24  # 3 bars -> smaller height

ax.barh(y - h, frozen_pct.values, height=h, label="Frozen",
        color=(0.80, 0.88, 0.98), edgecolor="black", linewidth=0.6)
ax.barh(y,      fixed_pct.values,  height=h, label="Fixed",
        color=(0.00, 0.45, 0.70), edgecolor="black", linewidth=0.6)
ax.barh(y + h,  mouse_pct.values,  height=h, label="Mouse proteome",
        color=(0.70, 0.70, 0.70), edgecolor="black", linewidth=0.6)

ax.set_yticks(y)
ax.set_yticklabels(display_labels)
ax.set_xlabel("% proteins")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", width=1.0, length=8)

ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
plt.show()
# plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS1-FF_CC.svg", bbox_inches="tight")


# # Fig 2 - region

# ## Cortex vs snpc

# ### pre-processing

# In[3]:


pdata_fig2_region_filter = pdata_all_region_filtered.filter_prot_significant()
pdata_fig2_region_filter = pdata_fig2_region_filter.filter_sample(condition="Grouping != 'ast'")
pdata_fig2_region_filter = pdata_fig2_region_filter.filter_sample(min_prot=1000)

pdata_fig2_region_filter.summary["region"] = pdata_fig2_region_filter.summary["Grouping"].replace({
    "snpc": "SNpc",
    "ctx": "Cortex"
})
pdata_fig2_region_filter.update_summary()

pdata_fig2_region_processed = pdata_fig2_region_filter.copy()
pdata_fig2_region_processed = pdata_fig2_region_processed.filter_prot_found(min_ratio=0.4, group = ['Grouping'], match_any=True)
pdata_fig2_region_processed = pdata_fig2_region_processed.filter_prot(valid_genes=True, unique_profiles=True)


# In[ ]:


cns_color={
    "Cortex": "#D19DCB",
    "SNpc": "#85BE9E"
}

text_kwargs=dict(fontsize=11,color='black',
    offset=1,         # vertical offset above anchor
)

bar_kwargs=dict(width=0.2, edgecolor="black", linewidth=0.6,)

strip_kwargs=dict(darken_factor=0.55, alpha=0.8)

fig, ax, df = pdata_fig2_region_processed.plot_abundance_boxgrid(namelist=['Th'], classes="region", fig_width=1.7, fig_height=3.18, 
    label_x=True, global_legend=False, plot_type="bar", show_n=False, text_kwargs=text_kwargs, bar_kwargs=bar_kwargs, strip_kwargs=strip_kwargs,
    palette=cns_color, order=['Cortex','SNpc'],log_scale=False, return_df=True)

# ymin, ymax = ax[0].get_ylim()
ax[0].set_ylim(0, 400000)

df


# In[ ]:


pdata_fig2_region_processed.get_abundance(namelist=['Th','Sncg','Camk2a']).to_csv('MANUSCRIPT/MANUSCRIPT_Fig2_Abundance-Th-Sncg-Camk2a.csv')


# In[ ]:


pdata_fig2_region_norm = pdata_fig2_region_processed.copy()
pdata_fig2_region_norm.normalize(method='directlfq')

pdata_fig2_region_norm_impute = pdata_fig2_region_processed.copy()
pdata_fig2_region_norm_impute.impute(method='min', min_scale=0.2)
pdata_fig2_region_norm_impute.normalize(method='median')


# In[ ]:


pdata_fig2_region_norm.summary.region.value_counts()


# ### plots

# In[4]:


cns_color={
    "Cortex": "#D19DCB",
    "SNpc": "#85BE9E"
}


# In[ ]:


order = ['Cortex', 'SNpc']

fig, ax = plt.subplots(figsize=(2, 3))
sns.barplot(x='region', y='protein_count',
            data=pdata_fig2_region_filter.summary, errorbar='sd', capsize=.1,
            saturation=1, alpha=0.5, palette=cns_color,width=0.75,
            order=order, ax=ax, edgecolor='black',      # outline color
            linewidth=0.8)
sns.swarmplot(x='region', y='protein_count',
                data=pdata_fig2_region_filter.summary, color='k', ax=ax)

from scipy.stats import ttest_ind
group_cortex = pdata_fig2_region_filter.summary[pdata_fig2_region_filter.summary['region'] == 'Cortex']['protein_count']
group_snpc = pdata_fig2_region_filter.summary[pdata_fig2_region_filter.summary['region'] == 'SNpc']['protein_count']
max_prot_count = pdata_fig2_region_filter.summary['protein_count'].max()

scplt.plot_significance(ax, max_prot_count+350, 100,
                        pval=ttest_ind(group_cortex, group_snpc).pvalue,
                        fontsize=11)

plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('Region', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 4200)
# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig2-CNS_protein-count.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


order = ['Cortex', 'SNpc']

fig, ax = plt.subplots(figsize=(1.76, 3))
sns.violinplot(x='region', y='protein_count',
            data=pdata_fig2_region_filter.summary,
            saturation=1, alpha=0.5, palette=cns_color,width=0.75,
            order=order, ax=ax, edgecolor='black', inner="points",     # outline color
            linewidth=0.8)
# sns.swarmplot(x='region', y='protein_count',
#                 data=pdata_fig2_region_filter.summary, color='k', ax=ax,)

from scipy.stats import ttest_ind
group_cortex = pdata_fig2_region_filter.summary[pdata_fig2_region_filter.summary['region'] == 'Cortex']['protein_count']
group_snpc = pdata_fig2_region_filter.summary[pdata_fig2_region_filter.summary['region'] == 'SNpc']['protein_count']
max_prot_count = pdata_fig2_region_filter.summary['protein_count'].max()

# scplt.plot_significance(ax, max_prot_count+550, 100,
#                         pval=ttest_ind(group_cortex, group_snpc).pvalue,
#                         fontsize=11)

plt.ylabel('Protein Count')
plt.xlabel('')
# plt.ylabel('Protein Count',fontsize=13, labelpad=4)
# plt.xlabel('Region', fontsize=13, labelpad=4)
# ax.tick_params(axis='x', labelsize=11)
# ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 4200)
# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig2-CNS_protein-count-violin-nosig.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


# CV

fig, ax = plt.subplots(figsize=(2.21,3))
cv_df = scplt.plot_cv(ax, pdata_fig2_region_filter, classes = 'region', return_df=True)

cv_df = cv_df.reset_index()
sns.violinplot(data=cv_df, x='Class', y='CV', palette=cns_color, linewidth=1, inner='quartile', saturation = 1, ax=ax)
# ax.set_title(classes)
# legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
# plt.legend(handles=legend_patches, bbox_to_anchor=(.75, 1), loc=2, borderaxespad=0., frameon=False)

from matplotlib.ticker import PercentFormatter
ax.yaxis.set_major_formatter(PercentFormatter(1))
ax.set_xlabel('')
ax.set_ylabel('Protein Abundance CV')

# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig2-CNS_CV.svg", bbox_inches="tight", pad_inches=0.1)


# In[5]:


fig, ax = plt.subplots(figsize=(3, 3))
plt.rcParams.update({'font.size': 12})
ax, venn_contents = scplt.plot_venn(ax, pdata_fig2_region_filter,
                                    classes='region',
                                    set_colors=["#D19DCB", '#85BE9E'],
                                    return_contents=True,
                                    weighted=True, fixed_subset_sizes=(1, 1, 3),)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Cortex", "SNpC"}]:
    x, y = t.get_position()
    if t.get_text() == "Cortex":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")

plt.show()
# plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig2-CNS_venn.svg", bbox_inches="tight", pad_inches=0.1)

setCortex = set(venn_contents["Cortex"])
setSNpc = set(venn_contents["SNpc"])
venn_df = pd.DataFrame({
    "cortex_only": pd.Series(sorted(setCortex - setSNpc)),
    "snpc_only": pd.Series(sorted(setSNpc - setCortex)),
    "Both":      pd.Series(sorted(setCortex & setSNpc))
})


# In[ ]:


cortex_list = venn_df['cortex_only'].unique().tolist()[:-1]
cortex_only = scutils.get_uniprot_fields(cortex_list)

cortex_only


# In[18]:


cortex_only['gene_primary'].to_list()


# #### directlfq

# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
ax, pca = scplt.plot_pca(ax, pdata_fig2_region_norm, classes='region', s=20, alpha=.8, return_fit=True, force=True, cmap=cns_color, add_ellipses=False)
# scplt.shift_legend(ax)

legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')
if legend is not None:
    legend.remove()

ax.set_xlabel(ax.get_xlabel(), fontsize=13, labelpad=4)
ax.set_ylabel(ax.get_ylabel(), fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_PCA.svg", bbox_inches="tight")
# plt.show()

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111, projection='3d')
ax = scplt.plot_pca(ax, pdata_fig2_region_norm, classes=['region'], force=True, plot_pc = [1,2,3], cmap=cns_color)
zlab = ax.get_zlabel()
ax.set_zlabel("")
ax.text2D(1.05, 0.55, zlab,
          transform=ax.transAxes, rotation=90,
          ha='left', va='center', fontsize=9.5)


# plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_PCA-3D.svg", bbox_inches="tight", pad_inches=0.1)
plt.show()


# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
umap_params = {'n_neighbors': 20, 'min_dist': 0.1}
ax, umap = scplt.plot_umap(ax, pdata_fig2_region_norm, classes = 'region', s=20, alpha=.8, force = True, umap_params=umap_params, return_fit=True, cmap=cns_color)
scplt.shift_legend(ax)

legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')

# plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_UMAP.svg", bbox_inches="tight", pad_inches=0.1)
plt.show()


# #### impute

# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
ax, pca = scplt.plot_pca(ax, pdata_fig2_region_norm_impute, classes='region', s=20, alpha=.8, return_fit=True, force=True, cmap=cns_color, add_ellipses=False)
scplt.shift_legend(ax)

legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')
if legend is not None:
    legend.remove()

ax.set_xlabel(ax.get_xlabel(), fontsize=13, labelpad=4)
ax.set_ylabel(ax.get_ylabel(), fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_PCA-impute.svg", bbox_inches="tight")
# plt.show()

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111, projection='3d')
ax = scplt.plot_pca(ax, pdata_fig2_region_norm_impute, classes=['region'], force=True, plot_pc = [1,2,3], cmap=cns_color)
zlab = ax.get_zlabel()
ax.set_zlabel("")
ax.text2D(1.05, 0.55, zlab,
          transform=ax.transAxes, rotation=90,
          ha='left', va='center', fontsize=9.5)


# plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_PCA-3D-impute.svg", bbox_inches="tight", pad_inches=0.1)
plt.show()


# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
umap_params = {'n_neighbors': 20, 'min_dist': 0.1}
ax, umap = scplt.plot_umap(ax, pdata_fig2_region_norm_impute, classes = 'region', s=20, alpha=.8, force = True, umap_params=umap_params, return_fit=True, cmap=cns_color)
scplt.shift_legend(ax)

legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_UMAP-impute.svg", bbox_inches="tight", pad_inches=0.1)
# plt.show()


# ### de

# #### directlfq

# In[ ]:


case_values = [{'region': 'Cortex'}, {'region': 'SNpc'}]

color_dict={'upregulated':"#D19DCB", 'downregulated': "#85BE9E",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(3, 3))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig2_region_norm, values=case_values, pval=0.05, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[0,0],
                                    group_annot_kwargs={"pos": {"group1_xy": (0.98, 1.09), "group2_xy": (0.02, 1.09)}},)

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_VOLCANO-directlfq.svg", bbox_inches="tight")
# volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_DE-directlfq.csv")


# In[ ]:


case_values = [{'region': 'Cortex'}, {'region': 'SNpc'}]

color_dict={'upregulated':"#D19DCB", 'downregulated': "#85BE9E",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig2_region_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', label=[0,0], color=color_dict)

ax, texts = scplt.mark_volcano(ax, volcano_df, label=["Sncg", "Ddc", "Th","Slc17a7", "Camk2a","Gria3", "Homer1", "Aldh1a1", "Calb2","Akr1b7"], label_color='k', return_texts=True)

# scplt.mark_volcano_by_significance(ax, volcano_df, label=["Camk2a", "Gria3", "Icam5", "Homer1", "Slc17a7", "Th", "Ddc", "Sncg", "Slc6a3", "Calb2", "Aldh1a1", "Akr1b7"], color=color_dict)
scplt.volcano_adjust_and_outline_texts(texts, expand=(3,3))

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_VOLCANO-directlfq-highlighted.svg", bbox_inches="tight", pad_inches=0.1)


# #### impute

# In[ ]:


case_values = [{'region': 'Cortex'}, {'region': 'SNpc'}]

color_dict={'upregulated':"#D19DCB", 'downregulated': "#85BE9E",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(3, 3))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig2_region_norm_impute,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[0,0],
                                    group_annot_kwargs={"pos": {"group1_xy": (0.98, 1.09), "group2_xy": (0.02, 1.09)}},)

# scplt.mark_volcano(ax, volcano_df, label=['Snca','Vamp2'], label_color='k')

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_VOLCANO-impute.svg", bbox_inches="tight", pad_inches=0.1)
# volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig2-CNS_DE-impute.csv")


# In[ ]:


pdata_fig2_region_norm.export_layer("X", "proteomics_formatted.csv", obs_names='region', var_names='Genes')


# ## PNS

# ### pre-processing

# In[ ]:


obs_columns = ['sample', 'region', 'size']
# obs_columns = ['sample', 'none', 'load', 'none', 'none', 'none', 'none', 'PBS', 'none']

pdata_pns = pAnnData.import_proteomeDiscoverer(prot_file = 'data/2401_all/Marion_20240116_OTE_Aur60min_FFmouse_Proteins.txt', pep_file = 'data/2401_all/Marion_20240116_OTE_Aur60min_FFmouse_PeptideGroups.txt', obs_columns = obs_columns)
pdata_pns


# In[ ]:


values_pns = [
    {'region': 'mp_axon'},
    {'region': 'mp_cellbody'},]

pdata_pns_filter = pdata_pns.filter_sample(values=values_pns, exact_cases=True)
pdata_pns_filter.summary["region"] = pdata_pns_filter.summary["region"].replace({
    "mp_axon": "NB",
    "mp_cellbody": "Ganglia"
})
pdata_pns_filter.update_summary()

pdata_pns_processed = pdata_pns_filter.copy()
pdata_pns_processed = pdata_pns_processed.filter_prot_found(min_ratio=0.4, group = ['region'], match_any=True)
pdata_pns_processed = pdata_pns_processed.filter_prot(valid_genes=True, unique_profiles=True)

pdata_pns_norm = pdata_pns_processed.copy()
pdata_pns_norm.normalize(method='directlfq')

pd.DataFrame(pdata_pns_norm.summary)


# In[ ]:


pdata_pns_norm.summary.region.value_counts()


# ### plots

# In[ ]:


pns_color={
    "Ganglia": "#F8AF88",
    "NB": "#A8A8FF"
}


# In[ ]:


order = ['NB', 'Ganglia']

fig, ax = plt.subplots(figsize=(2, 3))
sns.barplot(x='region', y='protein_count',
            data=pdata_pns_filter.summary, errorbar='sd', capsize=.1,
            saturation=1, alpha=0.5, palette=pns_color,width=0.75,
            order=order, ax=ax, edgecolor='black',      # outline color
            linewidth=0.8)
sns.swarmplot(x='region', y='protein_count',
                data=pdata_pns_filter.summary, color='k', ax=ax)

from scipy.stats import ttest_ind
group_NB = pdata_pns_filter.summary[pdata_pns_filter.summary['region'] == 'NB']['protein_count']
group_ganglia = pdata_pns_filter.summary[pdata_pns_filter.summary['region'] == 'Ganglia']['protein_count']
max_prot_count = pdata_pns_filter.summary['protein_count'].max()

scplt.plot_significance(ax, max_prot_count+350, 100,
                        pval=ttest_ind(group_NB, group_ganglia).pvalue,
                        fontsize=11)

plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('Region', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 2500)
# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_FigS2_protein-count.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


fig, ax = plt.subplots(figsize=(3, 3))
ax, venn_contents = scplt.plot_venn(ax, pdata_pns_filter,
                                    classes='region',
                                    set_colors=[pns_color['NB'], pns_color['Ganglia']],
                                    return_contents=True)
plt.rcParams.update({'font.size': 12})
# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_FigS2_venn.svg", bbox_inches="tight", pad_inches=0.1)

setN = set(venn_contents["NB"])
setY = set(venn_contents["Ganglia"])
venn_df = pd.DataFrame({
    "NB_only": pd.Series(sorted(setN - setY)),
    "Ganglia_only": pd.Series(sorted(setY - setN)),
    "Both":      pd.Series(sorted(setN & setY))
})


# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
ax, pca = scplt.plot_pca(ax, pdata_pns_norm, classes='region', s=20, alpha=.8, return_fit=True, force=True, cmap=pns_color, add_ellipses=False)
scplt.shift_legend(ax)

# legend=ax.get_legend()
# for patch in legend.legend_handles:
#     patch.set_edgecolor('black')

plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS2_PCA.svg", bbox_inches="tight", pad_inches=0.1)
# plt.show()


# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
umap_params={'min_dist': 0.1, 'n_neighbors': 6}
ax, umap = scplt.plot_umap(ax, pdata_pns_norm, classes = 'region', s=20, alpha=.8, force = True, umap_params=umap_params, return_fit=True, cmap=pns_color,)
scplt.shift_legend(ax)

# legend=ax.get_legend()
# for patch in legend.legend_handles:
#     patch.set_edgecolor('black')

plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS2_UMAP.svg", bbox_inches="tight", pad_inches=0.1)


# ### de

# In[ ]:


case_values = [{'region': 'Ganglia'}, {'region': 'NB'}]

color_dict={'upregulated':"#FC9D69", 'downregulated': "#8D8DFD",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_pns_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[0,0])

scplt.mark_volcano(ax, volcano_df, label=['Snca','Vamp2'], label_color='k')

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS2_VOLCANO.svg", bbox_inches="tight", pad_inches=0.1)
volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_FigS2_DE.csv")


# In[ ]:


text_kwargs=dict(
    fontsize=11,
    color='black',
    offset=1,         # vertical offset above anchor
)

order = ['Ganglia', 'NB']

# TODO: add option for order, ganglia first then NB
pdata_pns_filter.plot_abundance_boxgrid(namelist=['Vamp2'], classes="region", fig_width=1.7, fig_height=3.18, label_x=True, global_legend=False, box=True, show_n=False, text_kwargs=text_kwargs, y_max=11, palette=pns_color, order=order)
pdata_pns_filter.plot_abundance_boxgrid(namelist=['Snca'], classes="region", fig_width=1.7, fig_height=3.18, label_x=True, global_legend=False, box=True, show_n=False, text_kwargs=text_kwargs, y_max=11, palette=pns_color, order=order)

# pdata_pns_norm.plot_abundance_boxgrid(namelist=['Vamp2'], classes="region", fig_width=1.7, fig_height=3.18, label_x=True, global_legend=False, box=True, show_n=False, text_kwargs=text_kwargs, y_max=11, palette=pns_color, order=order)
# pdata_pns_norm.plot_abundance_boxgrid(namelist=['Snca'], classes="region", fig_width=1.7, fig_height=3.18, label_x=True, global_legend=False, box=True, show_n=False, text_kwargs=text_kwargs, y_max=11, palette=pns_color, order=order)


# # Fig 3 - Omics Analysis

# See `analysis_2512_omics_clean.ipynb`

# # Fig 5 - snpc subtypes

# ## anxa

# ### pre-processing

# In[3]:


pdata_fig3_snpc_filter = pdata_all_region_filtered.filter_sample(condition ='Grouping == "snpc"')
pdata_fig3_snpc_filter = pdata_fig3_snpc_filter.filter_sample(condition = 'Batch == "Apr-25"')
pdata_fig3_snpc_filter = pdata_fig3_snpc_filter.filter_prot_significant()
pdata_fig3_snpc_filter = pdata_fig3_snpc_filter.filter_sample(min_prot=1000)

pdata_fig3_snpc_filter.summary["cell_type"] = pdata_fig3_snpc_filter.summary["Sub grouping"].replace({
    "anxaY": "Anxa+",
    "anxaN": "Anxa-"
})
pdata_fig3_snpc_filter.update_summary()

pdata_fig3_snpc_processed = pdata_fig3_snpc_filter.copy()
pdata_fig3_snpc_processed = pdata_fig3_snpc_processed.filter_prot_found(min_ratio=0.4, group = ['Sub grouping'], match_any=True)
pdata_fig3_snpc_processed = pdata_fig3_snpc_processed.filter_prot(valid_genes=True, unique_profiles=True)

pdata_fig3_snpc_norm = pdata_fig3_snpc_processed.copy()
pdata_fig3_snpc_norm.normalize(method='directlfq')


# In[5]:


pdata_fig3_snpc_norm.summary.cell_type.value_counts()


# In[ ]:


pd.DataFrame(pdata_fig3_snpc_norm.summary)


# ### plots

# In[4]:


anxa_color={
    "Anxa-": "#A0FFA0",
    "Anxa+": "white"
}


# In[7]:


order = ['Anxa+', 'Anxa-']

dataset = pdata_fig3_snpc_filter

fig, ax = plt.subplots(figsize=(2, 3))
sns.barplot(x='cell_type', y='protein_count',
            data=dataset.summary, errorbar='sd', capsize=.1,
            saturation=1, alpha=0.5, palette=anxa_color,width=0.75,
            order=order, ax=ax,
            linewidth=0.8, edgecolor='k')
sns.swarmplot(x='cell_type', y='protein_count',
                data=dataset.summary, color='k', ax=ax)

from scipy.stats import ttest_ind
group_anxaY = dataset.summary[dataset.summary['cell_type'] == 'Anxa+']['protein_count']
group_anxaN = dataset.summary[dataset.summary['cell_type'] == 'Anxa-']['protein_count']
max_prot_count = dataset.summary['protein_count'].max()

scplt.plot_significance(ax, max_prot_count+350, 100,
                        pval=ttest_ind(group_anxaY, group_anxaN).pvalue,
                        fontsize=11)

plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('Region', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 4200)
# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_protein-count.svg", bbox_inches="tight", pad_inches=0.1)


# In[8]:


order = ['Anxa+', 'Anxa-']

dataset = pdata_fig3_snpc_filter

fig, ax = plt.subplots(figsize=(2, 3))
sns.violinplot(x='cell_type', y='protein_count',
            data=dataset.summary,saturation=1, alpha=0.5, palette=anxa_color,width=0.75,
            order=order, ax=ax, edgecolor='black', inner="points",
            linewidth=0.8)
# sns.swarmplot(x='cell_type', y='protein_count',
#                 data=dataset.summary, color='k', ax=ax)

from scipy.stats import ttest_ind
group_anxaY = dataset.summary[dataset.summary['cell_type'] == 'Anxa+']['protein_count']
group_anxaN = dataset.summary[dataset.summary['cell_type'] == 'Anxa-']['protein_count']
max_prot_count = dataset.summary['protein_count'].max()

scplt.plot_significance(ax, max_prot_count+550, 100,
                        pval=ttest_ind(group_anxaY, group_anxaN).pvalue,
                        fontsize=11)

plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('Region', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 4200)
plt.show()
# plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_protein-count-violin.svg", bbox_inches="tight", pad_inches=0.1)


# In[9]:


# CV

fig, ax = plt.subplots(figsize=(3,3))
cv_df = scplt.plot_cv(ax, dataset, classes = 'cell_type', return_df=True)

cv_df = cv_df.reset_index()
sns.violinplot(data=cv_df, x='Class', y='CV', palette=anxa_color, linewidth=1, inner='quartile', saturation = 1, ax=ax)
# ax.set_title(classes)
# legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, order)]
# plt.legend(handles=legend_patches, bbox_to_anchor=(.75, 1), loc=2, borderaxespad=0., frameon=False)

from matplotlib.ticker import PercentFormatter
ax.yaxis.set_major_formatter(PercentFormatter(1))
# ax.set_xlabel('Area ($\mu m^2$)')
ax.set_ylabel('Protein Abundance CV')
ax.figure.set_size_inches(2.5, 3)

# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_CV.svg", bbox_inches="tight", pad_inches=0.1)


# In[11]:


fig, ax = plt.subplots(figsize=(3, 3))
dataset = pdata_fig3_snpc_filter
plt.rcParams.update({'font.size': 15})
ax, venn_contents = scplt.plot_venn(ax, dataset,
                                    classes='cell_type',
                                    set_colors=["#A0FFA0", '#FFFFFF'],
                                    return_contents=True,
                                    label_order=['Anxa-','Anxa+'],
                                    weighted=True, fixed_subset_sizes=(1, 1, 2.5),)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Anxa+", "Anxa-"}]:
    x, y = t.get_position()
    if t.get_text() == "Anxa-":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")

# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_venn.svg", bbox_inches="tight", pad_inches=0.1)

setAnxaY = set(venn_contents["Anxa+"])
setAnxaN = set(venn_contents["Anxa-"])
venn_df = pd.DataFrame({
    "anxaY_only": pd.Series(sorted(setAnxaY - setAnxaN)),
    "anxaN_only": pd.Series(sorted(setAnxaN - setAnxaY)),
    "Both":      pd.Series(sorted(setAnxaY & setAnxaN))
})


# In[12]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
ax, pca = scplt.plot_pca(ax, pdata_fig3_snpc_norm, classes='cell_type', s=20, alpha=.8, return_fit=True, force=True, cmap=anxa_color, add_ellipses=False)
scplt.shift_legend(ax)

for coll in ax.collections:
    coll.set_edgecolor("black")
    coll.set_linewidth(0.4)

legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_PCA.svg", bbox_inches="tight", pad_inches=0.1)
# plt.show()

# make 3d figure 
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')
ax = scplt.plot_pca(ax, pdata_fig3_snpc_norm, classes=['cell_type'], force=True, plot_pc = [1,2,3], cmap=anxa_color, s=30)
zlab = ax.get_zlabel()
ax.set_zlabel("")
ax.text2D(1.1, 0.55, zlab,
          transform=ax.transAxes, rotation=90,
          ha='left', va='center', fontsize=9.5)

for coll in ax.collections:
    coll.set_edgecolor("black")
    coll.set_linewidth(0.4)

legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')
if legend is not None:
    legend.remove()

plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_PCA-3D.svg", bbox_inches="tight", pad_inches=0.1)
# plt.show()


# In[14]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
umap_params = {'n_neighbors': 5, 'min_dist': 0.1}
ax, umap = scplt.plot_umap(ax, pdata_fig3_snpc_norm, classes = 'cell_type', s=20, alpha=.8, force = True, umap_params=umap_params, return_fit=True, cmap=anxa_color)
scplt.shift_legend(ax)

legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')

for coll in ax.collections:
    coll.set_edgecolor("black")
    coll.set_linewidth(0.4)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_UMAP.svg", bbox_inches="tight", pad_inches=0.1)
# plt.show()


# ### de

# In[15]:


case_values = [{'cell_type': 'Anxa+'}, {'cell_type': 'Anxa-'}]

color_dict={'upregulated':anxa_color['Anxa+'], 'downregulated': anxa_color['Anxa-'],'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig3_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[5,5])

white = np.array([1.0, 1.0, 1.0])
green = np.array([160/255, 1.0, 160/255])  # #A0FFA0

for coll in ax.collections:
    fc = coll.get_facecolors()
    if fc is None or fc.size == 0:
        continue

    mask_white = (fc[:, :3] == white).all(axis=1)
    mask_green = (fc[:, :3] == green).all(axis=1)
    mask = mask_white | mask_green

    if mask.any():
        ec = np.zeros_like(fc)          # transparent edges by default
        ec[mask] = [0, 0, 0, 1]         # black outline
        coll.set_edgecolors(ec)

        lw = np.zeros(len(fc))
        lw[mask] = 0.4
        coll.set_linewidths(lw)

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_VOLCANO-directlfq-top5.svg", bbox_inches="tight", pad_inches=0.1)
# volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_DE-directlfq.csv")


# In[ ]:


# HIGHLIGHTED GENES

case_values = [{'cell_type': 'Anxa+'}, {'cell_type': 'Anxa-'}]

color_dict={'upregulated':anxa_color['Anxa+'], 'downregulated': anxa_color['Anxa-'],'not_significant': "#000000"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig3_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', no_marks=True)

scplt.mark_volcano_by_significance(ax, volcano_df, label=['Aldh1a1','Calb1','Anxa1','Th','Gfap'], color=color_dict)

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_VOLCANO-directlfq-highlighted.svg", bbox_inches="tight", pad_inches=0.1)


# In[7]:


case_values = [{'cell_type': 'Anxa+'}, {'cell_type': 'Anxa-'}]

color_dict={'upregulated':anxa_color['Anxa+'], 'downregulated': anxa_color['Anxa-'],'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig3_snpc_norm,
                                    values=case_values,
                                    pval=0.1, log2fc=0.5, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[0,0],
                                    up_kwargs={"color":"k"}, down_kwargs={"color":"k"}, group_annot=False)

white = np.array([1.0, 1.0, 1.0])
green = np.array([160/255, 1.0, 160/255])  # #A0FFA0

for coll in ax.collections:
    fc = coll.get_facecolors()
    if fc is None or fc.size == 0:
        continue

    mask_white = (fc[:, :3] == white).all(axis=1)
    mask_green = (fc[:, :3] == green).all(axis=1)
    mask = mask_white | mask_green

    if mask.any():
        ec = np.zeros_like(fc)          # transparent edges by default
        ec[mask] = [0, 0, 0, 1]         # black outline
        coll.set_edgecolors(ec)

        lw = np.zeros(len(fc))
        lw[mask] = 0.4
        coll.set_linewidths(lw)

sig_color_dict={'upregulated':'k', 'downregulated': 'k','not_significant': "#FFFFFF6A"}

texts = []

ax, t = scplt.mark_volcano_by_significance(ax, volcano_df, label=['Tmem14c','Aldh3b1','Baiap3','Strn4','Amy1','Fuca1','Plbd2','Rnf13'], color=color_dict, text_color='k', return_texts=True)
texts.extend(t)
ax, t = scplt.mark_volcano_by_significance(ax, volcano_df, label=['Calb1','Anxa1','Th','Aldh1a1'], color=sig_color_dict, return_texts=True, s=15)
texts.extend(t)

scplt.volcano_adjust_and_outline_texts(texts, expand=(2.5, 2.5))

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_VOLCANO-directlfq-top5.svg", bbox_inches="tight", pad_inches=0.1)
# volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_DE-directlfq-relaxed.csv")


# In[ ]:


genes = [
    "Calm1",
    "Th",
    "Map7d2",
    "Prrc2c",
    "Aldh1a1",
    "Rtn3",
    "Anks1b",
    "Srrm2",
    "Pak3",
    "Gaa",
    "Stmn3",
    "Timp2",
    "Ddb1",
    "Atp9a",
    "Tpi1",
    "Map2",
    "Cdk14",
    "Tsn",
    "Calb2",
    "Strn3",
    "Rit2",
    "Nwd2",
    "Apbb1",
    "Pitpna",
    "Gria3",
    "Strbp",
    "Naa15",
    "Scn2b",
    "Chmp4b",
    "Cadps2",
    "Rnf11",
    "Calb1",
    "Ids",
    "Tecpr1",
    "Pde2a",
    "Agap3",
    "Tnpo3",
    "Rasgrf2",
    "Hap1",
    "Madd",
    "Nop58",
    "Arhgap44",
    "Chl1",
    "Bcat1",
    "Ccdc136",
    "Cygb",
    "Upf2",
    "Cd99l2",
    "Usp7",
    "Dcaf7",
    "Rab1b",
    "Plxna1",
    "Adcy5",
    "Nrgn",
    "Hectd1",
    "Cmas",
    "Exoc5",
    "Ago2",
    "Cdc42bpb",
    "Hsdl1",
    "Wdfy3",
    "Rabgap1",
    "Setd7",
    "Tango2",
    "Phyhip",
    "Fbxo41",
    "Ckap4",
    "Tln2",
    "Lypla2",
    "Clcn6",
    "Ocrl",
    "Fxyd7",
    "Fxyd6",
    "Wipf2",
    "Pde1b",
    "Fbxo3",
    "Alg2",
    "Ralyl",
    "Nrbp1",
    "Pigk",
    "Wdr37",
    "Ipo9",
    "Csnk1d",
    "Atat1",
    "Xpo1",
    "Diras2",
    "Anxa7",
    "Psma3",
    "Asns",
    "Prmt5",
    "Nt5dc3",
    "Pnpo",
    "Uba2",
    "Abhd16a",
    "Sae1",
    "Nagk",
    "Igsf21",
    "Ppa1",
    "Rab3gap2",
    "Strn4",
    "Ssx2ip",
    "Trappc9",
    "Spast",
    "Ipo5",
    "Itpa",
    "Pde10a",
    "Pef1",
    "Spcs3",
    "Cap2",
    "Syt5",
    "Agk",
    "Camk2a",
    "Cstf2",
    "Arhgap26",
    "Fam169a",
    "Fbll1",
    "Srprb",
    "Vps18",
    "Vamp1",
    "Baiap3",
    "Ttc7b",
    "Rab9b",
    "Fsd1",
    "Scyl2",
    "Nutf2",
    "Ptpn9",
    "Impdh1",
    "Tmem263",
    "Slc9a1",
    "Dclk2",
    "Sdr39u1",
    "Camsap3",
    "Ccar2",
    "Rragb",
    "Sf3b3",
    "Prrc1",
    "Anxa1",
    "Lrrtm1",
    "Ccdc177",
    "Igf2r",
    "Capn1",
    "Focad",
    "Gba1",
]

annotate_genes=[
'Calb1',
'Calb2',
# 'Pak3',
# 'Scn2b',
'Aldh1a1',
'Anxa1',
]

genes2 = [
    "Itm2c",
    "Reep5",
    "Mobp",
    "Aplp1",
    "Ldhb",
    "Ralbp1",
    "Ubb",
    "Ptprd",
    "Eif5",
    "Metap2",
    "Ctsd",
    "Sparcl1",
    "Csnk1a1",
    "Eef2",
    "Rrbp1",
    "Agpat4",
    "Fkbp2",
    "Ociad1",
    "Prdx2",
    "Map7d1",
    "Capns1",
    "Lamp1",
    "Anp32a",
    "Psma7",
    "Pcp4",
    "Arpc2",
    "Hint1",
    "Gapvd1",
    "Psma4",
    "Park7",
    "Endod1",
    "Apod",
    "Ssr1",
    "Rps2",
    "Nudt3",
    "Rpn1",
    "Sec61a1",
    "Ubl3",
    "Ablim1",
    "Ube2k",
    "Npm1",
    "Tpm3",
    "Psmb4",
    "Actn4",
    "Phactr1",
    "Rac1",
    "Psmb3",
    "Tbcb",
    "Psma6",
    "Eif3f",
    "Ctnna1",
    "Ndufa11",
    "Pafah1b2",
    "Psmb1",
    "Tmed10",
    "Acox1",
    "Cnp",
    "Cldn11",
    "Trim9",
    "Gns",
    "Pa2g4",
    "Psmb6",
    "Pgp",
    "Pithd1",
    "Snx30",
    "Slc1a4",
    "Pitpnc1",
    "Tspan2",
    "Csnk2a1",
    "Dad1",
    "Prdx1",
    "Omg",
    "Psme1",
    "Mydgf",
    "Rpn2",
    "Ap1s2",
    "Safb",
    "Ablim2",
    "Bcas1",
    "Nf1",
    "Psma1",
    "Hip1r",
    "Pgm2",
    "Cluh",
    "Aimp1",
    "Bag6",
    "Sbds",
    "Dtd1",
    "Camk1",
    "Sbf1",
    "Blmh",
    "Ddost",
    "Ppp1r7",
    "Eif2a",
    "Rab5c",
    "Rnf13",
    "Fnbp1",
    "Rplp2",
    "Ahcyl2",
    "Arpp19",
    "Eif3h",
    "Sccpdh",
    "Snx17",
    "Xpo7",
    "Psma2",
    "Snrpd3",
    "Ndrg1",
    "Cap1",
    "Copb2",
    "Taldo1",
    "Abcd3",
    "Pgrmc2",
    "Slc25a1",
    "Kat6b",
    "Plec",
    "Adprh",
    "Cryzl1",
    "Flii",
    "Psmc6",
    "Ptma",
    "Acot13",
    "Eef1d",
    "Lta4h",
    "Psma5",
    "Psmb2",
    "Anp32e",
    "Pdk1",
    "Slc44a1",
    "Taf10",
    "Erlin2",
    "Eif3l",
    "Idh1",
    "Me1",
    "S100a10",
    "Arcn1",
    "Mpp2",
    "Rpl29",
    "Snrpf",
    "Upf1",
    "Aida",
    "Dnaja3",
    "Eif4e",
    "Fam50a",
    "Stt3a",
    "Tcea1",
    "Asah1",
    "Cirbp",
    "Emc2",
    "Fabp5",
    "Lrrc8a",
    "Lsm4",
    "Lsm6",
    "Smu1",
    "Vps26a",
    "Mogs",
    "Ermp1",
    "Gsn",
    "Myl6",
    "Rpl31",
    "Gad2",
    "Slc22a23",
    "Marcksl1",
    "Ide",
    "Rpl38",
    "Rab22a",
    "Vps11",
    "Asap2",
    "Eif4a3",
    "Gnpda2",
    "Lmna",
    "Smarce1",
    "Snrnp200",
    "Coro7",
    "Exoc7",
    "Lamp2",
    "Mpp1",
    "Opalin",
    "Sord",
    "Tmem14c",
    "Tubg1",
    "Actr3b",
    "Aldh7a1",
    "Ehd1",
    "Psmb5",
    "Psmb7",
    "Sorcs2",
    "Isyna1",
    "Man2c1",
    "Mrpl49",
    "Scpep1",
    "Baiap2",
    "Ago1",
    "Nt5c3b",
    "Chtop",
    "Plbd2",
    "Gtpbp1",
    "Magoh",
    "Plcl1",
    "Fuom",
    "Adcy9",
    "Aqp4",
    "Eml4",
    "Flnb",
    "Fuca1",
    "Ggct",
    "Hsd17b4",
    "Kcnj10",
    "Rab23",
    "Rap2b",
    "Siae",
    "Snrpe",
    "Ebp",
    "Eml1",
    "Fam120c",
    "Mrps21",
    "Nt5c",
    "Pdcd5",
    "Smad2",
    "Ubl7",
    "Cat",
    "Fam162a",
    "Gart",
    "Hnrnph2",
    "Idh2",
    "Imp4",
    "Jup",
    "Lnpep",
    "Mtap",
    "Myh9",
    "Pepd",
    "Rpa1",
    "Rpia",
    "Slc1a3",
    "Thop1",
    "Ak4",
    "Ddah2",
    "Flna",
    "Hspa2",
    "Lsm2",
    "Msra",
]


# In[ ]:


case_values = [{'cell_type': 'Anxa+'}, {'cell_type': 'Anxa-'}]

color_dict={'upregulated':anxa_color['Anxa+'], 'downregulated': anxa_color['Anxa-'],'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig3_snpc_norm,
                                    values=case_values,
                                    pval=0.1, log2fc=0.5, return_df=True,
                                    fold_change_mode='mean', no_marks=True)

# color_dict={'upregulated':'k', 'downregulated': 'k','not_significant': "#FFFFFF6A"}


scplt.mark_volcano_by_significance(ax, volcano_df, label=genes, color=color_dict, show_names=False)
scplt.mark_volcano(ax, volcano_df, label=annotate_genes, label_color='k', show_names=True)

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
# plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_VOLCANO-directlfq-top5.svg", bbox_inches="tight", pad_inches=0.1)
# volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_DE-directlfq-relaxed.csv")


# In[ ]:


case_values = [{'cell_type': 'Anxa+'}, {'cell_type': 'Anxa-'}]

color_dict={'upregulated':anxa_color['Anxa+'], 'downregulated': anxa_color['Anxa-'],'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig3_snpc_norm,
                                    values=case_values,
                                    pval=0.1, log2fc=0.5, return_df=True,
                                    fold_change_mode='mean', no_marks=True)

# color_dict={'upregulated':'k', 'downregulated': 'k','not_significant': "#FFFFFF6A"}


scplt.mark_volcano(ax, volcano_df, label=genes, label_color="red", show_names=False)
scplt.mark_volcano(ax, volcano_df, label=genes2, label_color="k", show_names=False)

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
# plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_VOLCANO-directlfq-top5.svg", bbox_inches="tight", pad_inches=0.1)
# volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_DE-directlfq-relaxed.csv")


# In[ ]:


volcano_df[volcano_df['Genes'].isin(genes)].significance.value_counts() # GREEN


# In[ ]:


volcano_df[volcano_df['Genes'].isin(genes2)].significance.value_counts() # YELLOW


# In[ ]:


text_kwargs=dict(
    fontsize=9,
    color='black',
    offset=0.35,         # vertical offset above anchor
)

fig, ax = pdata_fig3_snpc_norm.plot_abundance_boxgrid(namelist=['Calb1'], classes='cell_type', box=True, show_n=True, fig_width=1.7, fig_height=3.18,label_x=True, global_legend=False, palette=anxa_color, text_kwargs=text_kwargs)
fig, ax = pdata_fig3_snpc_norm.plot_abundance_boxgrid(namelist=['Anxa1'], classes='cell_type', box=True, show_n=True, fig_width=1.7, fig_height=3.18,label_x=True, global_legend=False, palette=anxa_color, text_kwargs=text_kwargs)
fig, ax = pdata_fig3_snpc_norm.plot_abundance_boxgrid(namelist=['Aldh1a1'], classes='cell_type', box=True, show_n=True, fig_width=1.7, fig_height=3.18,label_x=True, global_legend=False, palette=anxa_color, text_kwargs=text_kwargs)


# In[ ]:


# for transcriptomics
# try p=0.1, fc=0.5

# filter out Krt, filter out the "neuronal-like" proteins/genes compared to astrocyte
# need some ML algorithm to output a single confidence/quantity/marker


# In[ ]:


pdata_fig3_snpc_processed.export_layer("X", f"MANUSCRIPT/MANUSCRIPT_Fig3-Anxa_processed.csv", transpose=True, var_names='Genes')


# In[ ]:


pdata_fig3_snpc_norm.export_layer("X", "proteomics_formatted_anxa.csv", obs_names='cell_type', var_names='Genes')


# # Fig 4 - non-neuronal cells

# ## astrocyte vs neuron

# ### pre-processing

# In[49]:


# take subset of data from pdata_all_region_filtered

pdata_fig4_region_filtered = pdata_all_region_filtered.filter_prot_significant()
pdata_fig4_region_filtered = pdata_fig4_region_filtered.filter_sample(condition="Grouping != 'snpc'")
pdata_fig4_region_filtered = pdata_fig4_region_filtered.filter_sample(condition="Batch != 'Apr-25'")

mapping = {
    "ctx": "Neuron",
    "ast": "Astrocyte",
}
pdata_fig4_region_filtered.summary["cell_type"] = pdata_fig4_region_filtered.summary["Grouping"].replace(mapping)
pdata_fig4_region_filtered.update_summary()


# In[53]:


pdata_fig4_region_filtered.prot.obs.cell_type.value_counts()


# In[5]:


pdata_fig4_region_processed = pdata_fig4_region_filtered.copy()
pdata_fig4_region_processed = pdata_fig4_region_processed.filter_sample(min_prot=1800)

pdata_fig4_region_processed = pdata_fig4_region_processed.filter_prot_found(min_ratio=0.4, group = ['cell_type'], match_any=True)
pdata_fig4_region_processed = pdata_fig4_region_processed.filter_prot(valid_genes=True, unique_profiles=True)

pdata_fig4_region_norm = pdata_fig4_region_processed.copy()
pdata_fig4_region_norm.normalize(method='directlfq')

print(pdata_fig4_region_norm.summary['Grouping'].value_counts())
print(pdata_fig4_region_norm.summary['Batch'].value_counts())
print(pdata_fig4_region_norm.summary['cell_type'].value_counts())


# In[ ]:


pdata_fig4_region_processed.prot.obs.cell_type.value_counts()


# In[ ]:


pd.DataFrame(pdata_fig4_region_norm.summary)


# In[ ]:


pdata_fig4_region_processed.export_layer("X", f"MANUSCRIPT/MANUSCRIPT_Fig4-region_processed.csv", transpose=True, var_names='Genes')


# ### plots

# In[51]:


cell_type_color={
    "Neuron": "#D19DCB",
    "Astrocyte": "#84C8D5" # "#CAEEFB"
}


# In[ ]:


order = ['Neuron', 'Astrocyte']

dataset = pdata_fig4_region_filtered

fig, ax = plt.subplots(figsize=(2, 3))
sns.barplot(x='cell_type', y='protein_count',
            data=dataset.summary, errorbar='sd', capsize=.1,
            saturation=1, alpha=0.5, palette=cell_type_color,width=0.75,
            order=order, ax=ax, edgecolor='black',      # outline color
            linewidth=0.8)
sns.swarmplot(x='cell_type', y='protein_count',
                data=dataset.summary,order=order, color='k', ax=ax)

from scipy.stats import ttest_ind
neuronStats = dataset.summary[dataset.summary['cell_type'] == 'Neuron']['protein_count']
astrocyteStats = dataset.summary[dataset.summary['cell_type'] == 'Astrocyte']['protein_count']

max_prot_count = dataset.summary['protein_count'].max()

scplt.plot_significance(ax, max_prot_count+350, 100,
                        pval=ttest_ind(neuronStats, astrocyteStats).pvalue,
                        fontsize=11)

plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('Type', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 4200)
# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-region_protein-count.svg", bbox_inches="tight", pad_inches=0.1)


# In[54]:


dataset = pdata_fig4_region_filtered

fig, ax = plt.subplots(figsize=(3, 3))
plt.rcParams.update({'font.size': 12})
ax, venn_contents = scplt.plot_venn(ax, dataset,
                                    classes='cell_type',
                                    set_colors=["#D19DCB", "#84C8D5"],
                                    return_contents=True, weighted=True, fixed_subset_sizes=(1, 1, 3))

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Neuron", "Astrocyte"}]:
    x, y = t.get_position()
    if t.get_text() == "Neuron":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")

# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-region_venn.svg", bbox_inches="tight", pad_inches=0.1)

setNeuron = set(venn_contents["Neuron"])
setAstrocyte = set(venn_contents["Astrocyte"])
venn_df = pd.DataFrame({
    "Neuron_only": pd.Series(sorted(setNeuron - setAstrocyte)),
    "Astrocyte_only": pd.Series(sorted(setAstrocyte - setNeuron)),
    "Both":      pd.Series(sorted(setNeuron & setAstrocyte))
})


# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
ax, pca = scplt.plot_pca(ax, pdata_fig4_region_norm, classes=['cell_type'], s=20, alpha=.8,pca_params={'n_comps': 9}, force=True, cmap=cell_type_color, add_ellipses=False, return_fit=True)
scplt.shift_legend(ax)
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-region_PCA.svg", bbox_inches="tight", pad_inches=0.1)

fig, ax = plt.subplots(1,1, figsize = (3,3))
umap_params={'min_dist': 0.3, 'n_neighbors': 7}
ax = scplt.plot_umap(ax, pdata_fig4_region_norm, classes=['cell_type'], s=20, alpha=.8, force=True, umap_params=umap_params, cmap=cell_type_color)
scplt.shift_legend(ax)
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-region_UMAP.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


scutils.get_pca_importance(pca, pdata_fig4_region_norm.prot.var_names, 20)


# ### de

# In[7]:


case1 = {'cell_type': 'Neuron'}
case2 = {'cell_type': 'Astrocyte'}
case_values = [case1, case2]

cell_type_color={
    "Neuron": "#D19DCB",
    "Astrocyte": "#84C8D5" # "#CAEEFB"
}
color_dict={'upregulated': cell_type_color['Neuron'], 'downregulated': cell_type_color['Astrocyte'],'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(3.63, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig4_region_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[10,10])

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig4-region_VOLCANO-directlfq-top10.svg", bbox_inches="tight", pad_inches=0.1)
# volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig4-region_DE-directlfq.csv")


# In[ ]:


case1 = {'cell_type': 'Neuron'}
case2 = {'cell_type': 'Astrocyte'}
case_values = [case1, case2]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig4_region_norm, values = case_values, pval = 0.05, return_df=True, no_marks=True)
scplt.mark_volcano_by_significance(ax, volcano_df, label=['Syn1','Dlg4','Sncg','Sncb','Snca','Fabp7','Gfap'], color=color_dict)
# scplt.mark_volcano(ax, volcano_df, label=['Sncg','Sncb','Snca'],label_color='purple')
# scplt.mark_volcano(ax, volcano_df, label=['Fabp7','Gfap'],label_color='blue')

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig4-region_VOLCANO-directlfq-Highlighted.svg", bbox_inches="tight", pad_inches=0.1)


# ## SWI model

# ### pre-processing

# In[6]:


obs_columns = ['date', 'acquisition', 'size', 'cell_type', 'replicate']

pdata_fig4_swi = pAnnData.import_data(source_type='diann', report_file = 'data/2511_ast-SWI-only/report.parquet', obs_columns=obs_columns)


# In[7]:


pdata_fig4_swi_filtered = pdata_fig4_swi.filter_sample(condition="cell_type not in ['T2','TEAB35']")
pdata_fig4_swi_filtered = pdata_fig4_swi_filtered.filter_prot_significant()
pdata_fig4_swi_filtered = pdata_fig4_swi_filtered.filter_sample(min_prot=1000)

mapping = {
    "scartissue": "SWI",
    "astrocyte": "Uninjured",
}
pdata_fig4_swi_filtered.summary["condition"] = pdata_fig4_swi_filtered.summary["cell_type"].replace(mapping)
pdata_fig4_swi_filtered.update_summary()

pdata_fig4_swi_filtered.prot.obs.condition.value_counts()


# In[8]:


pdata_fig4_swi_processed = pdata_fig4_swi_filtered.copy()
pdata_fig4_swi_processed = pdata_fig4_swi_processed.filter_prot_found(min_ratio=0.4, group = ['condition'], match_any=True)
pdata_fig4_swi_processed = pdata_fig4_swi_processed.filter_prot(valid_genes=True, unique_profiles=True)

pdata_fig4_swi_norm = pdata_fig4_swi_processed.copy()
pdata_fig4_swi_norm.normalize(method='directlfq')


# In[ ]:


pdata_fig4_swi_processed.export_layer("X", f"MANUSCRIPT/MANUSCRIPT_Fig4-SWI_processed.csv", transpose=True, var_names='Genes')


# ### plots

# In[9]:


condition_color={
    "SWI": "#1F98B1",
    "Uninjured": "#84C8D5" # "#CAEEFB"
}


# In[ ]:


text_kwargs=dict(fontsize=11,color='black', offset=1,) # vertical offset above anchor
bar_kwargs=dict(width=0.25, edgecolor="black", linewidth=0.6,)
strip_kwargs=dict(darken_factor=0.55, alpha=0.8)

fig, ax, df = pdata_fig4_swi_processed.plot_abundance_boxgrid(namelist=['Cd68'], classes="condition", fig_width=2, fig_height=3.18, 
    label_x=True, global_legend=False, plot_type="bar", show_n=False, text_kwargs=text_kwargs, bar_kwargs=bar_kwargs, strip_kwargs=strip_kwargs,
    palette=condition_color, order=['Uninjured','SWI'],log_scale=False, return_df=True)

# ymin, ymax = ax[0].get_ylim()
ax[0].set_ylim(0, 900000)

fig, ax, df = pdata_fig4_swi_processed.plot_abundance_boxgrid(namelist=['Itgam'], classes="condition", fig_width=2, fig_height=3.18, 
    label_x=True, global_legend=False, plot_type="bar", show_n=False, text_kwargs=text_kwargs, bar_kwargs=bar_kwargs, strip_kwargs=strip_kwargs,
    palette=condition_color, order=['Uninjured','SWI'],log_scale=False, return_df=True)



# In[18]:


order = ['Uninjured', 'SWI']

dataset = pdata_fig4_swi_filtered

fig, ax = plt.subplots(figsize=(1.77, 3))
sns.barplot(x='condition', y='protein_count',
            data=dataset.summary, errorbar='sd', capsize=.1,
            saturation=1, alpha=0.5, palette=condition_color,width=0.75,
            order=order, ax=ax, edgecolor='black',      # outline color
            linewidth=0.8)
sns.swarmplot(x='condition', y='protein_count',
                data=dataset.summary,order=order, color='k', ax=ax)

from scipy.stats import ttest_ind
statsUninjured = dataset.summary[dataset.summary['condition'] == 'Uninjured']['protein_count']
statsSWI = dataset.summary[dataset.summary['condition'] == 'SWI']['protein_count']

max_prot_count = dataset.summary['protein_count'].max()

# scplt.plot_significance(ax, max_prot_count+350, 100,
#                         pval=ttest_ind(statsUninjured, statsSWI).pvalue,
#                         fontsize=11)

plt.ylabel('Protein Count')
plt.xlabel('')
# plt.ylabel('Protein Count',fontsize=13, labelpad=4)
# plt.xlabel('Type', fontsize=13, labelpad=4)
# ax.tick_params(axis='x', labelsize=11)
# ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 4800)
# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-SWI_protein-count.svg", bbox_inches="tight", pad_inches=0.1)


# In[28]:


fig, ax = plt.subplots(figsize=(3, 3))
plt.rcParams.update({'font.size': 12})
ax, venn_contents = scplt.plot_venn(ax, dataset,
                                    classes='condition',
                                    set_colors=["#84C8D5", "#1F98B1"],
                                    return_contents=True,
                                    weighted=True, fixed_subset_sizes=(1, 1, 3),)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Uninjured", "SWI"}]:
    x, y = t.get_position()
    if t.get_text() == "Uninjured":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")

plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-SWI_venn.svg", bbox_inches="tight", pad_inches=0.1)

setUninjured = set(venn_contents["Uninjured"])
setSWI = set(venn_contents["SWI"])
venn_df = pd.DataFrame({
    "Unijured_only": pd.Series(sorted(setUninjured - setSWI)),
    "SWI_only": pd.Series(sorted(setSWI - setUninjured)),
    "Both":      pd.Series(sorted(setSWI & setUninjured))
})


# In[32]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
ax = scplt.plot_pca(ax, pdata_fig4_swi_norm, classes=['condition'], s=20, alpha=.8,pca_params={'n_comps': 9}, force=True, cmap=condition_color, add_ellipses=False)
for coll in ax.collections:
    coll.set_edgecolor("black")
    coll.set_linewidth(0.4)

scplt.shift_legend(ax)
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-SWI_PCA.svg", bbox_inches="tight", pad_inches=0.1)

# make 3d figure 
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')
ax = scplt.plot_pca(ax, pdata_fig4_swi_norm, classes=['condition'], force=True, plot_pc = [1,2,3], cmap=condition_color)
zlab = ax.get_zlabel()
ax.set_zlabel("")
ax.text2D(1.1, 0.55, zlab,
          transform=ax.transAxes, rotation=90,
          ha='left', va='center', fontsize=9.5)
legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')
if legend is not None:
    legend.remove()
for coll in ax.collections:
    coll.set_edgecolor("black")
    coll.set_linewidth(0.4)

plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-SWI_PCA-3D.svg", bbox_inches="tight", pad_inches=0.1)

fig, ax = plt.subplots(1,1, figsize = (3,3))
umap_params={'min_dist': 0.3, 'n_neighbors': 7}
ax = scplt.plot_umap(ax, pdata_fig4_swi_norm, classes=['condition'], s=20, alpha=.8, force=True, umap_params=umap_params, cmap=condition_color)
scplt.shift_legend(ax)
# plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig4-SWI_UMAP.svg", bbox_inches="tight", pad_inches=0.1)


# ### de

# In[ ]:


case_values = [{'condition': 'SWI'}, {'condition': 'Uninjured'}]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig4_swi_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[10,10])

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig4-SWI_VOLCANO-directlfq-top10.svg", bbox_inches="tight", pad_inches=0.1)
volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig4-SWI_DE-directlfq.csv")


# In[37]:


case_values = [{'condition': 'SWI'}, {'condition': 'Uninjured'}]

color_dict={'upregulated':condition_color['SWI'], 'downregulated': condition_color['Uninjured'],'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(3.34, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig4_swi_norm, values = case_values, pval = 0.05, return_df=True, color=color_dict, label=[0,0])
scplt.mark_volcano(ax, volcano_df, label=['Aldh1l1','Gfap'], label_color='k')

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig4-SWI_VOLCANO-directlfq-Highlighted.svg", bbox_inches="tight", pad_inches=0.1)


# In[48]:


case_values = [{'condition': 'SWI'}, {'condition': 'Uninjured'}]

color_dict={'upregulated':condition_color['SWI'], 'downregulated': condition_color['Uninjured'],'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(3.356, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig4_swi_norm, values = case_values, pval = 0.05, return_df=True, color=color_dict, label=[0,0])
scplt.mark_volcano(ax, volcano_df, label=['C1qa','Itgam'], label_color='k')

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_FigS4-SWI_VOLCANO-directlfq-Highlighted.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


case_values = [{'condition': 'SWI'}, {'condition': 'Uninjured'}]

color_dict={'upregulated':condition_color['SWI'], 'downregulated': condition_color['Uninjured'],'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_fig4_swi_norm, values = case_values, pval = 0.05, return_df=True, color=color_dict, label=[0,0])

texts=[]
ax,t = scplt.mark_volcano(ax, volcano_df, label=['Aldh1l1','Slc1a3','Anxa3','Gfap','Vim'], label_color='blue', return_texts=True)
texts.extend(t)
ax,t = scplt.mark_volcano(ax, volcano_df, label=['Map2','Scn2a','Slc6a1'], label_color='purple', return_texts=True)
texts.extend(t)
ax,t = scplt.mark_volcano(ax, volcano_df, label=['C1qa','Itgam','Lgals1'], label_color='k', return_texts=True)
texts.extend(t)

scplt.volcano_adjust_and_outline_texts(texts, expand=(1.5, 3))

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig4-SWI_VOLCANO-Highlighted.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


text_kwargs=dict(
    fontsize=9,
    color='black',
    offset=0.35,         # vertical offset above anchor
)

fig, ax = pdata_fig4_swi_filtered.plot_abundance_boxgrid(namelist=['Cd68','Cd11b'], classes='condition', box=True, show_n=True, fig_width=2, fig_height=3,label_x=True, global_legend=False, palette=condition_color, text_kwargs=text_kwargs)

plt.savefig('MANUSCRIPT/MANUSCRIPT_Fig4-SWI_Abundance_CD68.svg')


# # Fig 6 - 6mo mouse

# ### pre-processing

# In[17]:


# pdata_6mo_diann_base = pAnnData.import_data(source_type='diann', report_file = 'data/2506_Mus1mo3mo6moAco/diann_6mo/report.parquet') #, obs_columns = obs_columns
pdata_6mo_diann_base = pAnnData.import_data(source_type='diann', report_file = 'data/2511_Mus6moSnpcONLY/report.parquet')


# In[18]:


colnames = [
    'date',
    'gradient',
    'sample_id',
    'size',
    'confirmation',
    'thickness',
    'sample',
    'organism',
    'region',
    'well_position'
]

# ignore 11-token parsing type on both
pdata_6mo_diann = pdata_6mo_diann_base.filter_sample(condition='parsingType == "10-tokens"')
pdata_6mo_diann.summary = scutils.parse_filename_index(pdata_6mo_diann.summary, colnames)

pdata_6mo_processed = pdata_6mo_diann.filter_prot_significant()


# In[19]:


pdata_6mo_snpc = pdata_6mo_processed.filter_sample(condition='region == "snpc"')
pdata_6mo_snpc = pdata_6mo_snpc.filter_sample(min_prot=1000)
pdata_6mo_snpc.summary["sample"] = pdata_6mo_snpc.summary["sample"].replace({
    "6mo-aggY": "Agg+",
    "6mo-aggN": "Agg-"
})
pdata_6mo_snpc.update_summary()

pdata_6mo_snpc_filter = pdata_6mo_snpc.copy()
pdata_6mo_snpc_filter = pdata_6mo_snpc_filter.filter_prot_found(min_ratio=0.4, group=['sample'],match_any=True)
pdata_6mo_snpc_filter = pdata_6mo_snpc_filter.filter_prot(valid_genes=True, unique_profiles=True)
pdata_6mo_snpc_filter = pdata_6mo_snpc_filter.filter_sample(min_prot=1000)
# pdata_6mo_snpc.impute(method='min', min_scale=0.2)

pdata_6mo_snpc_norm = pdata_6mo_snpc_filter.copy()
pdata_6mo_snpc_norm.normalize(method='directlfq')
# pdata_6mo_snpc_norm = pdata_6mo_snpc_norm.filter_sample(min_prot=1000)

pdata_6mo_snpc_norm.summary['sample'].value_counts()


# In[ ]:


pd.DataFrame(pdata_6mo_snpc_norm.summary)


# ### plots

# In[20]:


agg_color={
    "Agg+": "#C64D4A",
    "Agg-": "#BFBFBF"
}


# In[63]:


order = ['Agg+', 'Agg-']

dataset = pdata_6mo_snpc

fig, ax = plt.subplots(figsize=(2, 3))
sns.barplot(x='sample', y='protein_count',
            data=dataset.summary, errorbar='sd', capsize=.1,
            saturation=1, alpha=0.5, palette=agg_color,width=0.75,
            order=order, ax=ax, edgecolor='black',      # outline color
            linewidth=0.8)
sns.swarmplot(x='sample', y='protein_count',
                data=dataset.summary,order=order, color='k', ax=ax)

from scipy.stats import ttest_ind
aggN = dataset.summary[dataset.summary['sample'] == 'Agg-']['protein_count']
aggY = dataset.summary[dataset.summary['sample'] == 'Agg+']['protein_count']

max_prot_count = dataset.summary['protein_count'].max()

# scplt.plot_significance(ax, max_prot_count+350, 100,
#                         pval=ttest_ind(aggN, aggY).pvalue,
#                         fontsize=11)

plt.ylabel('Protein Count',fontsize=13, labelpad=4)
plt.xlabel('', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.ylim(0, 4200)
# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig5_protein-count.svg", bbox_inches="tight", pad_inches=0.1)


# In[22]:


plt.rcParams.update({'font.size': 12})

dataset = pdata_6mo_snpc
fig, ax = plt.subplots(figsize=(3, 3))
ax, venn_contents = scplt.plot_venn(ax, dataset,
                                    classes='sample',
                                    set_colors=["#C64D4A", "#BFBFBF"],
                                    return_contents=True,
                                    weighted=True, fixed_subset_sizes=(1, 1, 3),)

# Separate set labels
for t in [t for t in ax.texts if t.get_text() in {"Agg+", "Agg-"}]:
    x, y = t.get_position()
    if t.get_text() == "Agg+":
        t.set_position((x - 0.15, y - 0.02)); t.set_ha("right")
    else:
        t.set_position((x + 0.15, y - 0.02)); t.set_ha("left")

# plt.show()
plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig5_venn.svg", bbox_inches="tight", pad_inches=0.1)

setN = set(venn_contents["Agg+"])
setY = set(venn_contents["Agg-"])
venn_df = pd.DataFrame({
    "AggY_only": pd.Series(sorted(setN - setY)),
    "AggN_only": pd.Series(sorted(setY - setN)),
    "Both":      pd.Series(sorted(setN & setY))
})


# In[16]:


agg_color2={
    "Agg+": "#C64D4A",
    "Agg-": "#818181"
}

fig, ax = plt.subplots(1,1, figsize = (3.5,3.5))
ax = scplt.plot_pca(ax, pdata_6mo_snpc_norm, classes=['sample'], force=True, cmap=agg_color2)
scplt.shift_legend(ax)
# plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig5_PCA-2D.svg", bbox_inches="tight", pad_inches=0.1)

# make 3d figure 
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')
ax = scplt.plot_pca(ax, pdata_6mo_snpc_norm, classes=['sample'], force=True, plot_pc = [1,2,3], cmap=agg_color2)
zlab = ax.get_zlabel()
ax.set_zlabel("")
ax.text2D(1.1, 0.55, zlab,
          transform=ax.transAxes, rotation=90,
          ha='left', va='center', fontsize=9.5)
legend=ax.get_legend()
for patch in legend.legend_handles:
    patch.set_edgecolor('black')
if legend is not None:
    legend.remove()

plt.savefig(f"MANUSCRIPT/MANUSCRIPT_Fig5_PCA-3D.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


fig, ax = plt.subplots(1,1, figsize = (3,3))
umap_params={'min_dist': 0.1, 'n_neighbors': 5}
ax, umap = scplt.plot_umap(ax, pdata_6mo_snpc_norm, classes = 'sample', s=20, alpha=.8, force = True, umap_params=umap_params, return_fit=True, cmap=agg_color, )
scplt.shift_legend(ax)

# legend=ax.get_legend()
# for patch in legend.legend_handles:
#     patch.set_edgecolor('black')

# plt.savefig("MANUSCRIPT_Fig5_UMAP-impute.svg", bbox_inches="tight", pad_inches=0.1)


# ### de

# In[ ]:


case_values = [{'sample': 'Agg+'}, {'sample': 'Agg-'}]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_6mo_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[0,0])

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig5_VOLCANO-directlfq-empty.svg", bbox_inches="tight", pad_inches=0.1)
volcano_df.to_csv("MANUSCRIPT/MANUSCRIPT_Fig5_DE-directlfq.csv")


# In[82]:


case_values = [{'sample': 'Agg+'}, {'sample': 'Agg-'}]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_6mo_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', color=color_dict, label=[10,10])

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig5_VOLCANO-directlfq-top10.svg", bbox_inches="tight", pad_inches=0.1)


# In[71]:


def extract_gene_group(volcano_df, prefixes, top_n=5):
    """
    prefixes: list of lowercase prefixes, e.g. ['atp5']
    """
    mask_prefix = volcano_df["Genes"].str.lower().apply(
        lambda g: any(g.startswith(p) for p in prefixes)
    )
    sig_mask = (volcano_df["significance"] == 'upregulated') | \
               (volcano_df["significance"] == 'downregulated')

    gene_list = volcano_df.loc[mask_prefix & sig_mask, "Genes"].unique().tolist()

    if len(gene_list) == 0:
        return [], [], []

    top = (
        volcano_df[volcano_df["Genes"].isin(gene_list)]
        .sort_values("significance_score", ascending=True)
        .head(top_n)["Genes"]
        .tolist()
    )
    others = [g for g in gene_list if g not in top]

    return gene_list, top, others


# In[81]:


## RPL/RPS
rpl_list, rpl_top5, rpl_others = extract_gene_group(volcano_df, prefixes=["rpl"])
rps_list, rps_top5, rps_others = extract_gene_group(volcano_df, prefixes=["rps"])

print("rpl_list: ", len(rpl_list))
print("rps_list: ", len(rps_list))

case_values = [{'sample': 'Agg+'}, {'sample': 'Agg-'}]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}
rps_dict={'downregulated': '#5166FF'}
rpl_dict={'downregulated': '#1F2CCF'}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_6mo_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', no_marks=True)

texts = []
ax, t = scplt.mark_volcano(ax, volcano_df, label=rpl_top5, label_color='#1F2CCF',return_texts=True)
texts.extend(t)
ax, t = scplt.mark_volcano_by_significance(ax, volcano_df, label=rps_top5, color=rps_dict,return_texts=True)
texts.extend(t)
scplt.mark_volcano_by_significance(ax, volcano_df, label=rpl_others, color=rpl_dict, show_names=False)
scplt.mark_volcano_by_significance(ax, volcano_df, label=rps_others, color=rps_dict, show_names=False)

scplt.volcano_adjust_and_outline_texts(texts, expand=(2,2))

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

# plt.show()
plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig5_VOLCANO-directlfq-RplRps.svg", bbox_inches="tight", pad_inches=0.1)


# In[80]:


## Proteasome
proteasome_list, proteasome_top5, proteasome_others = extract_gene_group(volcano_df, prefixes=["psm"])
proteasome_list.append('Uchl1')

print("proteasome_list: ", len(proteasome_list))

case_values = [{'sample': 'Agg+'}, {'sample': 'Agg-'}]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_6mo_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', no_marks=True)

texts = []
ax, t = scplt.mark_volcano_by_significance(ax, volcano_df, label=proteasome_list, color=color_dict, return_texts=True)
texts.extend(t)
scplt.volcano_adjust_and_outline_texts(texts, expand=(2, 2))
# scplt.mark_volcano_by_significance(ax, volcano_df, label=proteasome_top5, color=color_dict)
# scplt.mark_volcano_by_significance(ax, volcano_df, label=proteasome_others, color=color_dict, show_names=False)

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig5_VOLCANO-directlfq-Proteasome.svg", bbox_inches="tight", pad_inches=0.1)


# In[76]:


## Mitochondria

mito_groups = {
    "complex5":   ["atp5"],
    "complex4":   ["cox"],
    "complex1":   ["nduf"],
    "complex2":   ["sdh"],
    "complex3":   ["uqcr"],
    "mito_other": ["cyc1", "cycs", "pgk", "vdac"],
}

mito_results = {}

for group_name, prefixes in mito_groups.items():
    gene_list, top5, others = extract_gene_group(volcano_df, prefixes)
    mito_results[group_name] = {
        "list": gene_list,
        "top5": top5,
        "others": others,
    }


# In[77]:


mito_color_dict = {
    "complex1": {     # NDUF*
        "upregulated":   "#D64545",   # strong red
        "downregulated": "#1E3A8A",   # deep indigo blue
    },
    "complex2": {     # SDH*
        "upregulated":   "#C03636",   # darker red
        "downregulated": "#1F4DA0",   # strong royal blue
    },
    "complex3": {     # UQCR*
        "upregulated":   "#B42E2E",   # deep cherry
        "downregulated": "#2665C2",   # vivid blue
    },
    "complex4": {     # COX*
        "upregulated":   "#A82828",   # wine red
        "downregulated": "#2D7AE6",   # bright azure
    },
    "complex5": {     # ATP5*
        "upregulated":   "#9E2020",   # maroon
        "downregulated": "#3D8BFF",   # clear blue
    },
    "mito_other": {
        "upregulated":   "#8C1A1A",   # dark brick red
        "downregulated": "#4C9CFF",   # lighter but still strong blue
    },
}


# In[78]:


## Mitochondria, all together

case_values = [{'sample': 'Agg+'}, {'sample': 'Agg-'}]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_6mo_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', no_marks=True)

texts = []
for group_name, d in mito_results.items():
    colors = mito_color_dict[group_name]

    if len(d["top5"]) > 0:
        # highlight top 5 (with names)
        ax, t = scplt.mark_volcano_by_significance(
            ax, volcano_df,
            label=d["top5"],
            color=colors,
            return_texts=True,
        )

        texts.extend(t)

    if len(d["others"]) > 0:
        # highlight the others (no names)
        scplt.mark_volcano_by_significance(
            ax, volcano_df,
            label=d["others"],
            color=colors,
            show_names=False,
        )

scplt.volcano_adjust_and_outline_texts(texts, expand=(2,2))

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig5_VOLCANO-directlfq-MitoAll.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


## Mitochondria, Ndufs separate, rest together
case_values = [{'sample': 'Agg+'}, {'sample': 'Agg-'}]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_6mo_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', no_marks=True)

d = mito_results["complex1"]
colors = mito_color_dict["complex1"]

if len(d["top5"]) > 0:
    # highlight top 5 (with names)
    scplt.mark_volcano_by_significance(
        ax, volcano_df,
        label=d["top5"],
        color=color_dict,
        show_names=True,
    )

if len(d["others"]) > 0:
    # highlight the others (no names)
    scplt.mark_volcano_by_significance(
        ax, volcano_df,
        label=d["others"],
        color=color_dict,
        show_names=False,
    )

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig5_VOLCANO-directlfq-MitoComplex1.svg", bbox_inches="tight", pad_inches=0.1)


# In[74]:


## RPL/RPS
kegg_list = [
    'H6pd', 'Psph', 'Inpp5f', 'Ndufb1', 'Ndufb11', 'Aldh6a1', 'Gmpr2', 'Aldh1l1',
    'Blvra', 'Pygl', 'Ggct', 'Cox6b1', 'Acat1', 'Pcca', 'Acadvl', 'Aco1', 'Aco2',
    'Impa2', 'Adsl', 'Adss2', 'Aldh1a1', 'Aldoc', 'Akr1b1', 'Alox12b', 'Auh', 'Bcat1', 'Cat', 'Ckb', 'Cmas',
    'Cox4i1', 'Cox6c', 'Dbt', 'Pcbd1', 'Eno2', 'Acsl1', 'Gad1', 'Gls',
    'Glud1', 'Got2', 'Gpx1', 'Gstm1', 'Gstm5', 'Gstz1', 'Hsd17b10', 'Hexa', 'Hexb',
    'Hsd17b4', 'Idh1', 'Idh3g', 'Ids', 'Inpp1', 'Itpa', 'Ldhb', 'Maoa', 'Mdh1',
    'mt-Atp8', 'mt-Co2', 'Ndufa2', 'Ndufa4', 'Ndufs4', 'Oat', 'Pafah1b2',
    'Pafah1b3', 'Pdha1', 'Pgk2', 'Prps1', 'Ptgs1', 'Abhd16a', 'Rpia', 'Spr',
    'Sucla2', 'Suox', 'Aldh5a1', 'Mars1', 'Acsl6', 'Ugdh', 'Uqcrq', 'Urod',
    'Pi4ka', 'Dlat', 'Aloxe3', 'Bpnt1', 'Bpnt2', 'Gnpda1', 'Abat', 'Asns', 'Pdhx',
    'Cdipt', 'Folh1', 'Atp6ap1', 'Nagk', 'Ivd', 'Alg2',
    'Hibadh', 'Smpd3', 'Ndufs5', 'Ndufb5', 'Ndufa9', 'Ndufa7', 'Asrgl1', 'Uqcrfs1',
    'Pccb', 'Oxct1', 'Ndufa6', 'Ndufb8', 'Ndufa10', 'Sdhb', 'Ndufc2',
    'Pdhb', 'Ndufb10', 'Ndufs3', 'Pycr2', 'Kdsr', 'Pgm2l1', 'Bdh1',
    'Oplah', 'Acyp2', 'Gns', 'Aacs', 'Sacm1l', 'Acsbg1'
]

print("kegg_list: ", len(kegg_list))

top10 = (
    volcano_df[volcano_df["Genes"].isin(kegg_list)]
    .assign(abs_score=lambda d: d["significance_score"].abs())
    .sort_values("abs_score", ascending=False)
    .head(10)
)

top10_kegg = top10["Genes"].tolist()

case_values = [{'sample': 'Agg+'}, {'sample': 'Agg-'}]

color_dict={'upregulated':"#C64D4A", 'downregulated': "#293DF5",'not_significant': "#FFFFFF6A"}

fig, ax = plt.subplots(figsize=(4, 4))
ax, volcano_df = scplt.plot_volcano(ax, pdata_6mo_snpc_norm,
                                    values=case_values,
                                    pval=0.05, return_df=True,
                                    fold_change_mode='mean', no_marks=True)

scplt.mark_volcano_by_significance(ax, volcano_df, label=top10_kegg, color=color_dict)

ax.set_ylabel('$log_{10}$ p value',fontsize=13, labelpad=4)
ax.set_xlabel('$log_{2}$ fold change', fontsize=13, labelpad=4)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)

plt.savefig("MANUSCRIPT/MANUSCRIPT_Fig5_VOLCANO-directlfq-KEGG.svg", bbox_inches="tight", pad_inches=0.1)


# In[ ]:


text_kwargs=dict(
    fontsize=11,
    color='black',
    offset=1,         # vertical offset above anchor
)

pdata_6mo_snpc_norm.plot_abundance_boxgrid(namelist=['Snca','Sncb','Sncg'], classes="sample", fig_width=2, fig_height=4, label_x=True, global_legend=True, box=True, show_n=True, text_kwargs=text_kwargs, y_max=11)


# In[ ]:


pdata_6mo_snpc_norm.annotate_found(on='peptide', classes='sample')


# In[ ]:


pdata_6mo_snpc_norm.pep.var


# ## check peptides

# In[ ]:


genes = ["Snca", "Sncb", "Sncg"]

idx = np.where(pdata_6mo_snpc_norm.prot.var["Genes"].isin(genes))[0]
idx


# In[ ]:


dense_rs = pdata_6mo_snpc_norm.rs[idx,:].toarray()


# In[ ]:


pep_mask = np.array(dense_rs.sum(axis=0)).ravel() > 0


# In[ ]:


pep_mask.sum()


# In[ ]:


dense_rs.sum()


# In[ ]:


dense_rs_filtered = dense_rs[:,pep_mask]


# In[ ]:


plt.figure(figsize=(3,3))

sns.heatmap(dense_rs_filtered, cmap='viridis', cbar=True, yticklabels=pdata_6mo_snpc_norm.prot.var["Genes"].iloc[idx].tolist(),)

plt.xlabel("Peptides")
plt.ylabel("Proteins")
plt.title("Protein–Peptide Mapping (RS Matrix)")
plt.show()


# In[ ]:


pdata_6mo_snpc_norm.rs[:,pep_mask].toarray()


# In[ ]:


# peptide sequence, number of samples it appears in, whether it's unique or not

