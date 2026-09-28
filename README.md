# Code repository for "Molecularly-guided spatial proteomics captures single-cell identity of the healthy and diseased nervous system"

This repository contains code associated with the manuscript:

**Molecularly-guided spatial proteomics captures single-cell identity of the healthy and diseased nervous system**  
*DOI: [10.1101/2025.02.10.637505](https://www.biorxiv.org/content/10.1101/2025.02.10.637505v2)*  
**Authors:**  
Sayan Dutta#, Marion Pang#, Gerard M. Coughlin, Sirisha Gudavalli,  
Michael L. Roukes, Tsui-Fen Chou, Viviana Gradinaru*

\# Equal contribution  
\* Corresponding author



This repository provides Python analysis code used to generate figures in the manuscript.  
The analyses rely on processed proteomics and transcriptomics data deposited in public repositories, together with the authors’ open-source analysis package [**scpviz**](https://github.com/gnaprs/scpviz).

## Contents

- `analysis_scSP_manuscript.ipynb`  
  Generates proteomics-based figures and analyses

- `analysis_scSP_manuscript_transcriptomics.ipynb`  
  Transcriptomics-informed analyses and cross-omics comparisons.

- `data/`  
  Directory containing processed input files required to run the notebooks. *At the time of this preprint and ongoing peer review, the processed proteomics reports are available under controlled access through the repository’s reviewer mechanism. Researchers requiring access during this period may contact the corresponding authors to request the temporary access credentials.*

Paths to downloaded data may need to be adjusted at the top of each notebook.

## Data availability

Transcriptomic reference data were obtained from the Allen Brain Atlas and Dropviz, as described in the manuscript.
Processed proteomics and transcriptomics reports (`.tsv`, `.csv`, `.parquet`) are found in the `data/` folder of this repository.

Raw data can be downloaded from:

- **MassIVE**: MSV000100704  
- **ProteomeXchange**: PXD074008  

For running the analysis within this directory, only processed reports found in `data/` are required; the raw MS files are not needed to reproduce figures. All report files should be placed in a local directory that will be referenced in the notebooks (paths are documented in the notebooks).

At the time of this preprint and ongoing peer review, the processed proteomics reports are available under controlled access through the repository’s reviewer mechanism. Researchers requiring access during this period may contact the corresponding authors to request the temporary access credentials.

Full public access to all processed datasets will be enabled upon formal publication.

## Setup

Analyses were performed using **Python 3.11**.  
The complete software environment, including all dependencies and exact package versions, is provided in `environment.yml`.

To recreate the environment:

```bash
conda env create -f environment.yml
conda activate scpviz
```

This environment includes the authors' open-source analysis package [`scpviz`](https://github.com/gnaprs/scpviz) and all additional dependencies required to run the notebooks.
