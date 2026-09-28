# Data directory for "Molecularly-guided spatial proteomics captures single-cell identity of the healthy and diseased nervous system"

This directory contains the processed input data required to run the analysis notebooks in the accompanying code repository.

**Manuscript:** Molecularly-guided spatial proteomics captures single-cell identity of the healthy and diseased nervous system  
*DOI: [10.1101/2025.02.10.637505](https://www.biorxiv.org/content/10.1101/2025.02.10.637505v2)*  
**Authors:**  
Sayan Dutta#, Marion Pang#, Gerard M. Coughlin, Sirisha Gudavalli,  
Michael L. Roukes, Tsui-Fen Chou, Viviana Gradinaru*

\# Equal contribution  
\* Corresponding author

## Data availability

Processed proteomics and transcriptomics reports (`.tsv`, `.csv`, `.parquet`) are organized by dataset in the folders listed below. Raw mass spectrometry files are not required to reproduce the figures and can be downloaded from:

- **MassIVE**: MSV000100704
- **ProteomeXchange**: PXD074008

At the time of this preprint and ongoing peer review, the processed proteomics reports are available under controlled access through the repository's reviewer mechanism. Researchers requiring access during this period may contact the corresponding authors to request the temporary access credentials. Full public access to all processed datasets will be enabled upon formal publication.

Transcriptomic reference data were obtained from the Allen Brain Atlas and Dropviz, as described in the manuscript.

## Directory structure

Each folder holds the processed report files for one dataset. File counts refer to the files deposited in each folder.

| Folder | Files | Figure(s) | Description |
|---|---|---|---|
| `2505_rna_prot_full` | 138 | Combined; 2, 4, 5 | Full multi-region dataset (`all_region_filtered`); see subsets below |
| `data2411_size` | 115 | 1 | Size series |
| `bulk-if` | 18 | 1 | Bulk – immunofluorescence (IF) |
| `bulk-hcr` | 6 | 1 | Bulk – hybridization chain reaction (HCR) |
| `bulk-ff` | 6 | 1 | Bulk – fixed-frozen |
| `2512_Ab-HCR` | 81 | 1 | Single-cell IF and HCR; see subsets below |
| `2511_ast-SWI-only` | 20 | 4 | Astrocyte SWI dataset; see subsets below |
| `2511_Mus6moSnpcONLY` | 36 | 6 | 6-month mouse SNpc |
| `2508_NHP` | 22 | 7 | Non-human primate |
| `2508_human` | 24 | 7 | Human |
| `2511_mus3mosnpc` | 28 | S15 | 3-month mouse SNpc |

## Subsets used in individual figures

Several figures use a filtered subset of a larger folder. Subsets are selected by filtering the files in the parent folder, as specified in the notebooks.

| Figure | Analysis | Parent folder | Files used |
|---|---|---|---|
| Combined | All regions | `2505_rna_prot_full` | 138 (all) |
| 2 | Cortex vs SNpc | `2505_rna_prot_full` | 100 |
| 4 | Astrocyte vs neuron | `2505_rna_prot_full` | 26 |
| 4 | SWI | `2511_ast-SWI-only` | 10 of 20 |
| 5 | ANXA+/- | `2505_rna_prot_full` | 33 |
| 1 | Single-cell IF | `2512_Ab-HCR` | 19 of 81 |
| 1 | Single-cell HCR | `2512_Ab-HCR` | 16 of 81 |

The bulk datasets in Figure 1 were processed with Proteome Discoverer (PD). All other datasets were processed with DIA-NN.

## Usage

Download the folders required for the figures of interest and place them in a single local directory. The path to this directory is set at the top of each notebook in the code repository and may need to be adjusted. See the code repository for the software environment and setup instructions.
