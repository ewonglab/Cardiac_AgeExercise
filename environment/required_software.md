# Reported software requirements

The manuscript reports the following analysis environment. This is a requirements record, not a lockfile; a fully executable environment cannot be constructed from the supplied sources alone because several external tools, models, and input datasets are not included.

| Component | Reported version or setting |
| --- | --- |
| R | 4.2.1 |
| Cell Ranger | 7.0.0 (snRNA-seq) |
| Cell Ranger ATAC | 1.2.0 |
| CellBender remove-background | 0.3.0 |
| Scrublet | 0.2.3 |
| Seurat | 4.4.0 |
| Signac | 1.15.0 |
| scDblFinder | 1.13.8 |
| MAST | 1.24.1 |
| edgeR | 3.40.2 |
| DESeq2 | 1.38.3 |
| clusterProfiler | 4.6.2 |
| CellChat | 1.6.1 |
| glmnet | 4.1.4 |
| TOBIAS | 0.17.0 |
| MACS3 | 3.0.0b1 |
| Cicero | 1.3.9 |
| Python | 3.12 (echocardiography analysis reported) |
| AlphaGenome MPRA workflow | `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `torch`, `biopython`, `pyjaspar` and `alphagenome-pytorch`; see `requirements_mpra_seqlet_grammar.txt` |

Additional R packages used by the selectively included plotting and module scripts include `BSgenome.Mmusculus.UCSC.mm10`, `GenomicRanges`, `VennDiagram`, `patchwork` and `scales`. Required command-line tools include nf-core/cutandrun, SEACR, HOMER, FIMO, and the TooManyCells Singularity image. The AlphaGenome MPRA implementation is included in `scripts/mpra_seqlet_grammar/`, but a runnable release still requires the locally supplied `alphagenome_encoder_ft` wheel, frozen encoder checkpoint, fold-specific MPRA heads, embeddings, motif panel and input data.
