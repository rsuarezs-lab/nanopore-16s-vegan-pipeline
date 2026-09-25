# Community analysis of full-length 16S rRNA data (R)

Reproducible R workflow for the community-level analysis of full-length 16S rRNA
(Oxford Nanopore) microbiome data. It takes per-rank count tables and a metadata
table and produces publication-ready figures and statistics: composition,
alpha/beta diversity, ordination (PCoA), PERMANOVA, PERMDISP, envfit, per-group
alpha-diversity tests, and a PERMANOVA power analysis.

The input count tables are produced upstream by a separate bioinformatic
pipeline (read subsampling with seqkit, chimera removal with VSEARCH, and
taxonomic classification with EMU); that stage is **not** part of this
repository. This repository begins at the classified count tables.

```
        [ upstream: seqkit -> VSEARCH -> EMU ]      (separate pipeline)
                          |
                          v
   EMU per-rank count tables  +  sample map  +  metadata
                          |
                          v
   [ community_analysis.R ]   decontam -> filtering -> diversity -> vegan stats
                          |
                          v
              figures/   results/   data/
```

---

## 1. Requirements

| Software | Version | Notes |
|----------|---------|-------|
| R | 4.3.3 (tested) | |
| vegan | 2.6-4 | Bray-Curtis, PCoA, PERMANOVA, PERMDISP, envfit |
| decontam | 1.22.0 | reagent-contaminant identification (Bioconductor) |
| phyloseq | 1.46.0 | object handling for decontam (Bioconductor) |
| ape | — | phylogenetic/ordination utilities |
| readxl | — | reads the input workbook |
| pheatmap | — | heatmap |
| RColorBrewer | — | palette utilities |

Install:

```r
# CRAN
install.packages(c("vegan", "ape", "readxl", "pheatmap", "RColorBrewer"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "decontam"))
```

---

## 2. Input

A single Excel workbook (`.xlsx`) with the following sheets:

**Count tables** — one sheet per taxonomic rank, as produced by EMU
(`emu combine-outputs ... --counts`):

| Sheet | Content |
|-------|---------|
| `counts_phylum` | phylum-level counts |
| `counts_class` | class-level counts |
| `counts_order` | order-level counts |
| `counts_family` | family-level counts |
| `counts_genus` | genus-level counts |
| `counts_species` | species-level counts |

Each table has a taxonomy column named exactly as the rank (e.g. `genus`), and
one column per sequencing library, named `run_barcode` (e.g. `run1_barcode03`).

**`sample_map`** — maps each library to a sample:

| Column | Description |
|--------|-------------|
| `run` | sequencing run identifier |
| `sample` | sample identifier (negative controls must contain `(-)`) |
| `barcode` | barcode identifier |
| `Total reads` | (optional) per-library read count; used to pick the deepest run |

**`metadata`** — one row per sample: the first column is the sample ID (matching
`sample` above), followed by the grouping variables to test. Each grouping
variable should be a two-level categorical column.

A minimal example workbook is provided (`example_input.xlsx`) so the pipeline can
be run end-to-end without real data.

---

## 3. Usage

1. Place `community_analysis.R` and your input workbook in the same folder.
2. Edit the configuration block at the top of the script:

   ```r
   INPUT_FILE   <- "input.xlsx"              # your workbook
   GROUP_VARS   <- c("group_1", "group_2")   # metadata columns to test
   GROUP_LABELS <- c(group_1 = "Group 1", group_2 = "Group 2")
   CONT_VAR     <- NA                        # optional continuous covariate
   ```

3. Run:

   ```bash
   Rscript community_analysis.R
   ```

   or open it in RStudio (`Session -> Set Working Directory -> To Source File
   Location`) and **Source** it.

---

## 4. What the script does

| Step | Description |
|------|-------------|
| 1 | Read the sample map; for samples sequenced in more than one run, keep only the **deepest run** (runs are not summed). |
| 2 | Load the six count tables and consolidate to one column per sample. |
| 3 | Identify and remove **reagent contaminants** with decontam (prevalence method, threshold 0.1), using the negative controls. |
| 4 | Convert to relative abundance within each rank. |
| 5 | Load and align the metadata. |
| 6 | **Depth diagnostics**: rarefaction curves and Good's coverage. |
| 7 | Mean composition (family/genus/species panel + phylum). |
| 8 | Alpha diversity (observed richness, Shannon, Simpson). |
| 9 | **Bray-Curtis** dissimilarity and **PCoA**. |
| 10 | Heatmap of the top-25 genera (UPGMA on Bray-Curtis). |
| 11 | Per-sample stacked bars. |
| 12 | **PERMANOVA** (`adonis2`) and **PERMDISP** (`betadisper`), 999 permutations. |
| 13 | **envfit** (vector fitting onto the ordination). |
| 14 | Per-group alpha-diversity tests (Mann-Whitney). |
| 15 | **PERMANOVA power analysis** by simulation across a range of effect sizes. |

All community-level statistics use genus-level Bray-Curtis dissimilarities.
Figures use a colorblind-safe scheme (Okabe-Ito for categories, viridis for
gradients).

---

## 5. Output

```
figures/    composition_family_genus_species, composition_phylum,
            alpha_diversity, heatmap_genera, stacked_bars_genus,
            pcoa_envfit, rarefaction_curves, permanova_power
            (each as PDF + TIFF + PNG, 600 dpi)

results/    permanova_permdisp.csv, envfit.csv, cohort.csv,
            depth_diagnostics.csv, alpha_diversity.csv,
            alpha_diversity_tests_by_group.csv, permanova_power.csv, braycurtis_matrix.csv, pcoa_coordinates.csv,
            decontam_result.csv

data/       genus_relative_abundance_pct.csv, genus_counts.csv,
            metadata_aligned.csv
```

---

## 6. Parameter summary

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `DECONTAM_THRESHOLD` | 0.1 | decontam prevalence threshold (conservative) |
| `N_PERM` | 999 | permutations (PERMANOVA / PERMDISP / envfit) |
| `POWER_NSIM` | 200 | simulations per effect size (power analysis) |
| Significance level | p < 0.05 | |

---

## 7. Notes

- Bray-Curtis dissimilarities are computed on relative abundances and are
  scale-invariant; sequencing-depth adequacy is documented with rarefaction
  curves and Good's coverage rather than by rarefying the data.
- Observed richness is sensitive to sequencing depth; alpha-diversity indices are
  reported alongside the depth diagnostics.
- Species-level assignments from full-length 16S should be treated as
  presumptive.
- PERMANOVA has no closed-form power; power is therefore estimated by simulation.

---

## 8. Citations

- Oksanen J, et al. (2022). vegan: Community Ecology Package (R package v2.6-4).
- Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ (2018). Simple
  statistical identification and removal of contaminant sequences in marker-gene
  and metagenomics data. *Microbiome* 6:226.
- Anderson MJ (2001). A new method for non-parametric multivariate analysis of
  variance. *Austral Ecology* 26:32-46.
- Anderson MJ (2006). Distance-based tests for homogeneity of multivariate
  dispersions. *Biometrics* 62:245-253.
- R Core Team (2024). R: A Language and Environment for Statistical Computing.

---

## 9. License

MIT.
