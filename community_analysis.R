#!/usr/bin/env Rscript
# ==============================================================================
#  community_analysis.R
#  Community-level analysis of full-length 16S rRNA (EMU) count tables
#  Bray-Curtis / PCoA / PERMANOVA / PERMDISP / envfit in R (vegan, decontam)
# ==============================================================================
#
#  Downstream analysis for EMU-classified full-length 16S data. Takes the six
#  per-rank EMU count tables plus a sample map and a metadata table, and runs:
#    - replicate consolidation (keeping the deepest run per sample)
#    - reagent-contaminant removal with decontam (negative controls)
#    - sequencing-depth diagnostics (rarefaction curves + Good's coverage)
#    - composition, alpha diversity, Bray-Curtis / PCoA
#    - PERMANOVA, PERMDISP, envfit, per-group alpha-diversity tests,
#      and a PERMANOVA power analysis by simulation
#  producing publication-ready figures (600 dpi) and result tables.
#
#  ---------------------------------------------------------------------------
#  INPUT  (single Excel workbook, see INPUT_FILE) with these sheets:
#
#    counts_phylum, counts_class, counts_order, counts_family,
#    counts_genus, counts_species
#        EMU combined count tables (one taxonomy column named as the rank, e.g.
#        "genus"; the remaining columns are libraries, named run + "_" + barcode).
#
#    sample_map      columns:  run | sample | barcode
#        Maps each library to a sample ID. Libraries of the same sample across
#        runs are consolidated (deepest run kept). Negative controls must contain
#        "(-)" in the `sample` field. A "Total reads" column (per-library depth)
#        is used to pick the deepest run; if absent, the first run is kept.
#
#    metadata        one row per sample (first column = sample ID) with the
#        grouping variables to test. Edit GROUP_VARS / CONT_VAR below.
#
#  OUTPUT
#    figures/   composition, alpha diversity, heatmap, stacked bars,
#               PCoA+envfit, rarefaction, PERMANOVA power  (PDF + TIFF/PNG 600dpi)
#    results/   PERMANOVA/PERMDISP, envfit, cohort, depth diagnostics,
#               alpha-diversity tests, power analysis, matrices
#    data/      relative-abundance and count matrices, aligned metadata
#  ---------------------------------------------------------------------------
#  DEPENDENCIES
#    CRAN:         vegan, ape, readxl, pheatmap, RColorBrewer
#    Bioconductor: phyloseq, decontam
#
#    install.packages(c("vegan","ape","readxl","pheatmap","RColorBrewer"))
#    if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#    BiocManager::install(c("phyloseq","decontam"))
#
#  Tested with R 4.3.3, vegan 2.6-4, decontam 1.22.0, phyloseq 1.46.0.
#  License: MIT.
# ==============================================================================

# ── 0. User configuration ─────────────────────────────────────────────────────
INPUT_FILE <- "input.xlsx"     # input workbook (see header)

# Grouping variables to test (must match column names in the `metadata` sheet).
# Each is treated as a two-level categorical variable.
GROUP_VARS <- c("group_1", "group_2", "group_3")

# Optional display labels for figures/tables (named by the GROUP_VARS above).
GROUP_LABELS <- c(
  group_1 = "Group 1",
  group_2 = "Group 2",
  group_3 = "Group 3"
)

# Optional continuous covariate in `metadata` to fit onto the ordination
# (set to NA to skip).
CONT_VAR <- NA

# Analysis parameters
N_PERM             <- 999   # permutations for PERMANOVA / PERMDISP / envfit
DECONTAM_THRESHOLD <- 0.1   # decontam prevalence threshold (0.1 = conservative)
SUBSAMPLE_CAP      <- 50000  # per-library depth cap used upstream (seqkit);
                            #   informs the deepest-run rule only.
POWER_NSIM         <- 200   # simulations per effect size in the power analysis

# ── Packages ──────────────────────────────────────────────────────────────────
.need <- function(pkgs, bioc = FALSE) for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(p, update = FALSE, ask = FALSE)
  } else install.packages(p, repos = "https://cloud.r-project.org")
}
.need(c("vegan","ape","readxl","pheatmap","RColorBrewer"))
.need(c("phyloseq","decontam"), bioc = TRUE)
suppressWarnings(suppressPackageStartupMessages({
  library(vegan); library(ape); library(readxl)
  library(pheatmap); library(RColorBrewer)
  library(phyloseq); library(decontam)
}))

# UTF-8 locale so any accented labels render correctly across platforms.
suppressWarnings(for (loc in c("C.UTF-8","en_US.UTF-8","C.utf8"))
  if (Sys.setlocale("LC_CTYPE", loc) != "") break)

set.seed(42)
for (d in c("figures","results","data")) dir.create(d, showWarnings = FALSE)

# ── Colorblind-safe palettes (Okabe-Ito qualitative + viridis sequential) ─────
OKABE_ITO <- c("#E69F00","#56B4E9","#009E73","#F0E442",
               "#0072B2","#D55E00","#CC79A7","#000000")
cb_palette <- function(n) {
  if (n <= length(OKABE_ITO)) OKABE_ITO[seq_len(n)] else colorRampPalette(OKABE_ITO)(n)
}
viridis_cb <- colorRampPalette(c("#440154","#3b528b","#21908c","#5dc863","#fde725"))

# ── Rank -> sheet / taxonomy-column mapping ───────────────────────────────────
LEVEL_SHEETS <- c(Phylum="counts_phylum", Class="counts_class", Order="counts_order",
                  Family="counts_family", Genus="counts_genus", Species="counts_species")
LEVEL_COLS   <- c(Phylum="phylum", Class="class", Order="order",
                  Family="family", Genus="genus", Species="species")

# ── Helpers ───────────────────────────────────────────────────────────────────
read_sheet <- function(sheet) as.data.frame(read_excel(INPUT_FILE, sheet = sheet))

# Natural (numeric-aware) ordering of sample IDs.
order_natural <- function(x) x[order(suppressWarnings(as.numeric(gsub("\\D","",x))), x)]

# Save a base-graphics figure in three publication formats (PDF + TIFF/PNG 600dpi).
has_cairo <- capabilities("cairo")
save_fig <- function(draw, name, w = 8, h = 6.5) {
  for (fmt in c("pdf","tiff","png")) {
    fn <- file.path("figures", paste0(name, ".", fmt))
    if (fmt=="pdf")  { if (has_cairo) cairo_pdf(fn,width=w,height=h) else pdf(fn,width=w,height=h) }
    if (fmt=="tiff") tiff(fn,width=w,height=h,units="in",res=600,compression="lzw",
                          type=if(has_cairo)"cairo" else NULL)
    if (fmt=="png")  png(fn,width=w,height=h,units="in",res=600,type=if(has_cairo)"cairo" else NULL)
    draw(); dev.off()
  }
  message("  figure: ", name)
}

message("== Community-level 16S analysis ==")

# ── 1. Sample map + deepest-run consolidation rule ────────────────────────────
message("[1] Reading sample map ...")
smap <- read_sheet("sample_map")
smap$library <- paste0(smap$run, "_", smap$barcode)
lib2sample   <- setNames(smap$sample, smap$library)

# For samples sequenced in more than one run, keep ONLY the deepest run (runs are
# not summed), so each sample derives from a single library. Depth comes from the
# "Total reads" column if present; otherwise the first library is kept.
depth_col <- grep("^total.?reads$", names(smap), ignore.case = TRUE, value = TRUE)[1]
smap$.depth <- if (!is.na(depth_col)) suppressWarnings(as.numeric(smap[[depth_col]])) else NA
drop_libs <- character(0)
real_samples <- unique(smap$sample[!grepl("\\(-\\)", smap$sample)])
for (s in real_samples) {
  rows <- smap[smap$sample == s, , drop = FALSE]
  if (nrow(rows) <= 1) next
  keep <- if (all(is.na(rows$.depth))) 1 else which.max(rows$.depth)
  drop_libs <- c(drop_libs, rows$library[-keep])
}
lib2sample_eff <- lib2sample
lib2sample_eff[drop_libs] <- NA
message("    libraries: ", nrow(smap),
        " | negative controls: ", sum(grepl("\\(-\\)", smap$sample)),
        " | dropped (shallower run of a resequenced sample): ", length(drop_libs))

# ── 2. Load count tables and consolidate to one column per sample ─────────────
message("[2] Loading count tables ...")
load_level <- function(sheet, level_col, keep_controls = TRUE) {
  df <- read_sheet(sheet)
  df <- df[!is.na(df[[level_col]]), ]                 # drop rows unassigned at rank
  lib_cols <- intersect(names(df), names(lib2sample))
  for (c in lib_cols) df[[c]] <- as.numeric(df[[c]])
  df[is.na(df)] <- 0
  mat <- rowsum(df[, lib_cols, drop=FALSE], group = df[[level_col]])  # taxa x libs
  mat <- t(as.matrix(mat))                            # libraries x taxa
  samp <- lib2sample_eff[rownames(mat)]               # effective map (deepest run)
  mat  <- mat[!is.na(samp), , drop=FALSE]; samp <- samp[!is.na(samp)]
  mat  <- rowsum(mat, group = samp)                   # one row per sample
  if (!keep_controls) mat <- mat[!grepl("\\(-\\)", rownames(mat)), , drop=FALSE]
  mat
}
LEVELS <- lapply(names(LEVEL_SHEETS), function(L)
  load_level(LEVEL_SHEETS[[L]], LEVEL_COLS[[L]], keep_controls = TRUE))
names(LEVELS) <- names(LEVEL_SHEETS)

# ── 3. decontam: reagent contaminants from negative controls (genus level) ────
message("[3] decontam: identifying reagent contaminants ...")
genus_with_ctrl <- LEVELS[["Genus"]]
ps <- phyloseq(
  otu_table(t(as.matrix(genus_with_ctrl)), taxa_are_rows = TRUE),
  sample_data(data.frame(is_neg = grepl("\\(-\\)", rownames(genus_with_ctrl)),
                         row.names = rownames(genus_with_ctrl)))
)
contam <- isContaminant(ps, method = "prevalence", neg = "is_neg",
                        threshold = DECONTAM_THRESHOLD)
contaminant_taxa <- rownames(contam)[contam$contaminant]
message("    contaminant genera removed: ",
        if (length(contaminant_taxa)) paste(contaminant_taxa, collapse=", ") else "none")
write.csv(contam, "results/decontam_result.csv")

# Samples = everything except negative controls.
samples <- order_natural(rownames(genus_with_ctrl)[!grepl("\\(-\\)", rownames(genus_with_ctrl))])
message("    samples retained: ", length(samples))

# ── 4. Build per-rank relative-abundance matrices ─────────────────────────────
# (no sample-level abundance filter here; add project-specific filters if needed)
message("[4] Building relative-abundance matrices ...")
REL <- list(); COUNTS <- list()
for (L in names(LEVELS)) {
  m <- LEVELS[[L]][samples, , drop=FALSE]
  m <- m[, !is.na(colnames(m)) & trimws(colnames(m)) != "", drop=FALSE]  # drop empty names
  if (L == "Genus") m <- m[, !colnames(m) %in% contaminant_taxa, drop=FALSE]
  m <- m[, colSums(m) > 0, drop=FALSE]
  COUNTS[[L]] <- m
  REL[[L]]    <- m / rowSums(m)
}
rel_g <- REL[["Genus"]]; counts_g <- COUNTS[["Genus"]]
write.csv(round(rel_g*100, 4), "data/genus_relative_abundance_pct.csv")
write.csv(counts_g,            "data/genus_counts.csv")

# Cohort summary
write.csv(data.frame(sample = samples, reads = rowSums(counts_g[samples,,drop=FALSE]),
                     genera = rowSums(rel_g[samples,,drop=FALSE] > 0)),
          "results/cohort.csv", row.names = FALSE)

# ── 5. Metadata ───────────────────────────────────────────────────────────────
message("[5] Loading metadata ...")
meta_all <- read_sheet("metadata")
rownames(meta_all) <- meta_all[[1]]; meta_all[[1]] <- NULL
metadata <- meta_all[samples, , drop = FALSE]
for (col in names(metadata)) if (is.character(metadata[[col]])) Encoding(metadata[[col]]) <- "UTF-8"
write.csv(metadata, "data/metadata_aligned.csv")

# ── 6. Sequencing-depth diagnostics (rarefaction + Good's coverage) ───────────
message("[6] Depth diagnostics (rarefaction) ...")
cnt <- counts_g; storage.mode(cnt) <- "integer"; depth <- rowSums(cnt)
save_fig(function() {
  par(mar=c(4.5,4.5,3,1))
  cols <- viridis_cb(nrow(cnt))[rank(depth, ties.method="first")]
  rarecurve(cnt, step=500, col=cols, lwd=1.8, label=FALSE,
            xlab="Number of reads (sequencing depth)",
            ylab="Observed genera (richness)",
            main=sprintf("Rarefaction curves - genus level (n=%d)", nrow(cnt)))
  low <- names(sort(depth)[1:min(3,length(depth))])
  for (b in low) text(depth[b], sum(cnt[b,]>0), labels=b, pos=4, cex=0.7, col="gray30")
}, "rarefaction_curves", w=9, h=6.5)
goods <- sapply(rownames(cnt), function(s){ x <- cnt[s,]; 1 - sum(x==1)/sum(x) })
write.csv(data.frame(sample=rownames(cnt), reads=depth, genera=rowSums(cnt>0),
                     goods_coverage_pct=round(goods*100,2)),
          "results/depth_diagnostics.csv", row.names=FALSE)

# ── 7. Mean composition (combined family/genus/species panel + phylum) ────────
message("[7] Composition figures ...")
bar_panel <- function(rel, tag, subtitle, top_n=12, italic=FALSE) {
  m <- sort(colMeans(rel), decreasing=TRUE)
  top <- m[seq_len(min(top_n,length(m)))]; rest <- sum(m[-seq_len(min(top_n,length(m)))])
  if (rest>0) top <- c(top, Other=rest)
  top <- rev(top); v <- as.numeric(top)*100
  bp <- barplot(v, horiz=TRUE, names.arg=names(top), las=1, col=cb_palette(length(top)),
                border="gray30", xlim=c(0,max(v)*1.18), xlab="Mean relative abundance (%)",
                cex.names=0.68, font.axis=ifelse(italic,3,1))
  text(v, bp, labels=sprintf("%.1f%%",v), pos=4, cex=0.6, xpd=TRUE)
  mtext(tag, side=3, adj=0, line=1, font=2, cex=1.1); mtext(subtitle, side=3, adj=0.5, line=1, cex=0.85)
}
save_fig(function() {
  par(mfrow=c(1,3), mar=c(4.5,9,3,2), oma=c(0,0,2,0))
  bar_panel(REL[["Family"]],  "a", "Family",  12, FALSE)
  bar_panel(REL[["Genus"]],   "b", "Genus",   12, FALSE)
  bar_panel(REL[["Species"]], "c", "Species", 12, TRUE)
  mtext(sprintf("Mean community composition (n=%d)", nrow(rel_g)), outer=TRUE, cex=1.05, font=2)
}, "composition_family_genus_species", w=15, h=5.5)
save_fig(function() {
  par(mar=c(4.5,12,3,2)); rel <- REL[["Phylum"]]
  top <- rev(sort(colMeans(rel), decreasing=TRUE)); v <- as.numeric(top)*100
  bp <- barplot(v, horiz=TRUE, names.arg=names(top), las=1, col=cb_palette(length(top)),
                border="gray30", xlim=c(0,max(v)*1.15), xlab="Mean relative abundance (%)",
                main=sprintf("Mean composition - Phylum (n=%d)", nrow(rel)))
  text(v, bp, labels=sprintf("%.1f%%",v), pos=4, cex=0.7, xpd=TRUE)
}, "composition_phylum", w=8, h=5.5)

# ── 8. Alpha diversity ────────────────────────────────────────────────────────
message("[8] Alpha diversity ...")
alpha <- data.frame(Observed=rowSums(rel_g>0),
                    Shannon=diversity(rel_g,"shannon"),
                    Simpson=diversity(rel_g,"simpson"))
write.csv(round(alpha,4), "results/alpha_diversity.csv")
save_fig(function() {
  par(mfrow=c(1,3), mar=c(3,4.5,3,1))
  ms <- c("Observed","Shannon","Simpson"); tt <- c("Observed richness","Shannon index","Simpson index (1-D)")
  cl <- c("#9ecae1","#a1d99b","#fcae91")
  for (i in seq_along(ms)) {
    boxplot(alpha[[ms[i]]], col=cl[i], ylab="Index value", main=tt[i],
            outline=FALSE, boxwex=0.5, ylim=range(alpha[[ms[i]]])*c(0.9,1.05))
    set.seed(1); points(jitter(rep(1,nrow(alpha)),amount=0.08), alpha[[ms[i]]],
                        pch=21, bg=cl[i], col="black", cex=1)
    md <- median(alpha[[ms[i]]]); abline(h=md, lty=2, col="gray40")
    text(1.35, md, sprintf("median = %.2f", md), pos=3, cex=0.85, col="gray20")
  }
}, "alpha_diversity", w=11, h=4)

# ── 9. Bray-Curtis + PCoA ─────────────────────────────────────────────────────
message("[9] Bray-Curtis + PCoA ...")
bc <- vegdist(rel_g, method="bray")
pcoa <- cmdscale(bc, k=nrow(rel_g)-1, eig=TRUE)
ev <- pcoa$eig[pcoa$eig>0]; var_exp <- ev/sum(ev)*100
write.csv(round(as.matrix(bc),6), "results/braycurtis_matrix.csv")
write.csv(data.frame(sample=rownames(pcoa$points),
                     PCo1=round(pcoa$points[,1],6), PCo2=round(pcoa$points[,2],6)),
          "results/pcoa_coordinates.csv", row.names=FALSE)

# ── 10. Heatmap of top-25 genera (UPGMA on Bray-Curtis) ───────────────────────
message("[10] Heatmap ...")
n_top <- min(25, ncol(rel_g))
top_idx <- order(colMeans(rel_g), decreasing=TRUE)[seq_len(n_top)]
sub_mat <- rel_g[, top_idx, drop=FALSE]
hc_s <- hclust(vegdist(rel_g,"bray"), method="average")
hc_t <- hclust(vegdist(t(sub_mat),"bray"), method="average")
dom  <- colnames(rel_g)[apply(rel_g,1,which.max)]; names(dom) <- rownames(rel_g)
topd <- names(sort(table(dom), decreasing=TRUE))[1:min(5,length(unique(dom)))]
grp  <- ifelse(dom %in% topd, dom, "Other")
ann  <- data.frame(Dominant=grp, row.names=names(dom))
annc <- list(Dominant=setNames(cb_palette(length(unique(grp))), unique(grp)))
save_fig(function() {
  pheatmap(t(log10(sub_mat + 1e-4)), cluster_rows=hc_t, cluster_cols=hc_s,
           color=viridis_cb(100), annotation_col=ann, annotation_colors=annc,
           fontsize_row=8, fontsize_col=8, border_color=NA,
           main="Top-25 genera (Bray-Curtis UPGMA)")
}, "heatmap_genera", w=13, h=8)

# ── 11. Per-sample stacked bars ───────────────────────────────────────────────
message("[11] Stacked bars ...")
save_fig(function(top_n=12) {
  m <- sort(colMeans(rel_g), decreasing=TRUE)
  keep <- names(m)[seq_len(min(top_n,length(m)))]; other <- setdiff(colnames(rel_g), keep)
  tt <- rel_g[, keep, drop=FALSE]
  if (length(other)>0) tt <- cbind(tt, Other=rowSums(rel_g[,other,drop=FALSE]))
  tt <- tt[order_natural(rownames(tt)), ]; M <- t(as.matrix(tt))
  cols <- c(cb_palette(nrow(M)-1), "gray70")
  par(mar=c(6,4.5,3,12), xpd=NA)
  barplot(M, col=cols, border=NA, las=2, cex.names=0.7, ylab="Relative abundance",
          ylim=c(0,1), main="Per-sample composition - Genus")
  usr <- par("usr")
  legend(x=usr[2]*1.02, y=usr[4], legend=rownames(M), fill=cols, border=NA,
         bty="n", cex=0.75, title="Genus", xjust=0)
}, "stacked_bars_genus", w=15, h=6.5)

# ── 12. PERMANOVA + PERMDISP per grouping variable ────────────────────────────
message("[12] PERMANOVA + PERMDISP ...")
res_tab <- data.frame()
for (v in GROUP_VARS) {
  if (!(v %in% names(metadata))) next
  idx <- !is.na(metadata[[v]]) & metadata[[v]]!=""
  if (sum(idx) < 4 || length(unique(metadata[[v]][idx])) != 2) next
  d  <- as.dist(as.matrix(bc)[idx,idx])
  ad <- adonis2(d ~ metadata[[v]][idx], permutations=N_PERM)
  bd <- betadisper(d, factor(metadata[[v]][idx]))
  pt <- permutest(bd, permutations=N_PERM)
  res_tab <- rbind(res_tab, data.frame(
    variable = if (v %in% names(GROUP_LABELS)) GROUP_LABELS[[v]] else v,
    n = sum(idx), R2 = round(ad$R2[1],3),
    PERMANOVA_p = ad[["Pr(>F)"]][1], PERMDISP_p = pt$tab$`Pr(>F)`[1]))
}
write.csv(res_tab, "results/permanova_permdisp.csv", row.names=FALSE)
print(res_tab, row.names=FALSE)

# ── 13. envfit + PCoA figure ──────────────────────────────────────────────────
message("[13] envfit + PCoA ...")
env_df <- data.frame(row.names=rownames(metadata))
if (!is.na(CONT_VAR) && CONT_VAR %in% names(metadata))
  env_df[[CONT_VAR]] <- as.numeric(metadata[[CONT_VAR]])
for (v in GROUP_VARS) if (v %in% names(metadata)) {
  lv <- sort(unique(na.omit(metadata[[v]])))
  if (length(lv)==2) env_df[[v]] <- as.integer(metadata[[v]] == lv[2])
}
ord <- pcoa$points[,1:2]
ef  <- if (ncol(env_df)) envfit(ord, env_df, permutations=N_PERM, na.rm=TRUE) else NULL
if (!is.null(ef)) write.csv(data.frame(variable=names(ef$vectors$r), R2=round(ef$vectors$r,4),
                     p=ef$vectors$pvals, row.names=NULL), "results/envfit.csv", row.names=FALSE)
save_fig(function() {
  tab <- sort(table(dom), decreasing=TRUE); pal <- setNames(cb_palette(length(tab)), names(tab))
  par(mar=c(4.5,4.5,3,2))
  plot(ord[,1], ord[,2], type="n",
       xlab=sprintf("PCo1 (%.1f%%)",var_exp[1]), ylab=sprintf("PCo2 (%.1f%%)",var_exp[2]),
       main="PCoA (Bray-Curtis, genus) + envfit")
  abline(h=0,v=0,col="gray70",lty=3)
  for (g in names(tab)) { ii <- dom==g; points(ord[ii,1],ord[ii,2],pch=21,bg=pal[g],col="black",cex=1.6) }
  if (!is.null(ef)) plot(ef, p.max=1, col=ifelse(ef$vectors$pvals<0.05,"#c0392b","#5f6a6a"), cex=0.85)
  legend("topleft", legend=paste0(names(tab)," (n=",as.integer(tab),")"),
         pt.bg=pal, pch=21, bty="n", title="Dominant genus", cex=0.85)
}, "pcoa_envfit", w=9, h=7.5)

# ── 14. Alpha-diversity tests by group (Mann-Whitney) ─────────────────────────
message("[14] Alpha-diversity tests by group ...")
alpha_tests <- data.frame()
for (v in GROUP_VARS) {
  if (!(v %in% names(metadata))) next
  idx <- !is.na(metadata[[v]]) & metadata[[v]]!=""
  g <- factor(metadata[[v]][idx]); if (nlevels(g)!=2) next
  for (metric in c("Observed","Shannon","Simpson")) {
    x <- alpha[[metric]][idx]; wt <- suppressWarnings(wilcox.test(x ~ g))
    alpha_tests <- rbind(alpha_tests, data.frame(
      variable = if (v %in% names(GROUP_LABELS)) GROUP_LABELS[[v]] else v,
      metric = metric,
      median_1 = round(median(x[g==levels(g)[1]]),2),
      median_2 = round(median(x[g==levels(g)[2]]),2),
      W = unname(wt$statistic), p_value = round(wt$p.value,4)))
  }
}
if (nrow(alpha_tests)) write.csv(alpha_tests, "results/alpha_diversity_tests_by_group.csv", row.names=FALSE)

# ── 15. PERMANOVA power analysis by simulation ────────────────────────────────
# PERMANOVA has no closed-form power. It is estimated by simulation: starting
# from the observed community pool (preserving natural dispersion), a location
# shift of increasing magnitude is induced between two groups, and the fraction
# of simulations reaching p<0.05 is recorded together with the achieved R^2.
message("[15] PERMANOVA power analysis (simulation) ...")
sim_power <- function(rel_base, n1, n2, shift, nsim=POWER_NSIM, nperm=199) {
  n <- n1+n2; mu <- colMeans(rel_base)
  resp <- order(mu, decreasing=TRUE)[seq_len(max(2, round(ncol(rel_base)*0.3)))]
  det <- 0; r2 <- numeric(nsim)
  for (s in seq_len(nsim)) {
    M <- rel_base[sample(nrow(rel_base), n, replace=TRUE), , drop=FALSE]
    g <- c(rep(0,n1), rep(1,n2))
    for (i in which(g==1)) M[i,resp] <- M[i,resp] + shift*mu[resp]
    M[M<0] <- 0; M <- M/rowSums(M)
    a <- adonis2(vegdist(M,"bray") ~ g, permutations=nperm)
    r2[s] <- a$R2[1]
    if (!is.na(a[["Pr(>F)"]][1]) && a[["Pr(>F)"]][1] < 0.05) det <- det+1
  }
  c(R2=mean(r2), power=det/nsim)
}
n_tot <- nrow(rel_g)
scenarios <- list(balanced = c(ceiling(n_tot/2), floor(n_tot/2)))
# add the observed group sizes of the first grouping variable, if available
if (length(GROUP_VARS) && GROUP_VARS[1] %in% names(metadata)) {
  gs <- as.integer(table(factor(metadata[[GROUP_VARS[1]]])))
  if (length(gs) == 2) scenarios$observed <- sort(gs, decreasing=TRUE)
}
power_tab <- data.frame(); set.seed(42); rel_mat <- as.matrix(rel_g)
for (nm in names(scenarios)) {
  n1 <- scenarios[[nm]][1]; n2 <- scenarios[[nm]][2]
  for (sh in c(0.5,1,2,3,5)) {
    r <- sim_power(rel_mat, n1, n2, sh)
    power_tab <- rbind(power_tab, data.frame(scenario=nm, groups=paste0(n1," vs ",n2),
                       effect_R2=round(r["R2"],3), power_pct=round(r["power"]*100,0)))
  }
}
write.csv(power_tab, "results/permanova_power.csv", row.names=FALSE)
save_fig(function() {
  par(mar=c(4.5,4.5,3,1)); esc <- unique(power_tab$scenario); cols <- cb_palette(length(esc))
  plot(NA, xlim=range(power_tab$effect_R2), ylim=c(0,100),
       xlab=expression("Effect size (PERMANOVA R"^2*")"), ylab="Statistical power (%)",
       main=sprintf("PERMANOVA power (n=%d, simulation)", n_tot))
  abline(h=80, lty=2, col="gray50"); text(max(power_tab$effect_R2), 82, "80% power", pos=2, cex=0.8, col="gray40")
  for (i in seq_along(esc)) { s <- power_tab[power_tab$scenario==esc[i],]; s <- s[order(s$effect_R2),]
    lines(s$effect_R2, s$power_pct, type="b", col=cols[i], pch=19, lwd=2) }
  legend("topleft", legend=esc, col=cols, pch=19, lwd=2, bty="n", cex=0.85)
}, "permanova_power", w=8, h=6)

# ── Done ──────────────────────────────────────────────────────────────────────
message("== Analysis complete ==")
message("  figures/  results/  data/  written to the working directory.")
message(sprintf("  vegan %s | decontam %s | R %s",
        packageVersion("vegan"), packageVersion("decontam"), getRversion()))
