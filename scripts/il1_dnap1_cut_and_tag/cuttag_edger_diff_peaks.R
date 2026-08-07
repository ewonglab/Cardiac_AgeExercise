# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Fit edgeR differential-peak models for IL-1/dnAP-1 CUT&Tag data.
# Steps: 1) Load peak counts and metadata; 2) fit contrasts; 3) write differential peaks and downstream annotations.
# =======================================
# Source reformatted from IL1_dnAP1_CUT&Tag/cuttag_edgeR_diff_peaks.R.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

#!/usr/bin/env Rscript

# CUT&Tag/CUT&RUN differential binding with edgeR (peak count matrix)
# Generalised from the original 10-week exercise script.
#
# Required args:
#   --run_dir        Root nf-core run directory (contains results_mm10/)
#   --out_dir        Output directory
#   --mark           Antibody/mark name exactly as in BAM filenames (e.g. C-FOS, H3K9Me3)
#   --cons_prefix    Prefix before condition name in consensus peak filenames (usually "")
#   --conditions     Comma-separated condition levels IN ORDER (e.g. Blank,IL_1_15m,IL_1_30m,IL_1_60m)
#   --contrasts      Comma-separated contrasts as 'A-B' (e.g. IL_1_15m-Blank,IL_1_30m-Blank)
#
# Optional args:
#   --covariate      Blocking factor column name extracted from sample names (e.g. sex)
#                    When provided the design is ~0 + condition + covariate.
#   --sample_filter  Substring that must appear in BAM sample ID to keep it (e.g. Male, Female).
#                    Also restricts which consensus peak files are loaded.
#   --support_min    Min samples a consensus peak must appear in (default 2)
#   --fdr            FDR threshold (default 0.05)
#   --lfc            |log2FC| threshold (default log2(1.5))

suppressPackageStartupMessages({
  library(edgeR)
  library(Rsamtools)
  library(GenomicAlignments)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
})

# ---------- argument parsing ----------

parse_args <- function(x) {
  out <- list()
  for (a in x) {
    if (!grepl('^--', a)) next
    a2 <- sub('^--', '', a)
    if (!grepl('=', a2, fixed=TRUE)) { out[[a2]] <- TRUE; next }
    kv <- strsplit(a2, '=', fixed=TRUE)[[1]]
    out[[kv[[1]]]] <- paste(kv[-1], collapse='=')
  }
  out
}

args <- parse_args(commandArgs(trailingOnly=TRUE))

need <- function(k) {
  if (is.null(args[[k]]) || identical(args[[k]], TRUE) || nchar(args[[k]]) == 0)
    stop('Missing required arg --', k, call.=FALSE)
  args[[k]]
}
opt <- function(k, default=NULL) if (!is.null(args[[k]]) && !identical(args[[k]], TRUE)) args[[k]] else default

run_dir       <- need('run_dir')
out_dir       <- need('out_dir')
mark          <- need('mark')
cons_prefix   <- opt('cons_prefix', '')
conditions    <- strsplit(need('conditions'), ',', fixed=TRUE)[[1]]
contrast_strs <- strsplit(need('contrasts'),  ',', fixed=TRUE)[[1]]

covariate       <- opt('covariate')
sample_filter   <- opt('sample_filter')
support_min     <- as.integer(opt('support_min', '2'))
pool_replicates <- isTRUE(opt('pool_replicates') == 'TRUE' || identical(opt('pool_replicates'), TRUE))
fdr_thr         <- as.numeric(opt('fdr', '0.05'))
lfc_thr         <- as.numeric(opt('lfc', as.character(log2(1.5))))
normalization   <- tolower(trimws(opt('normalization', 'tmm')))
if (!normalization %in% c('tmm', 'spikein'))
  stop('--normalization must be "tmm" or "spikein", got: ', normalization, call.=FALSE)
spikein_log_dir <- opt('spikein_log_dir')

if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

message('run_dir:      ', run_dir)
message('out_dir:      ', out_dir)
message('mark:         ', mark)
message('cons_prefix:  ', cons_prefix)
message('conditions:   ', paste(conditions, collapse=', '))
message('contrasts:    ', paste(contrast_strs, collapse=', '))
message('covariate:    ', if (is.null(covariate)) 'none' else covariate)
message('sample_filter:', if (is.null(sample_filter)) 'none' else sample_filter)
message('support_min:  ', support_min)
message('fdr_thr:      ', fdr_thr)
message('lfc_thr:      ', lfc_thr)
message('normalization:', normalization)

cons_dir <- file.path(run_dir, 'results_mm10', '03_peak_calling', '05_consensus_peaks')
aln_dir  <- file.path(run_dir, 'results_mm10', '02_alignment', 'bowtie2', 'target', 'markdup')
if (is.null(spikein_log_dir))
  spikein_log_dir <- file.path(run_dir, 'results_mm10', '02_alignment', 'bowtie2', 'spikein', 'log')

if (!dir.exists(cons_dir)) stop('Missing consensus dir: ', cons_dir)
if (!dir.exists(aln_dir))  stop('Missing alignment dir: ', aln_dir)

# ---------- helper functions ----------

escape_regex <- function(x) gsub('([][{}()+*^$\\\\.|?])', '\\\\\\1', x)

read_consensus_gr <- function(path, support_min) {
  con <- file(path, open='rt')
  on.exit(close(con), add=TRUE)
  chr <- character(); start0 <- integer(); end <- integer(); sup <- integer()
  n <- 0L; keep <- 0L
  while (length(line <- readLines(con, n=50000L, warn=FALSE)) > 0) {
    n <- n + length(line)
    parts <- strsplit(line, '\t', fixed=TRUE)
    for (p in parts) {
      if (length(p) < 10) next
      s <- suppressWarnings(as.integer(p[[10]]))
      if (is.na(s) || s < support_min) next
      chr    <- c(chr,    p[[1]])
      start0 <- c(start0, suppressWarnings(as.integer(p[[2]])))
      end    <- c(end,    suppressWarnings(as.integer(p[[3]])))
      sup    <- c(sup,    s)
      keep   <- keep + 1L
    }
  }
  message('read_consensus_gr kept ', keep, ' / streamed rows ', n, ' from ', basename(path))
  gr <- GRanges(seqnames=chr, ranges=IRanges(start=start0 + 1L, end=end), support=sup)
  gr <- keepStandardChromosomes(gr, pruning.mode='coarse')
  sort(gr)
}

# Build master peak set from consensus peak files.
# Handles three naming schemes produced by nf-core/cutandrun:
#   (a) {prefix}{cond}_{mark}.seacr.consensus...       (Arm 1: no sex in group name)
#   (b) {prefix}{cond}_{mark}_{sex}.seacr.consensus... (Arms 2/3: sex in group name)
# When sample_filter is set (e.g. "Male"), only files whose name contains that token are kept.
master_peaks_from_conditions <- function(prefix, conditions, mark, support_min,
                                         sample_filter=NULL) {
  re_mark <- escape_regex(mark)
  all_cons <- list.files(cons_dir,
    pattern='\\.seacr\\.consensus\\.peaks\\.awk\\.bed$', full.names=TRUE)

  grs <- list()
  for (cond in conditions) {
    # Pattern: file starts with {prefix}{cond}_ and contains _{mark} as a whole token
    # (i.e. followed by _, ., or end-of-name)
    cond_pat <- paste0('^', escape_regex(paste0(prefix, cond, '_')))
    mark_pat  <- paste0('_', re_mark, '(_|\\.|$)')
    matching  <- all_cons[
      grepl(cond_pat, basename(all_cons), perl=TRUE) &
      grepl(mark_pat, basename(all_cons), perl=TRUE)
    ]
    # Apply sample_filter if specified (e.g. keep only "Male" files)
    if (!is.null(sample_filter) && nzchar(sample_filter)) {
      sf_pat   <- paste0('_', escape_regex(sample_filter), '(_|\\.|$)')
      matching <- matching[grepl(sf_pat, basename(matching), perl=TRUE)]
    }
    if (length(matching) == 0) {
      avail <- paste(basename(all_cons), collapse='\n  ')
      stop('No consensus peak files found for condition "', cond, '" mark "', mark,
           if (!is.null(sample_filter)) paste0(' filter "', sample_filter, '"') else '',
           '\nSearched: ', cons_dir,
           '\nAvailable files:\n  ', avail)
    }
    message('Condition "', cond, '": using ', length(matching),
            ' consensus file(s):\n  ', paste(basename(matching), collapse='\n  '))
    grs <- c(grs, lapply(matching, read_consensus_gr, support_min=support_min))
  }
  red <- reduce(do.call(c, grs), ignore.strand=TRUE)
  red <- keepStandardChromosomes(red, pruning.mode='coarse')
  sort(red)
}

# Build master peak set from individual SEACR replicate BED files.
# Used for pooled (sex-ignored) analyses: counts support across ALL replicates
# in each condition regardless of sex, so n=1 females contribute equally.
# A peak passes if >= support_min replicates in ANY condition support it.
master_peaks_from_replicates <- function(conditions, mark, support_min,
                                         sample_filter=NULL) {
  re_mark  <- escape_regex(mark)
  seacr_dir <- file.path(run_dir, 'results_mm10', '03_peak_calling',
                          '04_called_peaks', 'seacr')
  if (!dir.exists(seacr_dir))
    stop('Missing SEACR dir: ', seacr_dir)

  all_beds <- list.files(seacr_dir,
    pattern='\\.seacr\\.peaks\\.stringent\\.bed$', full.names=TRUE)

  # Keep only files matching the mark
  mark_pat <- paste0('_', re_mark, '(_|\\.)')
  all_beds <- all_beds[grepl(mark_pat, basename(all_beds), perl=TRUE)]

  # Keep only files whose condition prefix matches one of the supplied conditions
  re_conds <- paste(sapply(conditions, escape_regex), collapse='|')
  all_beds <- all_beds[grepl(paste0('^(', re_conds, ')_'),
                              basename(all_beds), perl=TRUE)]

  # Optional sex filter (males-only or females-only sub-analysis)
  if (!is.null(sample_filter) && nzchar(sample_filter))
    all_beds <- all_beds[grepl(paste0('_', escape_regex(sample_filter), '_'),
                                basename(all_beds), perl=TRUE)]

  if (length(all_beds) == 0)
    stop('No SEACR replicate BEDs found for conditions "',
         paste(conditions, collapse=', '), '" mark "', mark, '"',
         '\nSearched: ', seacr_dir)

  message('Building master peaks from ', length(all_beds),
          ' replicate BED file(s) (pool_replicates mode):')
  message('  ', paste(basename(all_beds), collapse='\n  '))

  # Load each replicate as a GRanges, grouped by condition
  load_bed_gr <- function(fp) {
    df <- tryCatch(
      read.table(fp, sep='\t', header=FALSE, stringsAsFactors=FALSE),
      error=function(e) NULL)
    if (is.null(df) || nrow(df) < 1) return(GRanges())
    GRanges(seqnames=df[[1]], ranges=IRanges(start=as.integer(df[[2]])+1L,
                                              end=as.integer(df[[3]])))
  }

  grs_all <- lapply(all_beds, load_bed_gr)
  grs_all <- Filter(function(g) length(g) > 0, grs_all)

  if (length(grs_all) == 0)
    stop('All replicate BED files were empty.')

  # Union all peaks across replicates
  universe <- reduce(do.call(c, grs_all), ignore.strand=TRUE)

  # Count per-universe-peak how many replicates overlap it
  support_mat <- vapply(grs_all, function(gr)
    countOverlaps(universe, gr, ignore.strand=TRUE) > 0L,
    logical(length(universe)))

  # Support = number of replicates supporting each peak (summed within condition)
  # We need the max support within any condition (not global sum, to avoid
  # a peak being rescued only because it appears in many replicates of one cond)
  cond_labels <- sub(paste0('_', re_mark, '.*$'), '', basename(all_beds), perl=TRUE)

  max_support <- integer(length(universe))
  for (cond in conditions) {
    idx <- cond_labels == cond
    if (!any(idx)) next
    if (sum(idx) == 1) {
      cond_sup <- as.integer(support_mat[, idx])
    } else {
      cond_sup <- rowSums(support_mat[, idx, drop=FALSE])
    }
    max_support <- pmax(max_support, cond_sup)
  }

  keep <- max_support >= support_min
  message('support_min=', support_min, ' (out of max ',
          max(table(cond_labels)), ' replicates per condition)')
  message('Peaks passing: ', sum(keep), ' / ', length(universe),
          ' (', round(100*mean(keep), 1), '%)')

  out <- universe[keep]
  out <- keepStandardChromosomes(out, pruning.mode='coarse')
  sort(out)
}

# Find BAM files for the given mark.
# nf-core names BAMs as {group}_R{rep}.target.markdup.sorted.bam where {group}
# is the samplesheet group column.  The mark may appear anywhere within the
# group name, separated by underscores (e.g. Blank_C-FOS_R1 or
# AAV9_Control_C-FOS_Male_R1 or 1YS_H3K9Me3_Male_R1).
# The pattern uses _{mark}_ token matching to avoid partial matches
# (e.g. _C-JUN_ must not match _P-C-JUN_).
find_bams <- function(mark, sample_filter=NULL) {
  re_mark <- escape_regex(mark)
  # Match _{mark}_ or _{mark}_R (handles both Arm1 and Arm2/3 naming)
  pat <- paste0('_', re_mark, '(_.*)?_R[0-9]+\\.target\\.markdup\\.sorted\\.bam$')
  b   <- sort(list.files(aln_dir, pattern=pat, full.names=TRUE))
  if (!is.null(sample_filter) && nzchar(sample_filter))
    b <- b[grepl(paste0('_', escape_regex(sample_filter), '_'), basename(b), perl=TRUE)]
  if (length(b) == 0)
    stop('No BAMs found for mark "', mark, '"',
         if (!is.null(sample_filter)) paste0(' filter "', sample_filter, '"') else '',
         ' in ', aln_dir)
  bai <- paste0(b, '.bai')
  if (any(!file.exists(bai)))
    stop('Missing BAM index for: ', paste(b[!file.exists(bai)], collapse=', '))
  b
}

# Build a sample metadata table from BAM paths.
# BAM basename format: {group}_R{rep}.target.markdup.sorted.bam
# where {group} is the samplesheet group column, e.g.:
#   Arm 1: Blank_C-FOS_R1          → group=Blank_C-FOS,       cond=Blank
#   Arm 2: AAV9_Control_C-FOS_Male_R1 → group=AAV9_Control_C-FOS_Male, cond=AAV9_Control, sex=Male
#   Arm 3: 1YS_H3K9Me3_Male_R1     → group=1YS_H3K9Me3_Male,  cond=1YS, sex=Male
sample_table_from_bams <- function(bams, mark, conditions, covariate=NULL) {
  bn  <- basename(bams)
  sid <- sub('\\.target\\.markdup\\.sorted\\.bam$', '', bn)

  # Replicate number: last _R{digits} in the sample ID
  rep_str <- sub('^.*_R([0-9]+)$', '\\1', sid, perl=TRUE)
  rep     <- suppressWarnings(as.integer(rep_str))
  if (any(is.na(rep)))
    warning('Could not parse replicate number from: ',
            paste(sid[is.na(rep)], collapse=', '))

  # Group = sample ID minus the trailing _R{rep}
  group <- sub('_R[0-9]+$', '', sid, perl=TRUE)

  # Condition: strip _{mark}[_{anything}] from the right of the group name.
  # nf-core group names are always {condition}_{mark}[_{sex}], so stripping
  # from _{mark} onward gives the condition prefix directly.
  re_mark_strip <- escape_regex(mark)
  cond <- vapply(group, function(g) {
    candidate <- sub(paste0('_', re_mark_strip, '(_.+)?$'), '', g, perl=TRUE)
    matched <- conditions[conditions == candidate]
    if (length(matched) == 0) NA_character_ else matched[[1]]
  }, character(1))

  # Drop BAMs from other arms (condition prefix doesn't match any supplied condition)
  unmatched <- is.na(cond)
  if (any(unmatched)) {
    message('Dropping ', sum(unmatched), ' BAM(s) from other arms (no condition match): ',
            paste(sid[unmatched], collapse=', '))
    bams  <- bams[!unmatched]
    sid   <- sid[!unmatched]
    rep   <- rep[!unmatched]
    group <- group[!unmatched]
    cond  <- cond[!unmatched]
  }

  if (length(cond) == 0)
    stop('No BAMs matched any supplied condition. Expected: ',
         paste(conditions, collapse=', '))

  missing_conds <- setdiff(conditions, unique(cond))
  if (length(missing_conds) > 0)
    stop('No BAMs found for condition(s): ', paste(missing_conds, collapse=', '),
         '\nExpected: ', paste(conditions, collapse=', '))

  df <- data.frame(sample_id=sid, group=group, cond=cond, rep=rep, bam=bams,
                   stringsAsFactors=FALSE)

  # Extract covariate (e.g. sex) as whatever remains after stripping the
  # matched condition AND the mark from the group name.
  if (!is.null(covariate)) {
    re_mark_esc <- escape_regex(mark)
    df[[covariate]] <- vapply(seq_len(nrow(df)), function(i) {
      g <- df$group[i]; c <- df$cond[i]
      # Remove condition prefix (and trailing _)
      rem <- sub(paste0('^', escape_regex(c), '_?'), '', g, perl=TRUE)
      # Remove mark (and surrounding _)
      rem <- gsub(paste0('_?', re_mark_esc, '_?'), '_', rem, perl=TRUE)
      # Trim leading/trailing underscores
      rem <- gsub('^_+|_+$', '', rem, perl=TRUE)
      if (nzchar(rem)) rem else NA_character_
    }, character(1))
    n_na <- sum(is.na(df[[covariate]]))
    if (n_na > 0)
      warning(n_na, ' sample(s) could not have covariate "', covariate,
              '" extracted: ', paste(df$sample_id[is.na(df[[covariate]])], collapse=', '))
  }

  df
}

count_peaks <- function(peaks, bams) {
  bfs <- BamFileList(lapply(bams, function(p) BamFile(p, yieldSize=2000000L)))
  se  <- summarizeOverlaps(features=peaks, reads=bfs, mode='Union',
                           singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
  as.matrix(assay(se))
}

# Parse a contrast string like 'IL_1_15m-Blank' into c(numerator, denominator).
# Uses the condition vector to avoid splitting on hyphens inside condition names.
parse_contrast <- function(s, conditions) {
  conds_sorted <- conditions[order(nchar(conditions), decreasing=TRUE)]
  for (a in conds_sorted) {
    # Try: a is the numerator prefix
    m <- regmatches(s, regexec(paste0('^', escape_regex(a), '-(.+)$'), s, perl=TRUE))[[1]]
    if (length(m) == 3) return(c(a, m[[3]]))   # numerator=a, denominator=m[3]
    # Try: a is the denominator suffix
    m2 <- regmatches(s, regexec(paste0('^(.+)-', escape_regex(a), '$'), s, perl=TRUE))[[1]]
    if (length(m2) == 3) return(c(m2[[2]], a)) # numerator=m2[2], denominator=a
  }
  # Fallback: split on last hyphen
  idx <- max(gregexpr('-', s, fixed=TRUE)[[1]])
  if (idx < 1) stop('Cannot parse contrast: ', s)
  c(substr(s, 1, idx-1), substr(s, idx+1, nchar(s)))
}

# Parse spike-in aligned fragment counts from nf-core bowtie2 spike-in log files.
# Each log is at {spikein_log_dir}/{sample_id}.spikein.bowtie2.log
# Spike-in fragments = concordantly aligned exactly once + concordantly aligned >1 times.
# Returns a named numeric vector: sample_id -> spike_in_fragment_count.
get_spikein_counts <- function(sample_ids, spikein_log_dir) {
  counts <- sapply(sample_ids, function(sid) {
    path <- file.path(spikein_log_dir, paste0(sid, '.spikein.bowtie2.log'))
    if (!file.exists(path)) {
      warning('Spike-in log not found for ', sid, ': ', path)
      return(NA_real_)
    }
    txt <- readLines(path)
    parse_n <- function(pattern) {
      m <- grep(pattern, txt, value=TRUE)
      if (length(m) == 0) return(0)
      as.numeric(sub('^\\s*(\\d+).*', '\\1', m[1]))
    }
    n_unique <- parse_n('aligned concordantly exactly 1 time')
    n_multi  <- parse_n('aligned concordantly >1 times')
    total    <- n_unique + n_multi
    message('  ', sid, ': spike-in fragments = ', total,
            ' (unique=', n_unique, ', multi=', n_multi, ')')
    total
  })
  counts
}

run_edger <- function(counts, sample_df, conditions, contrast_strs,
                      covariate=NULL, normalization='tmm', spikein_log_dir=NULL) {
  cond <- factor(sample_df$cond, levels=conditions)

  # Check covariate is usable (>= 2 distinct non-NA levels)
  if (!is.null(covariate) && covariate %in% colnames(sample_df)) {
    cov_vals <- sample_df[[covariate]]
    n_cov    <- length(unique(na.omit(cov_vals)))
    if (n_cov < 2) {
      message('Covariate "', covariate, '" has < 2 distinct levels after filtering -- ignoring.')
      covariate <- NULL
    } else if (any(is.na(cov_vals))) {
      message('Some samples lack covariate "', covariate, '" -- setting covariate to NULL.')
      covariate <- NULL
    }
  } else {
    covariate <- NULL
  }

  if (!is.null(covariate)) {
    cov_fac <- factor(sample_df[[covariate]])
    design  <- model.matrix(~0 + cond + cov_fac)
    colnames(design) <- c(levels(cond),
                          paste0(covariate, '_', levels(cov_fac)[-1]))
  } else {
    design <- model.matrix(~0 + cond)
    colnames(design) <- levels(cond)
  }

  y    <- DGEList(counts=counts)
  keep <- filterByExpr(y, design=design)
  y    <- y[keep, , keep.lib.sizes=FALSE]

  if (normalization == 'spikein') {
    message('Reading spike-in counts from: ', spikein_log_dir)
    si <- get_spikein_counts(sample_df$sample_id, spikein_log_dir)
    if (any(is.na(si) | si == 0)) {
      warning('Some spike-in counts are missing or zero -- ',
              'falling back to TMM normalization.')
      y <- calcNormFactors(y)
    } else {
      # Spike-in scale factors: sf_i = si_i / lib_size_i, geometric-mean centred
      # so that norm.factors multiply to 1 across samples (edgeR convention).
      # Samples with more spike-in reads relative to their library get a larger
      # norm.factor, meaning edgeR scales them down — correct behaviour.
      sf <- si / y$samples$lib.size
      sf <- sf / exp(mean(log(sf)))
      y$samples$norm.factors <- sf
      message('Spike-in norm factors (geometric-mean centred):')
      for (i in seq_along(sf))
        message('  ', sample_df$sample_id[i], ': si=', si[i],
                '  lib=', y$samples$lib.size[i], '  nf=', round(sf[i], 4))
    }
  } else {
    y <- calcNormFactors(y)
    message('TMM norm factors: ', paste(round(y$samples$norm.factors, 4), collapse=', '))
  }

  n_per_cond <- table(sample_df$cond)
  is_underpowered <- any(n_per_cond < 2)

  if (is_underpowered) {
    message('WARNING: some conditions have < 2 samples.')
    message('Using prior.df=5 for dispersion estimation.')
    message('Results should be treated as DESCRIPTIVE/EXPLORATORY only.')
    y   <- estimateDisp(y, design, prior.df=5)
    fit <- glmQLFit(y, design, prior.count=2)
  } else {
    y   <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
  }

  results <- list(y=y, keep=keep)

  for (cs in contrast_strs) {
    parsed <- parse_contrast(cs, conditions)
    num    <- parsed[1]; den <- parsed[2]
    label  <- paste0(num, '_vs_', den)

    # Build contrast vector manually (handles condition names starting with digits
    # which makeContrasts cannot parse as R expressions).
    v <- numeric(ncol(design))
    names(v) <- colnames(design)
    if (!num %in% names(v))
      stop('Numerator "', num, '" not in design columns: ', paste(names(v), collapse=', '))
    if (!den %in% names(v))
      stop('Denominator "', den, '" not in design columns: ', paste(names(v), collapse=', '))
    v[num] <-  1
    v[den] <- -1
    ct <- matrix(v, ncol=1, dimnames=list(colnames(design), label))

    tt <- topTags(glmQLFTest(fit, contrast=ct[, 1]), n=Inf, sort.by="none")$table
    results[[label]] <- tt
    message('Contrast ', label, ': done')
  }

  results
}

write_qc_plots <- function(y, out_prefix, sample_df) {
  lab <- if ('cond' %in% colnames(sample_df)) sample_df$cond else sample_df$sample_id
  pdf(paste0(out_prefix, '.MDS.pdf'), width=7, height=6)
  plotMDS(y, labels=lab, col=as.integer(factor(lab)))
  dev.off()
  pdf(paste0(out_prefix, '.BCV.pdf'), width=7, height=6)
  plotBCV(y)
  dev.off()
}

annotate_peaks_df <- function(gr) {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  as.data.frame(annotatePeak(gr, TxDb=txdb, annoDb='org.Mm.eg.db', tssRegion=c(-3000, 3000)))
}

sig_filter <- function(tab, direction=c('any','pos','neg')) {
  direction <- match.arg(direction)
  ok <- (tab$FDR <= fdr_thr) & (abs(tab$logFC) >= lfc_thr)
  if (direction == 'pos') ok <- ok & (tab$logFC > 0)
  if (direction == 'neg') ok <- ok & (tab$logFC < 0)
  tab[ok, , drop=FALSE]
}

sig_to_gr <- function(df)
  GRanges(seqnames=df$chr, ranges=IRanges(start=df$start, end=df$end))

write_bed <- function(df, path) {
  if (nrow(df) == 0) return(invisible(NULL))
  bed <- data.frame(chr=df$chr, start0=df$start - 1L, end=df$end)
  write.table(bed, file=path, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# ---------- main ----------

bams      <- find_bams(mark, sample_filter)
sample_df <- sample_table_from_bams(bams, mark, conditions, covariate)
write.csv(sample_df, file=file.path(out_dir, paste0(mark, '.samples.csv')), row.names=FALSE)

message('\nSample table:')
print(sample_df[, intersect(c('sample_id','cond', covariate, 'rep'), colnames(sample_df))])

if (pool_replicates) {
  peaks_master <- master_peaks_from_replicates(conditions, mark, support_min, sample_filter)
} else {
  peaks_master <- master_peaks_from_conditions(cons_prefix, conditions, mark,
                                               support_min, sample_filter)
}
message('master peaks: ', length(peaks_master))

master_bed <- data.frame(chr=as.character(seqnames(peaks_master)),
                         start0=start(peaks_master)-1L, end=end(peaks_master))
write.table(master_bed,
  file=file.path(out_dir, paste0(mark, '.master_peaks.support', support_min, '.bed')),
  sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

message('Counting reads for ', nrow(sample_df), ' BAMs...')
counts <- count_peaks(peaks_master, sample_df$bam)
colnames(counts) <- sample_df$sample_id

                            
                            counts_df <- data.frame(
  peak_id = paste0(as.character(seqnames(peaks_master)), ":", start(peaks_master), "-", end(peaks_master)),
  chr = as.character(seqnames(peaks_master)),
  start = start(peaks_master),
  end = end(peaks_master),
  counts,
  check.names = FALSE
)

write.csv(
  counts_df,
  file = file.path(out_dir, paste0(mark, ".peak_counts.raw.support", support_min, ".csv")),
  row.names = FALSE
)
ed <- run_edger(counts, sample_df, conditions, contrast_strs, covariate, normalization, spikein_log_dir)

write_qc_plots(ed$y, file.path(out_dir, paste0(mark, '.edgeR')), sample_df)

peaks_kept <- peaks_master[ed$keep]
peak_id    <- paste0(as.character(seqnames(peaks_kept)), ':',
                     start(peaks_kept), '-', end(peaks_kept))

attach_coords <- function(tab) {
  tab$peak_id <- peak_id
  tab$chr     <- as.character(seqnames(peaks_kept))
  tab$start   <- start(peaks_kept)
  tab$end     <- end(peaks_kept)
  tab
}

contrast_labels <- setdiff(names(ed), c('y', 'keep'))
for (label in contrast_labels) {
  tab      <- attach_coords(ed[[label]])
  tab_gain <- sig_filter(tab, 'pos')
  tab_loss <- sig_filter(tab, 'neg')

  write.csv(tab,
    file=file.path(out_dir, paste0(mark, '.edgeR.', label, '.full.csv')),
    row.names=FALSE)
  write.csv(tab_gain,
    file=file.path(out_dir, paste0(mark, '.edgeR.', label, '.gain.sig.csv')),
    row.names=FALSE)
  write.csv(tab_loss,
    file=file.path(out_dir, paste0(mark, '.edgeR.', label, '.loss.sig.csv')),
    row.names=FALSE)
  write_bed(tab_gain, file.path(out_dir, paste0(mark, '.', label, '.gain.sig.bed')))
  write_bed(tab_loss, file.path(out_dir, paste0(mark, '.', label, '.loss.sig.bed')))

  if (nrow(tab_gain) > 0) {
    anno  <- annotate_peaks_df(sig_to_gr(tab_gain))
    anno2 <- anno[, setdiff(colnames(anno), c('seqnames','start','end','width')), drop=FALSE]
    write.csv(cbind(tab_gain, anno2),
      file=file.path(out_dir, paste0(mark, '.edgeR.', label, '.gain.sig.annotated.csv')),
      row.names=FALSE)
  }
  if (nrow(tab_loss) > 0) {
    anno  <- annotate_peaks_df(sig_to_gr(tab_loss))
    anno2 <- anno[, setdiff(colnames(anno), c('seqnames','start','end','width')), drop=FALSE]
    write.csv(cbind(tab_loss, anno2),
      file=file.path(out_dir, paste0(mark, '.edgeR.', label, '.loss.sig.annotated.csv')),
      row.names=FALSE)
  }
  message('Contrast ', label, ': ', nrow(tab_gain), ' gain, ',
          nrow(tab_loss), ' loss significant peaks')
}

message('Done: ', mark)
