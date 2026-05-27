## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(statgen)
set_verbosity("quiet")

## -----------------------------------------------------------------------------
reference_template <- file.path(
  dirname(system.file("extdata", "reference_chr1.bim", package = "statgen")),
  "reference_chr@.bim"
)
reference <- load_reference(reference_template)

sumstats_path <- system.file("extdata", "traits_complete.tsv.gz", package = "statgen")
sumstats <- load_sumstats(sumstats_path, reference)

annotation_paths <- system.file(
  "extdata",
  c("anno1.bed", "anno2.bed"),
  package = "statgen"
)
annotations <- load_annotations(annotation_paths, reference)

num_snp(reference)
head(logpvec(sumstats))
colnames(annomat(annotations))

## -----------------------------------------------------------------------------
sumstats_cache <- tempfile(fileext = ".rds")
save_sumstats_cache(sumstats, sumstats_cache)
sumstats_cached <- load_sumstats_cache(sumstats_cache)
identical(is_present(sumstats_cached), is_present(sumstats))

## -----------------------------------------------------------------------------
ld_root <- system.file("extdata", "ld", package = "statgen")
ld <- load_ld(ld_root)

num_snp(ld)
a1freq(ld)

## -----------------------------------------------------------------------------
scores <- seq_len(num_snp(ld))
multiply_r2(ld, scores)

logp <- rev(seq_len(num_snp(ld)))
fast_prune(logp, ld, r2_threshold = 0.35)

## -----------------------------------------------------------------------------
ref_chr1 <- load_reference(system.file("extdata", "reference_chr1.bim", package = "statgen"))
bed_chr1 <- system.file("extdata", "genotype_1.bed", package = "statgen")
genotype <- load_genotype(sub("\\.bed$", "", bed_chr1), ref_chr1)

geno <- fetch_genotypes_int8(genotype, c(1, 3))
dim(geno)
head(geno)

