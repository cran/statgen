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
binary_annotations <- load_annotations(annotation_paths, reference)
grouped_annotations <- load_annotation(
  system.file("extdata", "grouped.annot", package = "statgen"),
  reference,
  has_header = TRUE,
  group_column = "group"
)

continuous_path <- tempfile(fileext = ".annot")
writeLines(c(
  "chrom\tstart0\tend0\tscore",
  "1\t99\t200\t0.5",
  "1\t299\t400\t1.5",
  "X\t99\t301\t2.5"
), continuous_path)
continuous_annotations <- load_annotation(
  continuous_path,
  reference,
  has_header = TRUE,
  value_columns = "score"
)
annotations <- union_annotations(
  union_annotations(binary_annotations, grouped_annotations),
  continuous_annotations
)

num_snp(reference)
head(logpvec(sumstats))
colnames(annomat(annotations))
is_binary(annotations)
annotation_metadata(annotations)[!is_binary(annotations)]

## -----------------------------------------------------------------------------
sumstats_cache <- tempfile(fileext = ".rds")
save_sumstats_cache(sumstats, sumstats_cache)
sumstats_cached <- load_sumstats_cache(sumstats_cache)
identical(is_present(sumstats_cached), is_present(sumstats))

## -----------------------------------------------------------------------------
swap_root <- system.file("extdata", "allele_swaps", package = "statgen", mustWork = TRUE)
swap_reference <- load_reference(file.path(swap_root, "reference.bim"), shards = "X")
swap_stats <- load_sumstats(file.path(swap_root, "sumstats.tsv.gz"), swap_reference,
                            allow_indel_swaps = TRUE)
swap_genotype <- load_genotype(file.path(swap_root, "genotype", "@"), swap_reference,
                               allow_indel_swaps = TRUE)
eaf_vec(swap_stats)
is_allele_swapped(swap_genotype)
fetch_genotypes_int8(swap_genotype, c(1, 3, 1))

## ----eval=FALSE---------------------------------------------------------------
# trait_by_id <- load_sumstats(path, reference, match_by = "id")
# genotype_by_id <- load_genotype(bfile_prefix, reference, match_by = "id")

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

