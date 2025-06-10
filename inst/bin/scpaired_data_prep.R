#########################################################
# processes scATAC and scRNA data for scPaired numbat
#########################################################
suppressMessages({library(dplyr)          
library(data.table)     
library(tibble)        
library(optparse)
})

message('prepare scPaired numbat input')
option_list = list(
  make_option('--alleleCounts_RNA', type = "character", default = NULL,
              help = "RNA-based allele counts"),
  make_option('--alleleCounts_ATAC', type = "character", default = NULL,
              help = "ATAC-based allele counts"),
  make_option('--binCounts_RNA', type = "character", default = NULL,
              help = "RNA-based binned count matrix"),
  make_option('--binCounts_ATAC', type = "character", default = NULL,
              help = "ATAC-based binned count matrix"),
  make_option('--outAlleleCountsFile', type = "character", default = NULL,
              help = "output allele counts file"),
  make_option('--outBinCountsFile', type = "character", default = NULL,
              help = "output binned count matrix file"),
  make_option('--runMultiomeAsPaired', action = "store_true", default = FALSE,
              help = "Running multiome as paired... Such that the same barcode across two modalities have an additional 'atac' or 'rna' tag.")
)

args = parse_args(OptionParser(option_list = option_list))
invisible(list2env(args,environment()))


analyze_lost_bin_counts <- function(binCounts_RNA, binCounts_ATAC) {
  # analyzes lost bin counts when cbinding RNA and ATAC bin counts
  # based on how much of the chr p/q values are lost
  library(GenomicRanges)
  chr_pq_file <- 'data/cytoBand_pq.txt'
  #read in p/q ranges
  cyto_band = read.table(chr_pq_file, header = FALSE)
  cyto_band_gr = GRanges(
    seqnames = cyto_band[,1],
    ranges = IRanges(start = cyto_band[,2], end = cyto_band[,3]),
    pq = cyto_band[,4],
    name = paste0(cyto_band[,1], ':', cyto_band[,2], '-', cyto_band[,3], '_', cyto_band[,4])
  )
  cyto_band_gr$width = width(cyto_band_gr)

  rnavatac = setdiff(rownames(binCounts_RNA), rownames(binCounts_ATAC))  
  atacvrna = setdiff(rownames(binCounts_ATAC), rownames(binCounts_RNA))  
  bin_diff = c(rnavatac, atacvrna)

  rowname_to_gr <- function(rnames){
    split_rnames_chr <- strsplit(rnames, ":|-")[[1]]
    split_chr = list(chr = split_rnames_chr[1], 
                     start = as.numeric(split_rnames_chr[2]), 
                     end = as.numeric(split_rnames_chr[3]),
                     name = rnames)
    return(split_chr)
  }

  get_diff_gr <- function(split_chr_list){
    gr_gtf = GRanges(
        seqnames = split_chr_list$chr,
        ranges = IRanges(start = split_chr_list$start, 
                         end = split_chr_list$end),
        name = split_chr_list$name
        )
    return(gr_gtf)
  }
  split_chr = lapply(bin_diff, rowname_to_gr)
  split_chr_gr = lapply(split_chr, get_diff_gr)
  split_chr_gr = do.call(c, split_chr_gr)
  split_chr_gr$width = width(split_chr_gr)

  ov = findOverlaps(split_chr_gr, cyto_band_gr)
  ov_df = as.data.frame(ov)
  ov_df$pq = cyto_band_gr$pq[ov_df$subjectHits]
  ov_df$bin_diff = bin_diff[ov_df$queryHits]
  ov_df <- data.frame(
		bin_index = queryHits(ov),
		band_index  = subjectHits(ov),
		bin_name  = split_chr_gr$name[queryHits(ov)],
		band_name = cyto_band_gr$name[subjectHits(ov)],
		overlap_bp = split_chr_gr$width[queryHits(ov)]/cyto_band_gr$width[subjectHits(ov)],
		stringsAsFactors = FALSE
  )
  # Lost p/q counts
  prop_lost = ov_df %>% group_by(band_name) %>% summarise(prop_lost = sum(overlap_bp)) %>% arrange(desc(prop_lost))
  message("Lost p/q counts: ", prop_lost)

  return(list(split_chr_gr, cyto_band_gr, ov, ov_df, prop_lost))
}



combine_binned_counts <- function(binCounts_RNA, binCounts_ATAC, output_file, addBarcodeSuff = FALSE) {
  if (!all(file.exists(binCounts_RNA, binCounts_ATAC))) {
    missing_files <- c(binCounts_RNA, binCounts_ATAC)[!file.exists(c(binCounts_RNA, binCounts_ATAC))]
    stop("The following bin count files do not exist: ", 
         paste(missing_files, collapse = ", "))
  }
  # read RNA & ATAC bin counts
  message("Reading RNA bin counts from: ", binCounts_RNA)
  if (grepl("\\.tsv(\\.gz)?$", binCounts_RNA)) {
    rna_counts <- read.table(binCounts_RNA, header = TRUE, sep = "\t", row.names = 1, check.names = F) %>%
      as.matrix()
  } else if (grepl("\\.rds$", binCounts_RNA)) {
    rna_counts <- readRDS(binCounts_RNA)
    if (!is.matrix(rna_counts)) {
      rna_counts <- as.matrix(rna_counts)
    }
  } else {
    stop("Unsupported file format for RNA bin counts. Use .tsv, .tsv.gz, or .rds")
  }
  
  message("Reading ATAC bin counts from: ", binCounts_ATAC)
  if (grepl("\\.tsv(\\.gz)?$", binCounts_ATAC)) {
    #atac_counts <- data.table::fread(binCounts_ATAC) %>%
    atac_counts <- read.table(binCounts_ATAC, header = TRUE, sep = "\t", row.names = 1, check.names = F) %>%
      as.matrix()
  } else if (grepl("\\.rds$", binCounts_ATAC)) {
    atac_counts <- readRDS(binCounts_ATAC)
    if (!is.matrix(atac_counts)) {
      atac_counts <- as.matrix(atac_counts)
    }
  } else {
    stop("Unsupported file format for ATAC bin counts. Use .tsv, .tsv.gz, or .rds")
  }
  
  # row bind the matrices ? not sure why we are doing rowbinding here. it should be column binding?
  message("Combining bin count matrices on intersected bins")
  shared_features <- intersect(rownames(rna_counts),rownames(atac_counts))
  # change column names of rna_counts and atac_counts to have suffix:
  if(addBarcodeSuff == TRUE){
	  colnames(rna_counts) = paste0(colnames(rna_counts), '_RNA')
	  colnames(atac_counts) = paste0(colnames(atac_counts), '_ATAC')
  }
  combined_counts <- cbind(rna_counts[shared_features,], atac_counts[shared_features,])
  # Write output
  message("Writing combined bin counts to: ", output_file)
  if (grepl("\\.rds$", output_file)) {
    saveRDS(combined_counts, file = output_file)
  } else {
    data.table::fwrite(as.data.frame(combined_counts), output_file, 
                       sep = "\t", quote = FALSE, row.names = TRUE)
  }
  
  return(output_file)
}


#' Combine and harmonize allele counts from RNA and ATAC modalities
#' @param alleleCounts_RNA Character. Path to the RNA-based allele counts file
#' @param alleleCounts_ATAC Character. Path to the ATAC-based allele counts file
#' @param output_file Character. File for combined allele outputs 
#' @param compress Logical. Whether to compress the output file (default: TRUE)
#' @return Character. Path to the output file
combine_allele_counts <- function(alleleCounts_RNA, 
                                  alleleCounts_ATAC, output_file, 
				  addBarcodeSuff = FALSE,
                                  compress = TRUE) {
  geno_files <- c(alleleCounts_RNA, alleleCounts_ATAC)
  
  if (!all(file.exists(geno_files))) {
    missing_files <- geno_files[!file.exists(geno_files)]
    stop("The following genotype files do not exist: ", 
         paste(missing_files, collapse = ", "))
  }
  
  
#  message("Harmonizing allele counts for sample: ", sample_id)
  
  # Read first allele counts file (RNA)
  df1 <- tryCatch({
    data.table::fread(alleleCounts_RNA, header = TRUE, sep = "\t") %>% 
      filter(CHROM %in% 1:22) %>%
      mutate(CHROM = factor(CHROM, 1:22))
  }, error = function(e) {
    stop("Error reading RNA allele counts file: ", e$message)
  })
  
  # Extract unique SNPs and genotypes from RNA allele counts file
  geno1 <- df1 %>% distinct(snp_id, GT)
  
  # Read second allele counts file (ATAC)
  df2 <- tryCatch({
    data.table::fread(alleleCounts_ATAC, header = TRUE, sep = "\t") %>% 
      filter(CHROM %in% 1:22) %>%
      mutate(CHROM = factor(CHROM, 1:22))
  }, error = function(e) {
    stop("Error reading ATAC allele counts file: ", e$message)
  })
  
  # Extract unique SNPs and genotypes from ATAC allele counts file
  geno2 <- df2 %>% distinct(snp_id, GT)
  
  # Find shared SNPs between the two datasets
  sharedSNP <- intersect(geno1$snp_id, geno2$snp_id)
  message("Found ", length(sharedSNP), " shared SNPs between modalities")
  
  # Create reference dataframe from second file's genotypes
  snp_df <- geno2 %>% 
    filter(snp_id %in% sharedSNP) %>% 
    column_to_rownames("snp_id")
  
  # Update genotypes in first file to match second file for shared SNPs
  geno1 <- geno1 %>%
    mutate(GT = case_when(
      snp_id %in% sharedSNP ~ snp_df[snp_id, "GT"],
      TRUE ~ GT
    ))
  
  geno1_df <- geno1 %>% column_to_rownames("snp_id")
  df1$GT <- geno1_df[df1$snp_id, "GT"]

  if(addBarcodeSuff == TRUE){
	  df1$cell = paste0(df1$cell, '_RNA')
	  df2$cell = paste0(df2$cell, '_ATAC')
  }
  combined_df <- rbind(df1, df2)
  
#  output_file <- file.path(output_dir, paste0(sample_id, "_allele_counts_consistent.tsv"))
  
  tryCatch({
    data.table::fwrite(
      combined_df, 
      output_file, 
      quote = FALSE, 
      row.names = FALSE, 
      #nThread = 3, 
      sep = "\t"
    )
    message("Successfully wrote combined allele counts to: ", output_file)
  }, error = function(e) {
    stop("Error writing combined file: ", e$message)
  })
  
  if (compress) {
    tryCatch({
      system2("gzip", args = output_file)
      output_file <- paste0(output_file, ".gz")
      message("Successfully compressed output file")
    }, error = function(e) {
      warning("Failed to compress output file: ", e$message)
    })
  }
  
  return(output_file)
}

#########################################################
# executes function on args
#########################################################
combine_allele_counts(alleleCounts_RNA, 
		      alleleCounts_ATAC, outAlleleCountsFile, 
		      addBarcodeSuff = runMultiomeAsPaired,
		      compress = TRUE)
combine_binned_counts(binCounts_RNA,
		      binCounts_ATAC,
		      addBarcodeSuff = runMultiomeAsPaired,
		      outBinCountsFile)
