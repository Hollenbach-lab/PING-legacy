# Copyright 2016 Wesley Marin, Jill Hollenbach, Paul Norman
#
# This file is part of PING.
#
# PING is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PING is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PING.  If not, see <http://www.gnu.org/licenses/>.

ping_allele_caller <- function(
  sample.location = "PING_sequences/",
  fastq.pattern.1 = "_1.fastq.gz",
  fastq.pattern.2 = "_2.fastq.gz",
  bowtie.threads = 4,
  supported.loci = c("2DL1", "2DL23", "2DL4", "2DL5", "2DS3", "2DS4", "2DS5", "2DP1", "3DL1", "3DS1", "3DL2", "3DL3"),
  ping.gc.output = "Combined_results.csv",
  results.directory = ""
){
  
  library(data.table) ## Used for fread
  library(ape) ## Used for read.dna an seg.sites
  library(stringr)
  source("Resources/locus_functions.R", local = TRUE)
  
  ###################
  # Helper functions --------------------------------------------------------
  ###################
  
  
  # Creates results directory
  results_directory <- function() {
    cat("----- Getting PING ready -----\n\n")
    
    if(results.directory != ""){
      save_to <- results.directory
    }else{
      save_to <- paste0("Caller_results_", format(Sys.time(), "%Y_%m_%d_%H:%M"), "/")
      
      count <- 1
      while(file.exists(save_to)) {
        save_to <- paste0("Caller_results_", format(Sys.time(), "%Y_%m_%d_%H:%M"), "_", count, "/")
        count <- count + 1
      }
    }
    dir.create(save_to)
    
    cat(paste("Results being saved to", save_to, "\n\n"))
    
    return(save_to)
  }
  
  # Creates results folders
  ping.ready <- function() {
    dir.create(results.directory, showWarnings = F)
    dir.create(paste0(results.directory, "Vcf"), showWarnings = F)
    dir.create(paste0(results.directory, "KIRcaller"), showWarnings = F)
    dir.create(paste0(results.directory, "Fastq"), showWarnings = F)
    
    cat("Results directories created.\n\n")
  }
  
  # Finds sequences
  get_sequence_list <- function() {
    
    sequence_list = list.files(file.path(sample.location), pattern = fastq.pattern.1)
    
    if (is.na(sequence_list[1])) {
      string <- paste("No sequences found in", sample.location, "using fastq pattern", fastq.pattern.1)
      stop(string)
    } else {
      sequence_list <- gsub(fastq.pattern.1, "", sequence_list)
      cat(paste("Found sequences: ", paste(sequence_list, collapse = "\n"), sep = "\n"))
      cat("\n")
      
      return(sequence_list)
    }
    
  }
  
  # Reads in PING_GC results
  gc.results <- function(){
    gc.table <- read.csv(ping.gc.output, check.names = F)
    rownames(gc.table) <- gc.table[,1]
    gc.table <- gc.table[,-1, drop=F]
    return(gc.table)
  }
  
  # Caps sequences at 120,000 lines
  cut_fastq <- function(file.name, read.cap, post.file.name, is.gz) {
    
    # Uncompressing if is.gz and moving to post.file.name, this happens one at a time, so space should not be a concern
    if(is.gz){
      file_contents <- fread(paste("zcat", file.name), sep="\n", nrows = read.cap, header = F)
    }else{
      file_contents <- fread(file.name, sep="\n", nrows = read.cap, header = F)
    }
    
    write.table(file_contents, file = post.file.name, quote = F, row.names = F, col.names = F)
  }
  
  # Finds files that are smaller than the recommended size
  files_too_small <- function(sequence.list, is.gz){
    small_samples <- NA
    
    for(sample in sequence.list){
      
      # Get the line counts, both for gzip and regular files
      if(is.gz){
        line_count <- as.numeric(system2("zcat", c(paste0(sample.location, sample, fastq.pattern.1), "|", "wc", "-l"), stdout = T))
      }else{
        line_count <- as.numeric(system2("wc", c("-l", "<", paste0(sample.location, sample, fastq.pattern.1)), stdout = T))
      }
      
      if(line_count < 280000){
        small_samples <- c(small_samples, paste0(sample.location, sample, fastq.pattern.1))
      }
      
      # Same thing, but for _2 files
      if(is.gz){
        line_count <- as.numeric(system2("zcat", c(paste0(sample.location, sample, fastq.pattern.2), "|", "wc", "-l"), stdout = T))
      }else{
        line_count <- as.numeric(system2("wc", c("-l", "<", paste0(sample.location, sample, fastq.pattern.2)), stdout = T))
      }
      
      if(line_count < 280000){
        small_samples <- c(small_samples, paste0(sample.location, sample, fastq.pattern.2))
      }
    }
    
    small_samples <- small_samples[!is.na(small_samples)]
    
    if(length(small_samples) > 0){
      error.string <- "\n\nThese files have fewer than the recommended lines, results might be effected:\n"
      cat(error.string)
      error.log(error.string)
      
      cat(paste(small_samples), sep='\n')
      error.log(paste(small_samples, collapse='\n'))
      
      cat("\nContinuing PING analysis...\n\n")
    }
    
    return(small_samples)
  }
  
  # Create consensus sequence -- NOT WORKING
  # create_consensus <- function(sample, vcf.file, filter.directory, results.directory){
  #   
  #   
  # }
  
  # Permutates allele codes
  permutator <- function(codes) {
    sections <- length(codes)
    section_length <- length(codes) - 1
    total_length <- section_length * sections
    final_codes <- data.frame(matrix(0, nrow = total_length, ncol = 2))
    
    for(i in 0:(sections-1)){
      other_codes <- codes[!codes == codes[i + 1]]
      for(j in 1:section_length){
        final_codes[i*section_length + j, 1] <- codes[i + 1]
        final_codes[i*section_length + j, 2] <- other_codes[j]
      }
    }
    
    return(final_codes)
  }
  
  # Build primer table from forward and reverse primer files
  make_primer_table <- function(primerlist.file) {
    primer_table = read.delim(primerlist.file, header = FALSE)
    colnames(primer_table) <- c("Locus", "Primer")
    
    primer_table$Primer <- gsub(" ", "", primer_table[["Primer"]])
    
    return(primer_table)
  }
  
  # Counting primer matches to a table
  count_primer_matches <- function(primer.table, sample, is.gz) {
    primer_count_table <- data.frame(matrix(0, nrow=length(primer.table[,1]), ncol=length(sample)))
    rownames(primer_count_table) <- primer.table$Locus
    colnames(primer_count_table) <- sample
    
    counter <- 1
    while(counter <= 2) {
      
      if(counter == 1) {
        
        if(is.gz){
          cat(paste0(sample, no.gz.pattern.1, "\n"))
          inFile <- paste0(sample, no.gz.pattern.1)
        }else{
          cat(paste0(sample, fastq.pattern.1, "\n"))
          inFile <- paste0(sample, fastq.pattern.1)
        }
        
      }else{
        
        if(is.gz){
          cat(paste0(sample, no.gz.pattern.2, "\n"))
          inFile <- paste0(sample, no.gz.pattern.2)
        }else{
          cat(paste0(sample, fastq.pattern.2, "\n"))
          inFile <- paste0(sample, fastq.pattern.2)
        }
        
      }
      
      file_contents <- fread(inFile, sep="\n", header=FALSE)
      file_contents <- file_contents[seq(2, length(file_contents[[1]]), 4)]
      
      for(j in 1:length(primer.table[,"Primer"])) {
        primer_count_table[j, sample] <- sum(length(grep(primer.table[j, "Primer"], file_contents[[1]], fixed=TRUE)), primer_count_table[j, sample])
      }
      
      remove(file_contents)
      
      counter = counter + 1
    }
    
    return(primer_count_table)
  }
  
  # Normalize primer counts
  reduce_and_normalize_primer_counts <- function(primer.count.table, sample) {
    
    
    # Adding forward and reverse primer counts
    
    for(i in 1:length(primer.count.table[,1])) {
      primer.count.table[rownames(primer.count.table)[i],] <- primer.count.table[rownames(primer.count.table)[i],] + primer.count.table[paste0(rownames(primer.count.table)[i],"rc"),]
    }
    
    primer_rows <- rownames(primer.count.table)
    primer_rows <- primer_rows[1:(length(primer_rows)/2)]
    
    
    # Making a new frame from added values
    
    reduced_table <- as.data.frame(primer.count.table[1:(length(primer.count.table[,1])/2),], row.names = primer_rows)
    colnames(reduced_table) <- sample
    
    loci_list <- unlist(strsplit(rownames(reduced_table), "_"))
    loci_list <- unique(loci_list[grep(">", loci_list)])
    
    normalized_table <- data.frame(matrix(0, nrow=length(loci_list), ncol=length(sample)))
    rownames(normalized_table) <- loci_list
    colnames(normalized_table) <- sample
    
    
    # Build further reduced table by combining counts for each locus
    for(j in 1:length(loci_list)) {
      normalized_table[loci_list[j], sample] <- sum(reduced_table[grep(paste0(loci_list[j],"_"), rownames(reduced_table)), sample])
    }
    
    # Taking out all the arrows from the loci list
    new_list <- unique(unlist(strsplit(loci_list, ">")))
    new_list <- new_list[2:length(new_list)]
    rownames(normalized_table) <- sapply(strsplit(rownames(normalized_table), ">"), '[', 2)
    
    return(normalized_table)
  }
  
  # Check against threshold
  threshold_check <- function(normalized.table, sample, threshold.kff){
    
    # Checking to see if the results are above the set threshold value
    normalized.table[, sample] <- (normalized.table[, sample] > threshold.kff) * 1
    
    return(normalized.table)
  }
  
  # Allele generation for locus (relies on resource files)
  alleles_gen.vcf <- function(current.locus){
    all_alleles_preKFF <- read.delim(paste0("Resources/caller_resources/All_", current.locus, "_preKFF.fas"), header=F, colClasses="character")
    length_seq <- nchar(all_alleles_preKFF[1,2])
    num_seq <- length(all_alleles_preKFF[,1])
    top_row <- data.frame(cbind(num_seq, length_seq))
    names(top_row) <- names(all_alleles_preKFF)
    all_alleles <- rbind(top_row, all_alleles_preKFF)
    write.table(all_alleles, paste0(results.directory, "KIRcaller/All_", current.locus, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Snps generation for locus (relies on resource files)
  snps_gen.vcf <- function(current.locus){
    seq_locus <- read.dna(paste0(results.directory, "KIRcaller/All_", current.locus, ".txt"), format = "sequential")
    snps <- seg.sites(seq_locus)
    seq_locus.2 <- read.dna(paste0(results.directory, "KIRcaller/All_", current.locus, ".txt"), format = "sequential", as.character = TRUE)
    seq_locus.2 <- apply(seq_locus.2, MARGIN = c(1, 2), toupper)
    
    seq_locus.df <- data.frame(seq_locus.2, stringsAsFactors = F)
    snp.cols <- paste0("X", snps)
    snps_locus <- seq_locus.df[,snp.cols]
    allele <- gsub(">", "", row.names(seq_locus.2))
    snps_locus <- data.frame(cbind(allele, snps_locus), stringsAsFactors = F)
    write.table(snps_locus, paste0(results.directory, "KIRcaller/KIR_", current.locus, "_alleles.txt"), quote = FALSE, row.names = FALSE)
    
    SOS_locus_lookup <- snps_locus
    
    KIR_locus_gene_VScDNA <- read.delim(paste0("Resources/caller_resources/KIR", current.locus, "geneVScDNA.txt"))
    cdna <- paste0("X", KIR_locus_gene_VScDNA$cDNA)
    gene <- paste0("X", KIR_locus_gene_VScDNA$gene)
    map_pos <- data.frame(cbind(cdna, gene))
    pos_trans <- as.character(map_pos[match(names(SOS_locus_lookup[-1]), map_pos$cdna), 2])
    names(SOS_locus_lookup) <- c("allele", pos_trans)
    
    return(SOS_locus_lookup)
  }
  
  # Pull indels from VCF data frame
  indels.vcf <- function(x) {
    vcf_file <- data.frame(x)
    vcf_file_indels <- vcf_file[grep("INDEL;", vcf_file[,8]),]
    
    if(nrow(vcf_file_indels) > 0){
      vcf_file_indels[,1] <- paste0(sample, vcf_file_indels[,1], sep="_")
    }
    
    return(vcf_file_indels)
  }
  
  # Remove indels from VCF data frame
  nodels.vcf <- function(x) {
    vcf_file <- data.frame(x)
    vcf_file_nodels <- vcf_file[grep("INDEL;", vcf_file[,8], invert = TRUE),]
    return(vcf_file_nodels)
  }
  
  # Format VCF data frame for allele calling functions
  format.vcf <- function(x) {
    vcf_file <- data.frame(x)
    vcf_file[,10] <- gsub("0/0", "0", vcf_file[,10]) ### This is to convert between old and new samtools vcf format
    vcf_file <- vcf_file[grep("./.", vcf_file[,10], fixed = T, invert = T),]
    vcf_file <- vcf_file[grep("DP4=", vcf_file[,8]),]
    
    genotype <- sapply(strsplit(as.character(vcf_file[,10]), ":"), "[", 1)
    depth <- sapply(strsplit(as.character(vcf_file[,8]), ";"), "[", 1)
    
    # A little check to make sure we are actually getting the depth
    if(!isTRUE(length(grep("DP=", depth)) > 0)) {
      stop("Bad VCF format, DP is not where expected.")
    }
    
    depth <- as.numeric(gsub("DP=", "", depth))
    
    for_rev <- unlist(strsplit(as.character(vcf_file[,8]), ";"))
    for_rev <- for_rev[grep("DP4=", for_rev)]
    fw <- as.numeric(sapply(strsplit(as.character(for_rev), ","), "[", 3))
    rev <- as.numeric(sapply(strsplit(as.character(for_rev), ","), "[", 4))
    fw[is.na(fw)] <- 0
    rev[is.na(rev)] <- 0
    fw <- ifelse(fw > 3, 1, 0)
    rev <- ifelse(rev > 3, 1, 0)
    fwrev <- fw + rev
    
    KIR_sample <- data.frame(cbind(as.character(vcf_file[,2]), as.character(vcf_file[,4]), as.character(vcf_file[,5]), genotype, as.character(depth), fwrev), stringsAsFactors = FALSE)
    KIR_sample <- KIR_sample[!KIR_sample[,1] == "line", ]
    row.names(KIR_sample) <- KIR_sample[,1]
    
    KIR_sample_snps <- KIR_sample[as.character(positions),]
    
    KIR_sample_snps <- na.omit(KIR_sample_snps)
    names(KIR_sample_snps) <- c("position", "ref", "var", "genotype", "depth", "fwrev")
    
    KIR_sample_snps$depth <- as.numeric(KIR_sample_snps$depth)
    KIR_sample_snps <- KIR_sample_snps[!KIR_sample_snps$depth < 20,]
    
    if(length(KIR_sample_snps[,1]) == 0){
      return(NULL)
    }
    
    KIR_sample_snps$snp1 <- sapply(strsplit(as.character(KIR_sample_snps$genotype), "/"), "[", 1)
    KIR_sample_snps$snp2 <- sapply(strsplit(as.character(KIR_sample_snps$genotype), "/"), "[", 2)
    KIR_sample_snps[is.na(KIR_sample_snps)] <- 0
    
    KIR_sample_snps$var_called <- as.numeric(KIR_sample_snps$snp1)+as.numeric(KIR_sample_snps$snp2)
    KIR_sample_snps$var_called <- ifelse(KIR_sample_snps$var_called>0,1,0)
    KIR_sample_snps$var_called <- ifelse(KIR_sample_snps$var_called>0,1,0)
    
    KIR_sample_snps$fwrev <- ifelse(KIR_sample_snps$fwrev>0,1,0)
    
    return(KIR_sample_snps)
  }
  
  # Even more VCF formatting
  emformat.vcf <- function(x){
    KIR_sample_snps <- format.vcf(x)
    
    if(is.null(KIR_sample_snps)){
      return(NULL)
    }
    
    KIR_sample_snps$badcall <- KIR_sample_snps$var_called+KIR_sample_snps$fwrev
    KIR_sample_snps <- KIR_sample_snps[!KIR_sample_snps$badcall == 1, ]
    
    if(nrow(KIR_sample_snps) == 0){
      return(NULL)
    }
    
    KIR_sample_snps$genotype <- NULL
    
    KIR_sample_snps <- KIR_sample_snps[!KIR_sample_snps$snp1>2, ]
    KIR_sample_snps <- KIR_sample_snps[!KIR_sample_snps$snp2>2, ]
    
    for( i in 1:length(KIR_sample_snps$snp1) ) {
      KIR_sample_snps$snp1call[i] <- ifelse(KIR_sample_snps$snp1[i] == 0, as.character(KIR_sample_snps$ref[i]), as.character(str_split(KIR_sample_snps$var[i], ",")[[1]][as.numeric(KIR_sample_snps$snp1[i])]))
    }
    
    for( i in 1:length(KIR_sample_snps$snp2) ) {
      KIR_sample_snps$snp2call[i] <- ifelse(KIR_sample_snps$snp2[i] == 0, as.character(KIR_sample_snps$ref[i]), as.character(str_split(KIR_sample_snps$var[i], ",")[[1]][as.numeric(KIR_sample_snps$snp2[i])]))
    }
    
    return(KIR_sample_snps)
  }
  
  # Check VCF data frame for correct formatting
  check.vcf <- function(x){
    
    # if(length(x[,1]) == 383 && current.locus == "2DL23"){
    #   # Generating position information for the current locus
    #   SOS_locus_lookup <<- snps_gen.vcf("2DL2")
    #   
    #   positions <<- gsub("X", "", names(SOS_locus_lookup)[-1])
    # }
    
    KIR_sample_snps <- emformat.vcf(x)
    
    # if(length(x[,1]) == 383 && current.locus == "2DL23"){
    #   # Generating position information for the current locus
    #   SOS_locus_lookup <<- snps_gen.vcf(current.locus)
    #   
    #   positions <<- gsub("X", "", names(SOS_locus_lookup)[-1])
    # }
    
    good_vcf <- ifelse(dim(KIR_sample_snps)[1] > 2, "yes", "no")
    
    if(is.null(KIR_sample_snps)){
      good_vcf <- "no"
    }
    
    return(good_vcf)
  }
  
  # Function to enumerate all possible haplotypes of SNPs (all possible genotypes)
  # cannibalized from haplo.stats
  haplo.enum <- function(hmat, geno_call, reads_pos) {
    # Exact same as haplo.enum but made for handling one whole row of ugeno
    # instead of haploids
    
    # enumerate all possible haplotypes, given input vectors
    # h1 and h2 (one possible set of haplotypes)
    # and return in matrices h1 and h2 all possible haps.
    
    # This algorithm sets up the h1.mtx and h2.mtx matrices with NA
    # values, then moves across the loci which are heterozygous 
    # (after the 1st het locus), flipping alleles at heterozygous
    # locations, and saving results in rows moving down the matrices.
    
    # next three commands added 8/23 jps
    
    old_hmat <- hmat
    
    hmat <- matrix(hmat,nrow=1)
    h1 <- hmat[,seq(1,(ncol(hmat)-1),by=2)]
    h2 <- hmat[,seq(2,ncol(hmat),by=2)]
    # the rest is the same as before
    het  <- h1 != h2
    nhet <- sum(het)  # only enumerate if > 1 het loci
    
    # Trying to work out a shortcut for highly ambiguous calls
    if(nhet >= 10){
      lookup <- SOS_locus_lookup
      lookup <- lookup[reads_pos]
      for(i in 1:ncol(old_hmat)){
        lookup <- lookup[as.logical(lookup[i] == old_hmat[1,i] | lookup[i] == old_hmat[2,i]),]
        
        if(nrow(lookup) == 0 && geno_call == T){
          cat(paste("\n\nBroke out of hmat! This most likely saved a significant amount of time\n\n"))
          return(list(h1 = matrix(h1, nrow = 1), h2 = matrix(h2, nrow = 1
          )))
        }
      }
    }
    
    if(nhet <= 1) {
      return(list(h1 = matrix(h1, nrow = 1), h2 = matrix(h2, nrow = 1
      )))
    }
    
    # only need to flip at heterozygous loci, after 1st het locus:
    nloci <- length(h1)
    which <- (1:nloci)[het]
    which <- which[-1]
    h1.mtx <- h2.mtx <- matrix(NA, nrow = 2^(nhet - 1), ncol = nloci)
    h1.mtx[1,  ] <- h1
    h2.mtx[1,  ] <- h2
    indx.row <- 1
    for(i in which) {
      nr <- sum(!is.na(h1.mtx[, 1]))
      for(j in 1:nr) {
        indx.row <- indx.row + 1	
        # used to move down to the next row of matrix
        # now for flipping alleles across loci
        if(i < nloci) {
          h1.mtx[indx.row,  ] <- c(h1.mtx[j, 1:(i - 1)], 
                                   h2.mtx[j, i], h1.mtx[j, (i + 1):nloci])
          h2.mtx[indx.row,  ] <- c(h2.mtx[j, 1:(i - 1)], 
                                   h1.mtx[j, i], h2.mtx[j, (i + 1):nloci])
        }
        if(i == nloci) {
          h1.mtx[indx.row,  ] <- c(h1.mtx[j, 1:(i - 1)], 
                                   h2.mtx[j, i])
          h2.mtx[indx.row,  ] <- c(h2.mtx[j, 1:(i - 1)], 
                                   h1.mtx[j, i])
        }
      }
    }
    
    if(nhet >= 10 && geno_call == T){
      for(i in which[-1]){
        lookup_strings <- do.call(paste0, lookup[,1:i])
        
        good_haps_h1 <- apply(format(h1.mtx[,1:i]), 1, paste, collapse="") %in% lookup_strings
        good_haps_h2 <- apply(format(h2.mtx[,1:i]), 1, paste, collapse="") %in% lookup_strings
        h1.mtx <- h1.mtx[good_haps_h1 & good_haps_h2,, drop = FALSE]
        h2.mtx <- h2.mtx[good_haps_h1 & good_haps_h2,, drop = FALSE]
        
        if(nrow(h1.mtx) == 0){
          error_string <- paste("\n\nUh ohh, all the rows got taken out trying to find new alleles for", sample, "at", current.locus, "\n\n")
          cat(error_string)
          error.log(error_string)
          return(list(h1 = rbind(h1.mtx, "N"), h2 = rbind(h2.mtx, "N")))
        }
      }
    }else if(nhet >= 10){
      for(i in which[-1]){
        lookup_strings <- do.call(paste0, lookup[,1:i])
        
        good_haps_h1 <- apply(format(h1.mtx[,1:i]), 1, paste, collapse="") %in% lookup_strings
        good_haps_h2 <- apply(format(h2.mtx[,1:i]), 1, paste, collapse="") %in% lookup_strings
        h1.mtx <- h1.mtx[good_haps_h1 | good_haps_h2,, drop = FALSE]
        h2.mtx <- h2.mtx[good_haps_h1 | good_haps_h2,, drop = FALSE]
        
        if(nrow(h1.mtx) == 0){
          error_string <- paste("\n\nUh ohh, all the rows got taken out trying to find new alleles for", sample, "at", current.locus, "\n\n")
          cat(error_string)
          error.log(error_string)
          return(list(h1 = rbind(h1.mtx, "N"), h2 = rbind(h2.mtx, "N")))
        }
      }
    }
    
    return(list(h1 = h1.mtx, h2 = h2.mtx))
  }
  
  # Special 2DL23 genotyping
  # THIS SHOULD BE MOVED TO POST PROCESSING!!
  # genos_2DL23 <- function(genos.out){
  #   
  #   vector_call_2DL2 <- NA
  #   vector_call_2DL3 <- NA
  #   
  #   for(item in genos.out){
  #     
  #     item_frame <- as.data.frame(item)
  #     
  #     # Grabbing and reshaping each type of call into a vector
  #     call_2DL2 <- grep("KIR2DL2_", item_frame)
  #     if(length(call_2DL2) > 0){
  #       frame_call_2DL2 <- item_frame[call_2DL2]
  #       
  #       vector_call_2DL2 <- c(vector_call_2DL2, t(frame_call_2DL2))
  #       vector_call_2DL2 <- vector_call_2DL2[!is.na(vector_call_2DL2)]
  #     }
  #     
  #     call_2DL3 <- grep("KIR2DL3_", item_frame)
  #     if(length(call_2DL3) > 0){
  #       frame_call_2DL3 <- item_frame[call_2DL3]
  #       
  #       vector_call_2DL3 <- c(vector_call_2DL3, t(frame_call_2DL3))
  #       vector_call_2DL3 <- vector_call_2DL3[!is.na(vector_call_2DL3)]
  #     }
  #   }
  #   
  #   
  #   # Finding any 2DL2 alleles by finding the intersection of each 2DL2 call
  #   if(length(vector_call_2DL2) == 3){
  #     one <- unlist(strsplit(vector_call_2DL2, "/")[1])
  #     two <- unlist(strsplit(vector_call_2DL2, "/")[2])
  #     three <- unlist(strsplit(vector_call_2DL2, "/")[3])
  #     
  #     step_one <- intersect(two, three)
  #     allele_one <- intersect(one, step_one)
  #     
  #   }else if(length(vector_call_2DL2) == 4){
  #     one <- unlist(strsplit(vector_call_2DL2, "/")[1])
  #     two <- unlist(strsplit(vector_call_2DL2, "/")[2])
  #     three <- unlist(strsplit(vector_call_2DL2, "/")[3])
  #     four <- unlist(strsplit(vector_call_2DL2, "/")[4])
  #     
  #     allele_one <- intersect(one, three)
  #     allele_two <- intersect(two, four)
  #     
  #     if(length(allele_one) == 0 || length(allele_two) == 0){
  #       allele_one <- intersect(one, four)
  #       allele_two <- intersect(two, three)
  #     }
  #   }
  #   
  #   # Finding any 2DL3 alleles by finding the intersection of each 2DL3 call
  #   if(length(vector_call_2DL3) == 3){
  #     one <- unlist(strsplit(vector_call_2DL3, "/")[1])
  #     two <- unlist(strsplit(vector_call_2DL3, "/")[2])
  #     three <- unlist(strsplit(vector_call_2DL3, "/")[3])
  #     
  #     step_one <- intersect(two, three)
  #     allele_two <- intersect(one, step_one)
  #     
  #   }else if(length(vector_call_2DL3) == 4){
  #     one <- unlist(strsplit(vector_call_2DL3, "/")[1])
  #     two <- unlist(strsplit(vector_call_2DL3, "/")[2])
  #     three <- unlist(strsplit(vector_call_2DL3, "/")[3])
  #     four <- unlist(strsplit(vector_call_2DL3, "/")[4])
  #     
  #     allele_one <- intersect(one, three)
  #     allele_two <- intersect(two, four)
  #     
  #     if(length(allele_one) == 0 || length(allele_two) == 0){
  #       allele_one <- intersect(one, four)
  #       allele_two <- intersect(two, three)
  #     }
  #   }
  #   
  #   
  #   ## BUG FIX 7/18 new genotypes were not being caught by the tryCatch because they were character(0), so I am removing objects
  #   if(exists("allele_one") && length(allele_one) == 0){
  #     remove(allele_one)
  #   }
  #   if(exists("allele_two") && length(allele_two) == 0){
  #     remove(allele_two)
  #   }
  #   
  #   genos.out <- genos.out[1]
  #   genos.out[[1]]$X1 <- tryCatch(paste(allele_one, collapse = "/"), error=function(e) "new")
  #   genos.out[[1]]$X2 <- tryCatch(paste(allele_two, collapse = "/"), error=function(e) "new")
  #   
  #   return(genos.out)
  # }
  
  error.log <- function(error.string){
    error_file <- paste0(results.directory, "Error.log")
    
    if(file.exists(error_file)){
      error_log <- file(error_file, open = "at")
    }else{
      error_log <- file(error_file, open = "wt")
    }
    
    sink(error_log)
    cat(error.string)
    sink()
    close(error_log)
  }
  
  
  #################
  # PING functions ----------------------------------------------------------
  #################
  
  
  # Determine presence/absence by counting probe matches
  ping.kff <- function(sample, current.locus, primerlist, threshold) {
    
    # Run KFF -------------------------------------------------
    
    cat("\n\n")
    cat("----- Running KFF -----\n")
    
    cat("\n")
    cat("\n")
    cat("Building primer table. \n")
    primer_table <- make_primer_table(primerlist)
    
    cat("\n")
    cat("Counting primers in: \n")
    primer_match_table <- count_primer_matches(primer_table, sample, is_gz)
    
    cat("\n")
    cat("Reducing and normalizing primer counts. \n")
    normalized_results <- reduce_and_normalize_primer_counts(primer_match_table, sample)
    
    if(file.exists(paste0(results.directory, "KIRcaller/KFF_counts_", current.locus, ".txt"))){
      write.table(t(normalized_results), paste0(results.directory, "KIRcaller/KFF_counts_", current.locus, ".txt"), sep = " ", append = T, row.names = T, col.names = F, quote = F)
    }else{
      write.table(t(normalized_results), paste0(results.directory, "KIRcaller/KFF_counts_", current.locus, ".txt"), sep = " ", row.names = T, col.names = NA, quote = F)
    }
    
    checked_results <- threshold_check(normalized_results, sample, threshold)
    
    
    ## If the results file already exists, just add on to the bottom of it
    
    transposed_results <- data.frame(t(checked_results), check.names = F)  
    
    return(transposed_results)
  }
  
  # Integrate kff results with genotype calls
  kff_integrator <- function(sample, current.locus, lookitup, poss_genos){
    
    ###alle calling
    allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 1]),"new")
    allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 1]),"new")
    poss_genos <- data.frame(cbind(as.character(allele1),as.character(allele2)), stringsAsFactors = FALSE)
    poss_genos <- poss_genos[!poss_genos$X1 == "new",]
    poss_genos <- poss_genos[!poss_genos$X2 == "new",]
    
    kff_legend <- read.delim(paste0("Resources/caller_resources/KIR_", current.locus, "_kff_legend.txt"))
    kff_types <- kff_legend[grep("KIR", kff_legend[,1], invert=T),]
    kff_legend <- kff_legend[grep("KIR", kff_legend[,1]),]
    
    if(length(kff_legend$allele) != length(lookitup$string)) {
      cat("\n\n!!\nError: Number of kff legend alleles does not equal All_*_preKFF.fas alleles.\n!!\n\n")
    }
    
    
    ## Generate kff results to poss_genos table
    kff_results <- ping.kff(sample, current.locus, "Resources/caller_resources/2DL49or10_64bit.txt", 10)
    
    
    ## Get allele codes for kff positive results
    codes <- NULL
    rownum <- grep(sample, rownames(kff_results))
    
    for(i in 1:length(kff_results[rownum,])){
      if(as.numeric(kff_results[rownum, i]) == 1){
        codes <- c(codes, as.character(kff_types$code[grep(colnames(kff_results)[i], kff_types$allele)]))
      }
    }
    
    
    ## Permutate those codes
    if(length(codes) == 0){
      
      ## Kff found no results, falling back on the original allele calling method
      allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 3]),"new")
      allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 3]),"new")
      poss_genos <- data.frame(cbind(as.character(allele1),as.character(allele2)), stringsAsFactors = FALSE)
      
      error.message <- paste0("\n\nERROR: COULD NOT FIND ADEQUATE KFF MATCHES FOR ", sample, " at ", current.locus, ". FALLING BACK ON NON-KFF ALLELE CALLING.!!\n\n")
      cat(error.message)
      error.log(error.message)
      
      return(poss_genos)
      
    }else if(length(codes) == 1){
      code_frame <- data.frame(cbind(codes[1], codes[1]))
    }else{
      code_frame <- permutator(codes)
    }
    
    ## Add the codes to poss_genos
    frame_length <- length(code_frame$X1)*length(poss_genos$X1)
    new_genos <- data.frame(matrix(nrow = frame_length, ncol = 2))
    
    for(j in 0:(length(poss_genos$X1) - 1)){
      for(k in 1:length(code_frame$X1)){
        new_genos[j*length(code_frame$X1) + k, 1] <- paste0(poss_genos[j + 1, 1], code_frame[k, 1])
        new_genos[j*length(code_frame$X1) + k, 2] <- paste0(poss_genos[j + 1, 2], code_frame[k, 2])
      }
    }
    poss_genos <- new_genos
    
    
    ## Mutating the lookitup table
    for(i in 1:length(kff_legend$allele)){
      lookitup[grep(kff_legend[i, "allele"], rownames(lookitup)), "string"] <- paste0(lookitup[grep(kff_legend[i, "allele"], rownames(lookitup)), "string"], as.character(kff_legend[i, "code"]))
    }
    
    
    ## Running allele calling again with the new codes
    allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 3]),"new")
    allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 3]),"new")
    poss_genos <- data.frame(cbind(as.character(allele1),as.character(allele2)), stringsAsFactors=FALSE)
    
    return(poss_genos)
  }
  
  
  
  # Call genotypes from VCF data frame
  allele_call.vcf <- function(x, sample){
    
    # if(length(x[,1]) == 383 && current.locus == "2DL23"){
    #   # Generating position information for the current locus
    #   SOS_locus_lookup <<- snps_gen.vcf("2DL2")
    #   
    #   positions <<- gsub("X", "", names(SOS_locus_lookup)[-1])
    # }
    
    KIR_sample_snps <- emformat.vcf(x)
    
    if(is.null(KIR_sample_snps)){
      return(NULL)
    }
    
    reads_pos <- paste0("X", KIR_sample_snps$position)
    
    hmat <- as.matrix(rbind(KIR_sample_snps$snp1call,KIR_sample_snps$snp2call))
    poss_haps <- haplo.enum(hmat, T, reads_pos)
    poss_genos <- data.frame(cbind(apply(poss_haps$h1, MARGIN=1, paste, collapse=""), apply(poss_haps$h2, MARGIN=1, paste, collapse="")))
    poss_genos$X1 <- as.character(poss_genos$X1)
    poss_genos$X2 <- as.character(poss_genos$X2)
    
    
    
    ##get right loookup table based on which snps we have calls for in VCF file
    SOS_locus_lookup_reads<-SOS_locus_lookup[,c("allele", reads_pos)]
    lookitup <- data.frame(cbind(apply(SOS_locus_lookup_reads[,-1], MARGIN=1, paste, collapse=""), as.character(SOS_locus_lookup$allele)))
    names(lookitup) <- c("string", "allele")
    lookitup$allele <- as.character(lookitup$allele)
    lookitup$string <- as.character(lookitup$string)
    
    
    ##find allelic ambiguities in current lookup table
    
    all_amb <- c(1:length(lookitup))
    for(i in seq_along(lookitup$string)) {
      l <- str_detect(lookitup$string,lookitup$string[i])
      ll <- lookitup$allele[l]
      
      all_amb[i] <- paste(ll, collapse="/")
    }
    lookitup <- cbind(lookitup,as.character(all_amb), stringsAsFactors=FALSE)
    names(lookitup) <- c("string","allele","amb")
    
    
    ## Mutating the lookup table to accomodate kff results 
    if("2DL4" == current.locus){
      
      poss_genos <- kff_integrator(sample, current.locus, lookitup, poss_genos)
      
    }else if("2DS4" == current.locus){
      
      kff_results <- ping.kff(sample, current.locus, "Resources/caller_resources/KFF_2DS4.txt", 10)
      
      ###alle calling
      allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 3]),"new")
      allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 3]),"new")
      poss_genos <- data.frame(cbind(as.character(allele1),as.character(allele2)), stringsAsFactors=FALSE)
      poss_genos <- poss_genos[!poss_genos$X1 == "new",]
      poss_genos <- poss_genos[!poss_genos$X2 == "new",]
      
      
      ## Finding kff positive allele names
      kff_positive <- colnames(kff_results)[kff_results == 1]
      
      ## Mutating allele names to match poss_genos names
      kff_positive <- gsub(current.locus, paste0("KIR", current.locus, "_"), kff_positive)
      
      ## Check for del variant
      kff_del <- "KIR2DS4_del" %in% kff_positive
      
      if(kff_del){
        kff_positive <- kff_positive[kff_positive != "KIR2DS4_del"]
        if(!length(grep("del", kff_positive)) > 0){
          error.string <- paste0("\n\nERROR: Positive KFF del variant hit for ", sample, " at ", current.locus, ", but no del allele found.\n\n")
          cat(error.string)
          error.log(error.string)
        }
      }
      
      replacement_string <- "KIR2DS4_00101/KIR2DS4_00102/KIR2DS4_00103/KIR2DS4_003/KIR2DS4_006/KIR2DS4_009"
      
      if(length(kff_positive) == 0){
        error.string <- paste0("\n\nERROR: Copy number of 1 detected, but no kff matches found for ", sample, " at ", current.locus, "\n\n")
        cat(error.string)
        error.log(error.string)
      }else if(length(grep("1", copy_number)) > 0){
        if(length(grep(replacement_string, poss_genos$X1)) > 0){
          poss_genos$X1 <- gsub(replacement_string, kff_positive, poss_genos$X1)
          poss_genos$X2 <- gsub(replacement_string, kff_positive, poss_genos$X2)
        }
      }else if(length(grep("2", copy_number)) > 0){
        if(length(grep(replacement_string, poss_genos$X1)) > 0 || length(grep(replacement_string, poss_genos$X2)) > 0){
          
          if(length(kff_positive) == 1){
            poss_genos$X1 <- gsub(replacement_string, kff_positive, poss_genos$X1)
            poss_genos$X2 <- gsub(replacement_string, kff_positive, poss_genos$X2)
          }else if(length(kff_positive) == 2){
            poss_genos$X1 <- gsub(replacement_string, kff_positive[1], poss_genos$X1)
            poss_genos$X2 <- gsub(replacement_string, kff_positive[2], poss_genos$X2)
          }else{
            error.string <- paste0("\n\nERROR: Copy number of 2 detected, but KFF hits do not match for ", sample, " at ", current.locus, "\n\n")
            cat(error.string)
            error.log(error.string)
          }
        }
      }else{
        error.string <- paste0("\n\nERROR: Copy number > 2 detected, please perform KFF refinement by hand for ", sample, " at ", current.locus, "\n\n")
        cat(error.string)
        error.log(error.string)
      }
      
    }else if("3DL1" == current.locus || "3DS1" == current.locus){
      
      kff_results <- ping.kff(sample, current.locus, "Resources/caller_resources/KFF_3DL1.txt", 10)
      
      ###alle calling
      allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 3]),"new")
      allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 3]),"new")
      poss_genos <- data.frame(cbind(as.character(allele1),as.character(allele2)), stringsAsFactors=FALSE)
      poss_genos <- poss_genos[!poss_genos$X1 == "new",]
      poss_genos <- poss_genos[!poss_genos$X2 == "new",]
      
      kff_legend <- read.delim("Resources/caller_resources/KIR_3DL1S1_kff_legend.txt", header = F, fill = T)
      
      kff_positive <- colnames(kff_results[which(kff_results[1,] == 1)])
      
      if(length(kff_positive) == 0){
        genos <- poss_genos
        return(genos)
      }
      
      legend_rows <- grep(paste0(kff_positive, collapse = "|"), kff_legend[,1])
      legend_antirows <- grep(paste0(kff_positive, collapse = "|"), kff_legend[,1], invert = T)
      kff_poss_genos <- kff_legend[legend_rows,]
      
      
      if(length(unlist(poss_genos)) == 0){
        kff_poss_genos <- kff_poss_genos[,2:length(kff_poss_genos)]
        kff.string <- paste0(kff_poss_genos[!kff_poss_genos[,1:length(kff_poss_genos)] == ""], collapse = ", ")
        
        error.string <- paste0("\n\nERROR: KFF probe matches found, but no 3DL1/S1 genotype called for sample: ", sample, ". KFF probes matched for ", kff.string)
        cat(error.string)
        error.log(error.string)
      }else if(length(legend_antirows) > 0){
        kff_anti_genos <- kff_legend[legend_antirows,]
        kff_anti_genos <- kff_anti_genos[,2:length(kff_anti_genos)]
        
        if(any(unlist(poss_genos) %in% kff_anti_genos[kff_anti_genos[,1:length(kff_anti_genos)] != ""])){
          error.string <- paste0("\n\nERROR: 3DL1/S1 genotype called for alleles not found by KFF probes for sample: ", sample)
          cat(error.string)
          error.log(error.string)
        }
      }
      
    }else{
      ###alle calling
      allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 3]),"new")
      allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 3]),"new")
      poss_genos <- data.frame(cbind(as.character(allele1),as.character(allele2)), stringsAsFactors=FALSE)
      poss_genos <- poss_genos[!poss_genos$X1 == "new",]
      poss_genos <- poss_genos[!poss_genos$X2 == "new",]
    }
    
    ## If copy number is one, poss_genos exists, and X1 and X2 are equal, change X2 to a null allele
    # if( copy_number == 1 && length(poss_genos[,1]) !=0 && poss_genos$X1 == poss_genos$X2 && current.locus != "2DL23"){
    #   poss_genos$X2 <- paste0("KIR", current.locus, "_null")
    # }
    
    # if(length(x[,1]) == 383 && current.locus == "2DL23"){
    #   # Generating position information for the current locus
    #   SOS_locus_lookup <<- snps_gen.vcf(current.locus)
    #   
    #   positions <<- gsub("X", "", names(SOS_locus_lookup)[-1])
    # }
    
    genos <- poss_genos
    genos
  }
  
  # Call possible new snps from VCF data frame
  new_snps_call.vcf <- function(x){
    vcf_file <- data.frame(x)
    vcf_file[,10] <- gsub("0/0", "0", vcf_file[,10]) ### This is to convert between old and new samtools vcf format
    vcf_file <- vcf_file[grep("./.", vcf_file[,10], fixed = T, invert = T),]
    vcf_file <- vcf_file[grep("DP4=", vcf_file[,8]),]
    
    genotype <- sapply(strsplit(as.character(vcf_file[,10]), ":"), "[", 1)
    depth <- sapply(strsplit(as.character(vcf_file[,8]), ";"), "[", 1)
    
    # A little check to make sure we are actually getting the depth
    if(!isTRUE(length(grep("DP=", depth)) > 0)) {
      stop("Bad VCF format, DP is not where expected.")
    }
    
    depth <- as.numeric(gsub("DP=", "", depth))
    
    for_rev <- unlist(strsplit(as.character(vcf_file[,8]), ";"))
    for_rev <- for_rev[grep("DP4=", for_rev)]
    fw <- as.numeric(sapply(strsplit(as.character(for_rev), ","), "[", 3))
    rev <- as.numeric(sapply(strsplit(as.character(for_rev), ","), "[", 4))
    fw[is.na(fw)] <- 0
    rev[is.na(rev)] <- 0
    fw <- ifelse(fw > 3, 1, 0)
    rev <- ifelse(rev > 3, 1, 0)
    fwrev <- fw + rev
    
    KIR_sample <- data.frame(cbind(as.character(vcf_file[,2]), as.character(vcf_file[,4]), as.character(vcf_file[,5]), genotype, as.character(depth), fwrev), stringsAsFactors = FALSE)
    KIR_sample <- KIR_sample[!KIR_sample[,1] == "line", ]
    row.names(KIR_sample) <- KIR_sample[,1]
    
    names(KIR_sample)<-c("position","ref","var","genotype","depth","fwrev")
    
    KIR_sample$depth<-as.numeric(KIR_sample$depth)
    KIR_sample<-KIR_sample[!KIR_sample$depth < 20,]
    
    KIR_sample$snp1<-sapply(strsplit(as.character(KIR_sample$genotype), "/"), "[", 1)
    KIR_sample$snp2<-sapply(strsplit(as.character(KIR_sample$genotype), "/"), "[", 2)
    
    KIR_sample[ is.na(KIR_sample) ] <- 0
    
    KIR_sample$var_called<-as.numeric(KIR_sample$snp1)+as.numeric(KIR_sample$snp2)
    KIR_sample$var_called<-ifelse(KIR_sample$var_called>0,1,0)
    KIR_sample$var_called<-ifelse(KIR_sample$var_called>0,1,0)
    
    KIR_sample<-KIR_sample[KIR_sample$var_called==1,]
    
    KIR_sample_variants<-KIR_sample[KIR_sample$var_called==1,]
    
    new_snps<-KIR_sample_variants[!KIR_sample_variants$position %in% positions,]
    new_snps$fwrev<-NULL
    new_snps$snp1<-NULL
    new_snps$snp2<-NULL
    new_snps$var_called<-NULL
    
    dp4_lines <- str_split(data.frame(x)[data.frame(x)[,2] %in% new_snps$position,8], ";")
    dp4_values <- unlist(lapply(dp4_lines, grep, pattern="DP4", value = TRUE))
    new_snps$DP4 <- dp4_values
    
    new_snps
  }
  
  # Call possible new alleles from VCF data frame
  new_alleles.vcf <- function(x){
    KIR_sample_snps <- emformat.vcf(x)
    
    if(is.null(KIR_sample_snps)){
      return(NULL)
    }
    
    reads_pos <- paste0("X", KIR_sample_snps$position)
    
    ##enumerate all possible genotypes
    hmat <- as.matrix(rbind(KIR_sample_snps$snp1call, KIR_sample_snps$snp2call))
    poss_haps <- haplo.enum(hmat, F, reads_pos)
    poss_genos <- data.frame(cbind(apply(poss_haps$h1, MARGIN=1, paste, collapse=""), apply(poss_haps$h2, MARGIN=1, paste, collapse="")))
    poss_genos$X1 <- as.character(poss_genos$X1)
    poss_genos$X2 <- as.character(poss_genos$X2)
    
    ##get right loookup table based on which snps we have calls for in VCF file
    SOS_locus_lookup_reads <- SOS_locus_lookup[,c("allele", reads_pos)]
    lookitup <- data.frame(cbind(apply(SOS_locus_lookup_reads[,-1], MARGIN=1, paste, collapse=""), as.character(SOS_locus_lookup$allele)))
    names(lookitup) <- c("string","allele")
    lookitup$allele <- as.character(lookitup$allele)
    lookitup$string <- as.character(lookitup$string)
    ##find allelic ambiguities in current lookup table
    
    all_amb<-c(1:length(lookitup))
    for(i in seq_along(lookitup$string)) {
      l <- str_detect(lookitup$string,lookitup$string[i])
      ll <- lookitup$allele[l]
      
      all_amb[i] <- paste(ll, collapse="/")
    }
    
    lookitup <- cbind(lookitup,as.character(all_amb), stringsAsFactors=FALSE)
    names(lookitup) <- c("string", "allele", "amb")
    
    ###allle calling
    allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 3]),"new")
    allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 3]),"new")
    genos <- data.frame(cbind(as.character(allele1), as.character(allele2)))
    
    ## EXACTLY THE SAME AS KIR.CALL UNTIL THIS POINT ^^^^
    
    genos$X1X2 <- apply(genos, 1, paste, collapse="")
    genos$string1 <- poss_genos$X1
    genos$string2 <- poss_genos$X2
    genos <- genos[!genos$X1X2=="newnew",]
    
    if(is.data.frame(genos) && nrow(genos)==0) {
      return(NULL)
    }
    
    genos$X1X2 <- NULL
    genos$newallele <- ifelse(genos$X1=="new", as.character(genos$string1), as.character(genos$string2))
    genos$calledallele <- ifelse(genos$X1=="new", as.character(genos$X2), as.character(genos$X1))
    genos$X1 <- NULL
    genos$X2 <- NULL
    genos$string1 <- NULL
    genos$string2 <- NULL
    res.df <- data.frame(adist(genos$newallele,lookitup$string))
    names(res.df) <- lookitup$allele
    oneoff.df <- data.frame(which(res.df==1,arr.ind=TRUE))
    oneoff.df <- oneoff.df[order(oneoff.df$row),]
    
    for(i in seq_along(oneoff.df$row)) {
      x <- oneoff.df$row[i]
      oneoff.df$known[i] <- genos$calledallele[x]
    }
    
    for(i in seq_along(oneoff.df$col)) {
      x <- oneoff.df$col[i]
      oneoff.df$near[i] <- names(res.df)[x]
    }
    
    knownwithsomething <- ifelse(length(oneoff.df$col)==0,genos$calledallele, "whatev")
    something <- rep("something", length(knownwithsomething))
    something.df <- cbind(knownwithsomething, something)
    ifelse(length(oneoff.df$col)==0, oneoff.df <- something.df, oneoff.df <- oneoff.df)
    oneoff.df <- data.frame(oneoff.df)
    oneoff.df$col <- NULL
    oneoff.df$row <- NULL
    names(oneoff.df) <- c("known.allele", "new.allele.near")
    oneoff.df
  }
  
  # Converts calls to Ucode
  get_calls.vcf <- function(x){
    KIR_sample_snps <- emformat.vcf(x)
    
    KIR_sample_snps$snpcalls <- ifelse(KIR_sample_snps$snp1call == KIR_sample_snps$snp2call, as.character(KIR_sample_snps$snp1call), paste(KIR_sample_snps$snp1call, KIR_sample_snps$snp2call, sep="/"))
    codes <- c("A/G","G/A","C/T", "T/C","A/C","C/A","A/T","T/A","C/G","G/C","T/G","G/T")
    replacement <- c("R","R","Y","Y","M","M","W","W","S","S","K","K")
    
    allcalls.df <- KIR_sample_snps[, c("position", "snpcalls")]   
    
    for( i in 1:length(codes) ) {
      allcalls.df$snpcalls <- replace(allcalls.df$snpcalls, grep(codes[i], allcalls.df$snpcalls), replacement[i])
    }
    
    position.df <- data.frame(as.numeric(positions))
    colnames(position.df) <- "position"
    
    allpos.df <- merge(position.df, allcalls.df, by = "position", sort = T, all.x = TRUE)
    
    allpos.df
  }
  
  
  ######################
  # Execution functions -----------------------------------------------------
  ######################
  
  # Create a master list of samples in both the Sequence folder and PING_gc results
  master.creator <- function(sequence_list, gc_results){
    # Create master list to check against
    #master_list <- sub('\\_.*', '', sequence_list)
    master_list <- sequence_list
    
    # Cutting down master list to only include samples that have been run through PING_gc
    gc_matches <- master_list[pmatch(master_list, colnames(gc_results), nomatch = 0) != 0]
    sequence_matches <- master_list[pmatch(master_list, sequence_list, nomatch = 0) != 0]
    
    master_list <- intersect(gc_matches, sequence_matches)
    
    # Cutting sequence list and gc results table to match the master list
    sequence_list <- sequence_list[pmatch(master_list, sequence_list)]
    gc_results <- gc_results[,pmatch(master_list, colnames(gc_results))]
    
    return(master_list)
  }
  
  # VCF creation (single sample)
  ping.locus <- function(sample, current.locus, is_gz) {
    
    if("2DL1" == current.locus) {
      KIR_2DL1(sample, is_gz)
    }
    
    if("2DL23" == current.locus) {
      KIR_2DL23(sample, is_gz)
    }
    
    if("2DL4" == current.locus) {
      KIR_2DL4(sample, is_gz)
    }
    
    if("2DL5" == current.locus) {
      KIR_2DL5(sample, is_gz)
    }
    
    if("2DP1" == current.locus) {
      KIR_2DP1(sample, is_gz)
    }
    
    if("2DS3" == current.locus && !any("2DS35" %in% loci.list)){
      KIR_2DS3(sample, is_gz)
    }else if("2DS35" == current.locus){
      KIR_2DS35(sample, is_gz)
    }else if("2DS3" == current.locus){
      
    }
    
    if("2DS4" == current.locus) {
      KIR_2DS4(sample, is_gz)
    }
    
    if("3DL1" == current.locus && any("3DS1" %in% loci.list)){
      KIR_3DL1S1(sample, is_gz)
    }else if("3DL1" == current.locus){
      KIR_3DL1(sample, is_gz)
    }else if("3DS1" == current.locus && !any("3DL1" %in% loci.list)){
      KIR_3DS1(sample, is_gz)
    }
    
    if("3DL2" == current.locus) {
      KIR_3DL2(sample, is_gz)
    }
    
    if("3DL3" == current.locus) {
      KIR_3DL3(sample, is_gz)
    }
  }
  
  # Genotype caller (multi sample)
  ping.caller <- function (sample, current.locus) {
    
    # Read in vcf files and check if they are good ----------------------------
    
    vcf.location <- paste0(results.directory, "Vcf/")
    
    current_wd <- getwd()
    setwd(vcf.location)
    
    vcf_files <- list.files(pattern = paste0(sample, "_", current.locus, "nuc.vcf"))
    
    ## 2DL23 specific vcf_file finding
    # if(current.locus == "2DL23"){
    #   if(presence2DL2or3$two){
    #     vcf_files <- c(vcf_files, list.files(pattern = paste0(sample, "_", "2DL2", "nuc.vcf")))
    #   }
    #   if(presence2DL2or3$three){
    #     vcf_files <- c(vcf_files, list.files(pattern = paste0(sample, "_", "2DL3", "nuc.vcf")))
    #   }
    # }
    
    setwd(current_wd)
    
    # Return NULL if VCF files are empty (VCF files with header info only)
    for (i in vcf_files) {
      vcf_test <- tryCatch(read.table(paste0(vcf.location, i), header = F), error=function(e) NULL)
      if (is.null(vcf_test)){
        return(NULL)
      }
    }
    
    vcf_list = lapply(paste0(vcf.location, vcf_files), read.table, header = FALSE)
    names(vcf_list) <- vcf_files
    
    vcf_list_indels <- lapply(vcf_list, indels.vcf)
    
    vcf_list <- lapply(vcf_list, nodels.vcf)
    
    check_vcf.out <- lapply(vcf_list, check.vcf)
    vcf_check.df <- data.frame(as.matrix(check_vcf.out))
    
    good <- data.frame(unlist(vcf_check.df))[,1] == "yes"
    vcf_files.good = vcf_files[good]
    
    if(length(vcf_files.good) == 0){
      return(NULL)
    }
    
    vcf_list.good = lapply(paste0(vcf.location, vcf_files.good), read.table, header=F)
    vcf_list.good <- lapply(vcf_list.good, nodels.vcf)
    names(vcf_list.good) <- vcf_files.good
    
    vcf_files.bad <- data.frame(vcf_files[!good])
    
    print(class(vcf_files.bad))
    print(current.locus)
    
    write.table(vcf_files.bad, paste0(results.directory, "KIRcaller/bad_", current.locus, "_files.txt"), append = TRUE, row.names = F, col.names = F, quote = F)
    
    indel_presence <- any(lapply(vcf_list_indels, nrow) > 0)
    
    if(indel_presence) {
      if(length(vcf_list_indels) != 10){
        vcf_list_indels <- do.call(rbind, vcf_list_indels)
      }
      write.table(vcf_list_indels, paste0(results.directory, "KIRcaller/indels_", current.locus, ".txt"), append = TRUE, row.names = F, col.names = F)
    }
    
    # Run genotype caller -----------------------------------------------------
    
    genos.out <- lapply(vcf_list.good, allele_call.vcf, sample)
    
    if(all(unlist(lapply(genos.out, is.null))) || all(unlist(lapply(genos.out, nrow)) == 0)){
      cat("\nNo genotype found, moving on to look for new alleles.\n")
      genos.df <- as.data.frame(do.call(rbind, genos.out))
      genos.df <- genos.df[!genos.df$X1 == "new",]
      genos.df <- genos.df[!genos.df$X2 == "new",]
      genos.df$X1 <- NULL
      genos.df$X2 <- NULL
    }else{
      
      # if(current.locus == "2DL23"){
      #   
      #   genos.out <- genos.out[!sapply(genos.out, is.null)]
      #   
      #   # To make sure new alleles don't break anything
      #   if(!any(lapply(genos.out, nrow) == 0) && length(genos.out) > 1){
      #     
      #     # Intersecting results
      #     genos.out <- genos_2DL23(genos.out)
      #   }else if(all(lapply(genos.out, nrow)) == 0){
      #     return(NULL)
      #   }
      # }
      
      genos.df <- as.data.frame(do.call(rbind, genos.out))
      genos.df <- genos.df[!genos.df$X1 == "new",]
      genos.df <- genos.df[!genos.df$X2 == "new",]
      genos.df$genotype <- paste(genos.df$X1,genos.df$X2,sep="+")
      genos.df$X1 <- NULL
      genos.df$X2 <- NULL
      genos.df$sample <- row.names(genos.df)
      genos.df$sample <- sapply(strsplit(genos.df$sample, "_"), "[", 1)
      genos.df <- genos.df[,c(2,1)]
      var <- gsub("vcf","1",row.names(genos.df))
      genos.df$var <- sapply(strsplit(var, "nuc."), "[", 2)
      genos.wide <- reshape(genos.df, idvar = "sample", timevar="var", direction = "wide")
      genos.wide[is.na(genos.wide)]   <- " "
      
      # Asterix is used as a marker for conflicting GC results within correctly called allele output files (ex: genos_2DL4.txt)
      if(grepl("_", copy_number)){
        genos.wide[1] <- paste0(genos.wide[1], "*")
      }
      
      write.table(genos.wide, paste0(results.directory, "KIRcaller/genos_", current.locus, ".txt"), append = TRUE, row.names=F, col.names=F, quote=F)
    }
    # Run new allele caller ---------------------------------------------------
    
    got_types <- row.names(genos.df)
    good_types <- paste((unlist(lapply(strsplit(as.character(got_types), "vcf"), "[", 1))),"vcf", sep="")
    
    no_types.list <- vcf_list.good[!names(vcf_list.good) %in% good_types]
    
    if(length(no_types.list) != 0){
      
      newalleles.out <- lapply(no_types.list, new_alleles.vcf)
      
      if(!is.null(unlist(newalleles.out))){
        names(newalleles.out) <- names(no_types.list)
        newalleles.df <- as.data.frame(do.call(rbind, newalleles.out))
        newalleles.df$sample <- row.names(newalleles.df)
        newalleles.df$sample <- row.names(newalleles.df)
        newalleles.df$sample <- sapply(strsplit(newalleles.df$sample, "_"), "[", 1)
        newalleles.df <- newalleles.df[,c(3,1,2)]
        
        if(file.exists(paste0(results.directory, "KIRcaller/newalleles_", current.locus, ".txt"))) {
          write.table(newalleles.df, paste0(results.directory, "KIRcaller/newalleles_", current.locus, ".txt"), append = TRUE, row.names=F, col.names=F, quote=F)
        }else{
          write.table(newalleles.df, paste0(results.directory, "KIRcaller/newalleles_", current.locus, ".txt"), row.names=F, col.names=T, quote=F)
        }
      }
    }
    
    newsnps.out <- lapply(vcf_list.good, new_snps_call.vcf)
    newsnps.df <- as.data.frame(do.call(rbind, newsnps.out))
    
    if(file.exists(paste0(results.directory, "KIRcaller/newsnps_", current.locus, ".txt"))) {
      write.table(newsnps.df, paste0(results.directory, "KIRcaller/newsnps_", current.locus, ".txt"), append = TRUE, row.names=T, col.names=F, quote=F)
    }else{
      write.table(newsnps.df, paste0(results.directory, "KIRcaller/newsnps_", current.locus, ".txt"), row.names=T, col.names=T, quote=F)
    }
    
    
    snps.out <- lapply(vcf_list.good, get_calls.vcf)
    snps.df <- as.data.frame(do.call(cbind, snps.out))
    
    if(length(snps.df[1,]) > 2) {
      oddcols <- seq(3, ncol(snps.df), by = 2)
      snpsonly.df <- snps.df[,-oddcols]
    }else{
      snpsonly.df <- snps.df
    }
    
    colnames(snpsonly.df)[1] <- "position"
    
    if(file.exists(paste0(results.directory, "KIRcaller/snps_", current.locus, ".txt"))) {
      previous_results <- read.table(paste0(results.directory, "KIRcaller/snps_", current.locus, ".txt"), header = TRUE, sep = " ")
      snpsonly.df <- merge(previous_results, snpsonly.df, by = "position")
    }
    
    write.table(snpsonly.df, file=paste0(results.directory, "KIRcaller/snps_", current.locus, ".txt"), col.names=T, row.names=F, quote=F)
    
    return(1)
  }
  
  # Create a results directory
  results.directory <- results_directory()
  
  # Create results subfolders
  ping.ready()
  
  # Find sequences
  sequence_list <- get_sequence_list()
  
  # Are the files gzipped?
  is_gz <- last(unlist(strsplit(fastq.pattern.1, ".", fixed = T))) == "gz"
  
  # Are any sequence files too small?
  #too_small <- files_too_small(sequence_list, is_gz)
  
  # Setting a no.gz fastq pattern
  if(is_gz){
    no.gz.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz"))
    no.gz.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz"))
  }
  
  # Pull in the GC results
  gc_results <- gc.results()
  
  # The master list is the intersection of sample name matches in sequence_list and gc_results
  master_list <- master.creator(sequence_list, gc_results)
  
  for(sample in master_list){
    
    # Checking to see what loci are present for each sample based on PING_gc results
    locus_presence <- gc_results[,pmatch(sample, colnames(gc_results))]
    
    # As long as the PING_gc locus result isn't 0, it will be run, this includes discrepencies
    loci.list <- rownames(gc_results)[locus_presence != 0]
    
    # Cutting down the loci.list to only include supported loci
    loci.list <- loci.list[loci.list %in% supported.loci]
    
    # Mutating 2DS5 to 2DS35 for allele calling
    if("2DS5" %in% loci.list){
      loci.list[loci.list == "2DS5"] <- "2DS35"
    }
    
    if("2DL23" %in% loci.list){
      loci.list <- c(loci.list, "2DL2")
    }
    
    # Going backwards from master_list sample name to sequence_list sample name
    sample <- sequence_list[pmatch(sample, sequence_list)]
    
    for(current.locus in loci.list){
      
      if(current.locus == "2DS3" && any("2DS35" %in% loci.list)){
        next
      }
      
      # Setting up for genotype caller
      cat(paste("\n\n---> Starting genotype calling for", sample, "at", current.locus, "<---\n\n"))
      
      # Get sample copy number information for current.locus
      copy_number <- gc_results[pmatch(current.locus, rownames(gc_results)),pmatch(sample, colnames(gc_results))]
      
      if("2DS35" == current.locus){
        current.locus <- "2DS5"
        copy_number <- gc_results[pmatch(current.locus, rownames(gc_results)),pmatch(sample, colnames(gc_results))]
        current.locus <- "2DS35"
      }
      
      # Generating known alleles for the current locus
      alleles_gen.vcf(current.locus)
      
      # Generating position information for the current locus
      SOS_locus_lookup <- snps_gen.vcf(current.locus)
      
      positions <- gsub("X", "", names(SOS_locus_lookup)[-1])
      
      
      # Running the bowtie2 scripts
      ping.locus(sample, current.locus, is_gz)
      
      
      # 2DL23 specific logic
#      if("2DL23" == current.locus){
        # ## Determine if 2DL2+ 2DL3+
        # # 2DL2+ test
        # pos2DL2 <- FALSE
        # vcf_table <- tryCatch(read.table(paste0(results.directory, "Vcf/", sample, "_2DL2nuc.vcf"), header = F), error=function(e) NULL)
        # 
        # if(!is.null(vcf_table)){
        #   vcf_table <- nodels.vcf(vcf_table)
        # }
        # 
        # if(length(vcf_table[,1]) == 383) {
        #   pos2DL2 <- TRUE
        #   alleles_gen.vcf("2DL2")
        # }
        # 
        # # 2DL3+ test
        # pos2DL3 <- FALSE
        # vcf_table <- tryCatch(read.table(paste0(results.directory, "Vcf/", sample, "_2DL3nuc.vcf"), header = F), error=function(e) NULL)
        # 
        # if(!is.null(vcf_table)){
        #   vcf_table <- nodels.vcf(vcf_table)
        # }
        # 
        # if(length(vcf_table[,1]) == 386) {
        #   pos2DL3 <- TRUE
        # }
        # 
        # 
        # presence2DL2or3 <- list("two" = pos2DL2, "three" = pos2DL3)
      #}
      
      # Running the allele calling scripts
      has_genotype <- ping.caller(sample, current.locus)
      
      if(is.null(has_genotype)){
        cat(paste0("\nNo genotype or new alleles found for ", sample, " at ", current.locus, "\n"))
        write.table(sample, paste0(results.directory, "KIRcaller/bad_", current.locus, "_files.txt"), append = TRUE, row.names = F, col.names = F, quote = F)
      }
      
      cat(paste("\n---> Finished genotype calling for", sample, "at", current.locus, "<---\n"))
    }
    
    if(is_gz){
      file.remove(paste0(sample, no.gz.pattern.1))
      file.remove(paste0(sample, no.gz.pattern.2))
    }else{
      file.remove(paste0(sample, fastq.pattern.1))
      file.remove(paste0(sample, fastq.pattern.2))
    }
    
    cat(paste("\n\n\n*****---------- Finished with genotype calling for", sample, "----------*****\n\n"))
  }
  
  cat("\n\nfread warnings are usually due to sequence files being too small, please see Error.log if this shows up.")
  cat("\n\nAll finished! Please look for the final genotype calls in the Results/KIRcaller folder.")
}
