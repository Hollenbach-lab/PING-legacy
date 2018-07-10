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



###################
# CALLER FUNCTIONS -------------------------------------------------
###################

##==================================================================
## HAPLO caller functions ------------------------------------------
##==================================================================

## Support funciton for make_bed, this function turns a msf file into a dataframe of aligned alleles
msf_to_allele_frame <- function(haplo_directory, current_locus){
  ## Converts MSF files from IPD-KIR to a table of allele coding sequences.
  ## This is used for allele calling by comparing variant positions. We want
  ## to update this code to be able to handle copy and pasted MSF output.
  
  if(current_locus == '2DL5A' | current_locus == '2DL5B'){
    old_locus <- current_locus
    current_locus <- '2DL5'
  }
  ## Read in MSF file for current_locus in msf_directory
  msf_file_path <- file.path(haplo_directory,paste0('KIR',current_locus,'_nuc.msf'))
  msf_file <- read.table(msf_file_path, sep='\n', stringsAsFactors = F)
  
  ## Look at the second line, pull out the number of bases following 'MSF: '
  number_of_bases_raw <- unlist(strsplit(msf_file[2,], ' '))
  number_of_bases <- number_of_bases_raw[grep('MSF:', number_of_bases_raw)+1]
  number_of_bases <- as.integer(number_of_bases)
  
  ## Grab current_locus allele names from the header of the MSF file
  allele_names_raw <- msf_file[grep('Name:', msf_file[,1]),]
  allele_names_raw <- tstrsplit(allele_names_raw, 'Name: ')[[2]]
  allele_names <- tstrsplit(allele_names_raw, ' ')[[1]]
  
  ## Initialize data.frame for all alleles and all positions
  allele_frame <- data.frame(matrix(nrow=length(allele_names), ncol=number_of_bases), row.names = allele_names)
  
  ## Cut off the header from the MSF file
  msf_file <- msf_file[grep('//', msf_file[,1])+1:nrow(msf_file),]
  
  for(allele in allele_names){
    allele_lines_raw <- msf_file[grep(paste0(allele, ' '), msf_file, fixed=T)]
    allele_lines <- tstrsplit(allele_lines_raw, allele, fixed=T)[[2]]
    allele_string <- paste0(unlist(strsplit(allele_lines, ' ')), collapse='')
    allele_vector <- strsplit(allele_string, '')[[1]]
    allele_frame[allele,] <- allele_vector
  }
  return(allele_frame)
}

## Returns a list that converts between bed positions and allele positions
make_bed_to_pos_conv <- function(bed_file_path){
  ## This function should work to convert between bed genomic coordinates 
  ## and CDS coordinates. It outputs this conversion as a list, we should
  ## change this to data.frame format. This function is meant to handle
  ## multiple alleles at once, but should work for single alleles just fine.
  ## The allele naming might be weird for normal PING bed files
  
  bed_file_frame <- data.frame(read.table(bed_file_path, header=F, stringsAsFactors = F))
  bed_allele_list <- unique(bed_file_frame[,'V1'])
  list_bed_to_pos <- list()
  for(bed_allele in bed_allele_list){
    list_bed_to_pos[[bed_allele]] <- c()
    bed_allele_frame <- bed_file_frame[bed_file_frame[,'V1'] == bed_allele,]
    
    for(i in 1:nrow(bed_allele_frame)){
      first_value <- bed_allele_frame[i,2]+1
      last_value <- bed_allele_frame[i,3]
      iter_values <- iterate(first_value, last_value)
      list_bed_to_pos[[bed_allele]] <- c(list_bed_to_pos[[bed_allele]], iter_values)
    }
  }
  return(list_bed_to_pos)
}

## Gotta iterate
iterate <- function(first_value, last_value){
  ## Support function for make_bed_to_pos_conv
  return(as.integer(first_value):as.integer(last_value))
}

## Returns boolean
is_nuc <- function(chr){
  nuc_list <- c('A', 'T', 'C', 'G')
  return(as.character(chr) %in% nuc_list)
}

## Counts how many unique nucleotides are in a vector
num_unique_nuc <- function(vector_to_check){
  return(sum(is_nuc(names(table(vector_to_check)))))
}


##==================================================================
## General Functions -----------------------------------------------
##==================================================================

# Creates results directory
results_directory <- function(results.directory) {
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


# Allele generation for locus (relies on resource files)
alleles_gen.vcf <- function(current.locus,resource.location){
  all_alleles_preKFF <- read.delim(paste0(resource.location,"/All_", current.locus, "_preKFF.fas"), header=F, colClasses="character")
  length_seq <- nchar(all_alleles_preKFF[1,2])
  num_seq <- length(all_alleles_preKFF[,1])
  top_row <- data.frame(cbind(num_seq, length_seq))
  names(top_row) <- names(all_alleles_preKFF)
  all_alleles <- rbind(top_row, all_alleles_preKFF)
  write.table(all_alleles, paste0(results.directory, "KIRcaller/All_", current.locus, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Snps generation for locus (relies on resource files) --> SOS lookup table
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
  #write.table(snps_locus, paste0(results.directory, "KIRcaller/KIR_", current.locus, "_alleles.txt"), quote = FALSE, row.names = FALSE)
  
  SOS_locus_lookup <- snps_locus
  
  ### Remove positions that are only variable due to unknown/deletion sequence from SOS lookup table
  col_names_SOS_lookup <- names(SOS_locus_lookup)
  pos_to_exclude <- c()
  for(i in 2:length(col_names_SOS_lookup)){
    if (num_unique_nuc(SOS_locus_lookup[,i]) == 1){
      pos_to_exclude <- c(pos_to_exclude,col_names_SOS_lookup[i])
    }
  }
  # Remove these postions
  SOS_locus_lookup[,pos_to_exclude] <- NULL
  
  # Output Variable Positions prior to conversion to genomic positions
  write.table(SOS_locus_lookup, paste0(results.directory, "KIRcaller/KIR_", current.locus, "_alleles.txt"), quote = FALSE, row.names = FALSE)
  
  # Convert CDS positions to genomic positions
  KIR_locus_gene_VScDNA <- read.delim(paste0("Resources/caller_resources/KIR", current.locus, "geneVScDNA.txt"))
  cdna <- paste0("X", KIR_locus_gene_VScDNA$cDNA)
  gene <- paste0("X", KIR_locus_gene_VScDNA$gene)
  map_pos <- data.frame(cbind(cdna, gene))
  pos_trans <- as.character(map_pos[match(names(SOS_locus_lookup[-1]), map_pos$cdna), 2])
  names(SOS_locus_lookup) <- c("allele", pos_trans)
  
  return(SOS_locus_lookup)
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



### NESTED GENERAL FUNCTIONS

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
  
  if("2DL1" == current.locus)  {KIR_2DL1(sample, is_gz)}
  
  if("2DL23" == current.locus) {KIR_2DL23(sample, is_gz)}
  
  if("2DL4" == current.locus)  {KIR_2DL4(sample, is_gz)}
  
  if("2DL5" == current.locus)  {KIR_2DL5(sample, is_gz)}
  
  if("2DP1" == current.locus)  {KIR_2DP1(sample, is_gz)}
  
  if("2DS3" == current.locus && !any("2DS35" %in% loci.list)){
    KIR_2DS3(sample, is_gz)
  }else if("2DS35" == current.locus){
    KIR_2DS35(sample, is_gz)
  }else if("2DS3" == current.locus){
    
  }
  
  if("2DS4" == current.locus) {KIR_2DS4(sample, is_gz)}
  
  if("3DL1" == current.locus && any("3DS1" %in% loci.list)){
    KIR_3DL1S1(sample, is_gz)
  }else if("3DL1" == current.locus){
    KIR_3DL1(sample, is_gz)
  }else if("3DS1" == current.locus && !any("3DL1" %in% loci.list)){
    KIR_3DS1(sample, is_gz)
  }
  
  if("3DL2" == current.locus) {KIR_3DL2(sample, is_gz)}
  
  if("3DL3" == current.locus) {KIR_3DL3(sample, is_gz)}
}


##==================================================================
## Remaining PING Functions ----------------------------------------
##==================================================================

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

# Special 2DL23 genotyping
# THIS SHOULD BE MOVED TO POST PROCESSING!!
genos_2DL23 <- function(genos.out){
  
  vector_call_2DL2 <- NA
  vector_call_2DL3 <- NA
  
  for(item in genos.out){
    
    item_frame <- as.data.frame(item)
    
    # Grabbing and reshaping each type of call into a vector
    call_2DL2 <- grep("KIR2DL2_", item_frame)
    if(length(call_2DL2) > 0){
      frame_call_2DL2 <- item_frame[call_2DL2]
      
      vector_call_2DL2 <- c(vector_call_2DL2, t(frame_call_2DL2))
      vector_call_2DL2 <- vector_call_2DL2[!is.na(vector_call_2DL2)]
    }
    
    call_2DL3 <- grep("KIR2DL3_", item_frame)
    if(length(call_2DL3) > 0){
      frame_call_2DL3 <- item_frame[call_2DL3]
      
      vector_call_2DL3 <- c(vector_call_2DL3, t(frame_call_2DL3))
      vector_call_2DL3 <- vector_call_2DL3[!is.na(vector_call_2DL3)]
    }
  }
  
  
  # Finding any 2DL2 alleles by finding the intersection of each 2DL2 call
  if(length(vector_call_2DL2) == 3){
    one <- unlist(strsplit(vector_call_2DL2, "/")[1])
    two <- unlist(strsplit(vector_call_2DL2, "/")[2])
    three <- unlist(strsplit(vector_call_2DL2, "/")[3])
    
    step_one <- intersect(two, three)
    allele_one <- intersect(one, step_one)
    
  }else if(length(vector_call_2DL2) == 4){
    one <- unlist(strsplit(vector_call_2DL2, "/")[1])
    two <- unlist(strsplit(vector_call_2DL2, "/")[2])
    three <- unlist(strsplit(vector_call_2DL2, "/")[3])
    four <- unlist(strsplit(vector_call_2DL2, "/")[4])
    
    allele_one <- intersect(one, three)
    allele_two <- intersect(two, four)
    
    if(length(allele_one) == 0 || length(allele_two) == 0){
      allele_one <- intersect(one, four)
      allele_two <- intersect(two, three)
    }
  }
  
  # Finding any 2DL3 alleles by finding the intersection of each 2DL3 call
  if(length(vector_call_2DL3) == 3){
    one <- unlist(strsplit(vector_call_2DL3, "/")[1])
    two <- unlist(strsplit(vector_call_2DL3, "/")[2])
    three <- unlist(strsplit(vector_call_2DL3, "/")[3])
    
    step_one <- intersect(two, three)
    allele_two <- intersect(one, step_one)
    
  }else if(length(vector_call_2DL3) == 4){
    one <- unlist(strsplit(vector_call_2DL3, "/")[1])
    two <- unlist(strsplit(vector_call_2DL3, "/")[2])
    three <- unlist(strsplit(vector_call_2DL3, "/")[3])
    four <- unlist(strsplit(vector_call_2DL3, "/")[4])
    
    allele_one <- intersect(one, three)
    allele_two <- intersect(two, four)
    
    if(length(allele_one) == 0 || length(allele_two) == 0){
      allele_one <- intersect(one, four)
      allele_two <- intersect(two, three)
    }
  }
  
  
  ## BUG FIX 7/18 new genotypes were not being caught by the tryCatch because they were character(0), so I am removing objects
  if(exists("allele_one") && length(allele_one) == 0){
    remove(allele_one)
  }
  if(exists("allele_two") && length(allele_two) == 0){
    remove(allele_two)
  }
  
  genos.out <- genos.out[1]
  genos.out[[1]]$X1 <- tryCatch(paste(allele_one, collapse = "/"), error=function(e) "new")
  genos.out[[1]]$X2 <- tryCatch(paste(allele_two, collapse = "/"), error=function(e) "new")
  
  return(genos.out)
}


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
  
  kff_legend <- read.delim(paste0("Resources/caller_resources/KIR_", current.locus, "_kff_legend.txt"), sep = "")
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
    lookitup[grep(kff_legend[i, "allele"], lookitup$allele), "string"] <- paste0(lookitup[grep(kff_legend[i, "allele"], lookitup$allele), "string"], as.character(kff_legend[i, "code"]))
  }
  
  
  ## Running allele calling again with the new codes
  allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 3]),"new")
  allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 3]),"new")
  poss_genos <- data.frame(cbind(as.character(allele1),as.character(allele2)), stringsAsFactors=FALSE)
  
  return(poss_genos)
}



##==================================================================
## PING New Functions ----------------------------------------------
##==================================================================

## INITIALIZE locus dataframes for new allele variants
initialize_new_alleles_df <- function(resource.location){
  # Only works if resource files have this naming scheme: All_[LOCUS]_preKFF.fas
  fasta.files     <- dir(resource.location, pattern = "preKFF.fas")
  locus_list      <- unlist(lapply(strsplit(fasta.files,"_"),'[[',2))
  
  # Loci not supported by PING: 2DS1,2DS2,3DP1
  ### DEBUG THIS!! 2DL23 currently not working when inputted into snps_gen.vcf, it is because not all alleles are equal in length
  loci_to_exclude <- c("2DS1","2DS2","3DP1","2DL2","2DL3","2DS3","2DS5","2DL23")
  locus_list      <- locus_list[! locus_list %in% loci_to_exclude]
  
  # This list will contain new alleles data frames across all KIR loci
  newAllelesDataFrames.list <- list()
  
  for(i in 1:length(locus_list)){
    locus              <- locus_list[i]
    locus_alleles_file <- fasta.files[grepl(paste0(locus,"_"),fasta.files)]
    locus_alleles_fh   <- paste0(resource.location,"/",locus_alleles_file)
    
    # General SOS lookup table for each locus
    # - This will make it so the 'KIRcaller' results directory will always have allele FASTAs for all loci, not just what is being run 
    alleles_gen.vcf(locus, resource.location)
    SOS_locus_lookup   <- snps_gen.vcf(locus)
    
    newVariantsDataFrame <- data.frame(matrix(nrow=0, ncol=length(names(SOS_locus_lookup))))
    names(newVariantsDataFrame) <- names(SOS_locus_lookup)
    newVariantsDataFrame$allele <- NULL
    
    # Add initializes Variant Dataframe into a nested list
    #newAllelesDataFrames.list[[i]] <- list()
    #newAllelesDataFrames.list[[i]][[1]] <- locus
    #newAllelesDataFrames.list[[i]][[2]] <- newVariantsDataFrame
    newAllelesDataFrames.list[[locus]] <- newVariantsDataFrame
  }
  
  return(newAllelesDataFrames.list)
}


## Formatting name for New allele
format_new_allele_name <- function(defining_new_allele_pos,Aref){
  new_allele_name_list <- c(Aref,"_new")
  new_allele_pos_vec   <- names(defining_new_allele_pos)
  for(pos in new_allele_pos_vec){
    nuc <- defining_new_allele_pos[,pos]
    conv_pos <- gsub("X","",pos)
    new_allele_name_list <- c(new_allele_name_list,conv_pos,nuc,".")
  }
  new_allele_name <- paste(new_allele_name_list,collapse = "")
  new_allele_name <- substr(new_allele_name,1,nchar(new_allele_name)-1) # remove last char "."
  
  return(new_allele_name)
}


## Create Genotype string for New alleles Results
create_new_alleles_genos <- function(new_genos){
  new_genos$string1 <- NULL
  new_genos$string2 <- NULL
  
  # Collapse Genotype combinations that overlap (Ex: same first allele, different second allele)
  new_genos_list <- list()
  for(i in 1:length(row.names(new_genos))){
    A1 <- new_genos[i,1]
    A2 <- new_genos[i,2]
    
    geno_combo <- c()
    if (grepl("_new",A1) & grepl("_new",A2)){
      geno_combo <- sort(new_genos[i,]) 
    } else if (grepl("_new",A1)){
      geno_combo <- c(A1,A2)
    } else if (grepl("_new",A2)){
      geno_combo <- c(A2,A1)
    }
    
    A1 <- geno_combo[1]
    A2 <- geno_combo[2]
    new_genos_list[[A1]] <- c(new_genos_list[[A1]],A2,"/")
    
  }
  
  # Create Final genotype string
  final_new_genos_list  <- c()
  col_names_new_alleles <- names(new_genos_list)
  for (i in 1:length(names(new_genos_list))){
    new_allele <- col_names_new_alleles[i]
    A2_vec <- unlist(new_genos_list[[i]])
    A2_vec <- A2_vec[1:length(A2_vec)-1]  # get rid of the "/" that will always be on the end
    A2 <- paste(A2_vec,collapse = "")
    new_geno <- paste0(new_allele,"+",A2)
    final_new_genos_list <- c(final_new_genos_list,new_geno)
  }
  
  new_geno_final <- paste(final_new_genos_list,collapse = " ")
  return(new_geno_final)
}



## ALTERED NEW ALLELES FUNCTION: creates genotype string from new allele calls, each new allele will have a unique ID
##  - mismatchThresh = how many mismatches are allowed between a new allele and a known allele to be considered for genos output
##                     new alleles that fail to meet the threshold with be assigned 'new' in genos output
new_alleles_genos.vcf <- function(x, current.locus,DPthresh=6, mismatchThresh = 1){
  KIR_sample_snps <- emformat.vcf(x,DPthresh)
  
  if(is.null(KIR_sample_snps)){
    return(NULL)
  }
  
  # CDS Read Positions (prior to conversion to genomic positions) 
  SOS_locus_lookup_cds <- read.delim(paste0(results.directory,"KIRcaller/KIR_",current.locus,"_alleles.txt"), sep = " ")
  reads_pos_cds        <- names(SOS_locus_lookup_cds)
  reads_pos_cds        <- reads_pos_cds[grepl("X",reads_pos_cds)]
  reads_pos_genomic    <- names(SOS_locus_lookup)
  reads_pos_genomic    <- reads_pos_genomic[grepl("X",reads_pos_genomic)]
  
  
  # Genomic Read Positions
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
  genos <- data.frame(cbind(as.character(allele1), as.character(allele2)),stringsAsFactors = F)
  
  ## EXACTLY THE SAME AS KIR.CALL UNTIL THIS POINT ^^^^
  
  #genos$X1X2 <- apply(genos, 1, paste, collapse="")
  genos$string1 <- poss_genos$X1
  genos$string2 <- poss_genos$X2
  #genos <- genos[!genos$X1X2=="newnew",]  # new+new should not be removed
  
  # This is in cases where there is new+new
  #if(is.data.frame(genos) && nrow(genos)==0) {return(NULL)}
  
  #### CREATE NEW ALLELE DATA FRAME HERE!!!  ################################################
  pos_vec <- reads_pos
  # Conversion from Genomic to CDS positions
  pos_vec_conv        <- reads_pos_cds
  names(pos_vec_conv) <- reads_pos_genomic
  
  new_alleles_vec      <- c()
  new_alleles_seq_vec  <- c()
  new_alleles_geno_vec <- c()
  new_alleles_categ    <- c() # possible genotypes will always be either new+known(allele) or new+new
  # Create Data Frame of all variable positions for new allele variants
  new_alleles_pos_var_df <- data.frame(matrix(nrow=0, ncol=length(pos_vec_conv)))
  names(new_alleles_pos_var_df) <- pos_vec_conv
  
  for(i in 1:length(row.names(genos))){
    A1                <- genos[i,1]
    A2                <- genos[i,2]
    A1_seq            <- genos[i,3]
    A2_seq            <- genos[i,4]
    A1_seq_vec        <- unlist(strsplit(A1_seq,""))
    A2_seq_vec        <- unlist(strsplit(A2_seq,""))
    names(A1_seq_vec) <- pos_vec
    names(A2_seq_vec) <- pos_vec
    
    Aref <- ref_allele_names[grep(current.locus,ref_allele_names)]
    Aref_seq_vec <- SOS_locus_lookup_reads[SOS_locus_lookup_reads$allele == Aref,]
    Aref_seq_vec$allele <- NULL
    
    # Category of possible genotype
    categ_geno <- "new+known"
    if (A1 =="new" && A2 == "new"){categ_geno <- "new+new"}
    
    if (categ_geno == "new+known"){
      # Only one new allele detected in this possible genotype
      if (A1 == "new"){
        comp_Aref_new_df <- rbind(Aref_seq_vec,A1_seq_vec)
        unique_nuc_vec   <- apply(comp_Aref_new_df,2, num_unique_nuc)
        unique_nuc_vec   <- unique_nuc_vec[unique_nuc_vec > 1]
        
        # Add reformatted name back into the genos data frame to replace the 'new' flags
        defining_new_allele_pos            <- as.data.frame(comp_Aref_new_df[2,names(unique_nuc_vec)])
        names(defining_new_allele_pos)     <- names(unique_nuc_vec)
        row.names(defining_new_allele_pos) <- "2"
        # Convert Genomic positions to CDS positions
        names(defining_new_allele_pos) <- pos_vec_conv[names(defining_new_allele_pos)]
        # Generate new allele name and substitute it to replace 'new' flags
        new_allele_name                <- format_new_allele_name(defining_new_allele_pos,Aref)
        genos[i,1]                     <- new_allele_name
        poss_genotype                  <- c(genos[i,1],genos[i,2])
        poss_genotype                  <- paste(poss_genotype,collapse = "+")
        
        # Store new allele variant in a vector
        new_alleles_vec      <- c(new_alleles_vec,new_allele_name)
        new_alleles_geno_vec <- c(new_alleles_geno_vec,poss_genotype)
        new_alleles_seq_vec  <- c(new_alleles_seq_vec, A1_seq)
        new_alleles_categ    <- c(new_alleles_categ,categ_geno)
        
        # Format Data Frame entry
        #   - If there are missing variable positions that have been filtered due to low depth, subsitute in an N
        new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec_conv)))
        names(new_allele_entry)              <- reads_pos_genomic
        new_allele_entry[1,]                 <- "N"
        new_allele_entry[,names(A1_seq_vec)] <- A1_seq_vec
        row.names(new_allele_entry)          <- new_allele_name
        names(new_allele_entry)              <- as.character(pos_vec_conv)
        new_alleles_pos_var_df               <- rbind(new_alleles_pos_var_df, new_allele_entry)
        
      }
      
      if (A2 == "new"){
        comp_Aref_new_df <- rbind(Aref_seq_vec,A2_seq_vec)
        unique_nuc_vec   <- apply(comp_Aref_new_df,2, num_unique_nuc)
        unique_nuc_vec   <- unique_nuc_vec[unique_nuc_vec > 1]
        
        # Add reformatted name back into the genos data frame to replace the 'new' flags
        defining_new_allele_pos            <- as.data.frame(comp_Aref_new_df[2,names(unique_nuc_vec)])
        names(defining_new_allele_pos)     <- names(unique_nuc_vec)
        row.names(defining_new_allele_pos) <- "2"
        # Convert Genomic positions to CDS positions
        names(defining_new_allele_pos) <- pos_vec_conv[names(defining_new_allele_pos)]
        # Generate new allele name and substitute it to replace 'new' flags
        new_allele_name                <- format_new_allele_name(defining_new_allele_pos,Aref)
        genos[i,2]                     <- new_allele_name
        poss_genotype                  <- c(genos[i,2],genos[i,1])
        poss_genotype                  <- paste(poss_genotype,collapse = "+")
        
        # Store new allele variant in a vector
        new_alleles_vec      <- c(new_alleles_vec,new_allele_name)
        new_alleles_geno_vec <- c(new_alleles_geno_vec,poss_genotype)
        new_alleles_seq_vec  <- c(new_alleles_seq_vec, A2_seq)
        new_alleles_categ    <- c(new_alleles_categ,categ_geno)
        
        # Format Data Frame entry
        #   - If there are missing variable positions that have been filtered due to low depth, subsitute in an N
        new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec_conv)))
        names(new_allele_entry)              <- reads_pos_genomic
        new_allele_entry[1,]                 <- "N"
        new_allele_entry[,names(A2_seq_vec)] <- A2_seq_vec
        row.names(new_allele_entry)          <- new_allele_name
        names(new_allele_entry)              <- as.character(pos_vec_conv)
        new_alleles_pos_var_df               <- rbind(new_alleles_pos_var_df, new_allele_entry)
      }
    } else{
      # No known allele detected in this possible genotype
      a1_comp_Aref_new_df <- rbind(Aref_seq_vec,A1_seq_vec)
      a1_unique_nuc_vec   <- apply(a1_comp_Aref_new_df,2, num_unique_nuc)
      a1_unique_nuc_vec   <- a1_unique_nuc_vec[a1_unique_nuc_vec > 1]
      
      a2_comp_Aref_new_df <- rbind(Aref_seq_vec,A2_seq_vec)
      a2_unique_nuc_vec   <- apply(a2_comp_Aref_new_df,2, num_unique_nuc)
      a2_unique_nuc_vec   <- a2_unique_nuc_vec[a2_unique_nuc_vec > 1]
      
      
      # Add reformatted name back into the genos data frame to replace the 'new' flags
      a1_defining_new_allele_pos            <- as.data.frame(a1_comp_Aref_new_df[2,names(a1_unique_nuc_vec)])
      names(a1_defining_new_allele_pos)     <- names(a1_unique_nuc_vec)
      row.names(a1_defining_new_allele_pos) <- "2"
      a2_defining_new_allele_pos            <- as.data.frame(a2_comp_Aref_new_df[2,names(a2_unique_nuc_vec)])
      names(a2_defining_new_allele_pos)     <- names(a2_unique_nuc_vec)
      row.names(a2_defining_new_allele_pos) <- "2"
      
      # Convert Genomic positions to CDS positions
      names(a1_defining_new_allele_pos) <- pos_vec_conv[names(a1_defining_new_allele_pos)]
      names(a2_defining_new_allele_pos) <- pos_vec_conv[names(a2_defining_new_allele_pos)]
      # Generate new allele name and substitute it to replace 'new' flags
      a1_new_allele_name                <- format_new_allele_name(a1_defining_new_allele_pos,Aref)
      a2_new_allele_name                <- format_new_allele_name(a2_defining_new_allele_pos,Aref)
      genos[i,1]                        <- a1_new_allele_name
      genos[i,2]                        <- a2_new_allele_name
      
      # CHANGE this part: add conditions to how possible genotypes are formed for new+new
      # new allele variants are the same
      if (a1_new_allele_name == a2_new_allele_name && A1_seq == A2_seq){
        poss_genotype                  <- c(genos[i,1],genos[i,2])
        poss_genotype                  <- paste(poss_genotype,collapse = "+")
        
        # Store new allele variant in a vector
        new_alleles_vec      <- c(new_alleles_vec,a1_new_allele_name)
        new_alleles_geno_vec <- c(new_alleles_geno_vec,poss_genotype)
        new_alleles_seq_vec  <- c(new_alleles_seq_vec, A1_seq)
        new_alleles_categ    <- c(new_alleles_categ,categ_geno)
        
        # Format Data Frame entry
        #   - If there are missing variable positions that have been filtered due to low depth, subsitute in an N
        new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec_conv)))
        names(new_allele_entry)              <- reads_pos_genomic
        new_allele_entry[1,]                 <- "N"
        new_allele_entry[,names(A1_seq_vec)] <- A1_seq_vec
        row.names(new_allele_entry)          <- a1_new_allele_name
        names(new_allele_entry)              <- as.character(pos_vec_conv)
        new_alleles_pos_var_df               <- rbind(new_alleles_pos_var_df, new_allele_entry)
        
      } else{
        # Create new alleles entry for each distinct variant
        a1_poss_genotype                  <- c(genos[i,1],genos[i,2])
        a1_poss_genotype                  <- paste(a1_poss_genotype,collapse = "+")
        
        a2_poss_genotype                  <- c(genos[i,2],genos[i,1])
        a2_poss_genotype                  <- paste(a2_poss_genotype,collapse = "+")
        
        # Store new allele variant in a vector
        new_alleles_vec      <- c(new_alleles_vec,a1_new_allele_name)
        new_alleles_geno_vec <- c(new_alleles_geno_vec,a1_poss_genotype)
        new_alleles_seq_vec  <- c(new_alleles_seq_vec, A1_seq)
        new_alleles_categ    <- c(new_alleles_categ,categ_geno)
        
        new_alleles_vec      <- c(new_alleles_vec,a2_new_allele_name)
        new_alleles_geno_vec <- c(new_alleles_geno_vec,a2_poss_genotype)
        new_alleles_seq_vec  <- c(new_alleles_seq_vec, A2_seq)
        new_alleles_categ    <- c(new_alleles_categ,categ_geno)
        
        
        # Format Data Frame entry
        a1_new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec_conv)))
        names(a1_new_allele_entry)              <- reads_pos_genomic
        a1_new_allele_entry[1,]                 <- "N"
        a1_new_allele_entry[,names(A1_seq_vec)] <- A1_seq_vec
        row.names(a1_new_allele_entry)          <- a1_new_allele_name
        names(a1_new_allele_entry)              <- as.character(pos_vec_conv)
        new_alleles_pos_var_df                  <- rbind(new_alleles_pos_var_df, a1_new_allele_entry)
        
        a2_new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec_conv)))
        names(a2_new_allele_entry)              <- reads_pos_genomic
        a2_new_allele_entry[1,]                 <- "N"
        a2_new_allele_entry[,names(A2_seq_vec)] <- A2_seq_vec
        row.names(a2_new_allele_entry)          <- a2_new_allele_name
        names(a2_new_allele_entry)              <- as.character(pos_vec_conv)
        new_alleles_pos_var_df               <- rbind(new_alleles_pos_var_df, a2_new_allele_entry)
      }
    }
    
  }
  ####
  
  names(new_alleles_seq_vec)  <- new_alleles_vec
  names(new_alleles_geno_vec) <- new_alleles_vec
  names(new_alleles_categ)    <- new_alleles_vec
  
  
  # Compute difference (mismatches) of new allele compared to known alleles for the KIR locus specified  
  res.df <- data.frame(adist(new_alleles_seq_vec,lookitup$string))
  names(res.df) <- lookitup$allele
  
  # Start formatting Data frame of new variants
  new_alleles_df      <- as.data.frame(new_alleles_seq_vec)
  new_alleles_df$geno <- new_alleles_geno_vec
  
  # NO filter: Create genotype string from genos data frame  --> not being used
  #new_alleles.geno <- create_new_alleles_genos(genos)
  
  # FILTER: only keep new alleles combos that are only one mismatch away from a known allele 
  oneoff.alleles.df            <- data.frame(matrix(nrow = length(row.names(res.df)),ncol = 1))
  row.names(oneoff.alleles.df) <- row.names(res.df)
  names(oneoff.alleles.df)     <- "new.allele.near"
  
  oneoff.coords  <- data.frame(which(res.df==mismatchThresh,arr.ind=TRUE))
  oneoff.rows    <- unique(oneoff.coords$row)
  oneoff.df      <- res.df[oneoff.rows,]
  
  known_alleles.vec <- names(res.df)
  variants.vec      <- row.names(res.df)
  for (r in oneoff.rows){
    coord_subset    <- subset(oneoff.coords,oneoff.coords$row == r)
    col_of_interest <- coord_subset$col
    variant_name    <- variants.vec[r]
    oneoff_alleles  <- known_alleles.vec[col_of_interest] 
    
    all_nearest_alleles <- ""
    if (length(oneoff_alleles) == 1){
      all_nearest_alleles <- oneoff_alleles
    } else{
      all_nearest_alleles <- paste(oneoff_alleles,collapse = "/")
    }
    oneoff.alleles.df[variant_name,] <- all_nearest_alleles
  }
  
  # Remove NAs for One off alleles data frame
  oneoff.alleles.df <- subset(oneoff.alleles.df, !is.na(oneoff.alleles.df$new.allele.near))

  # If any new allele is near a known allele, include that info
  new_alleles_df$new.allele.near <- NA
  if (length(row.names(oneoff.alleles.df)) != 0){
    for (target_allele in row.names(oneoff.alleles.df)){
      curr_geno           <- unlist(strsplit(as.character(new_alleles_df[target_allele,2]),"+",fixed = T))
      nearest_allele      <- oneoff.alleles.df[target_allele,]
      curr_geno[1]        <- nearest_allele
      curr_geno           <- paste(curr_geno,collapse = "+")
      new_alleles_df[target_allele,2] <- curr_geno
      new_alleles_df[target_allele,3] <- nearest_allele 
    }
  }
  names(new_alleles_df) <- c("variable.pos.bases","geno","new.allele.near")
  new_alleles_df$categ  <- new_alleles_categ
  new_alleles_df$sample <- sample
  new_alleles_df$new.allele.name <- row.names(new_alleles_df)
  
  # Re-arrange into proper output form 
  #new_alleles_df$new.allele.near <- as.character(new_alleles_df$new.allele.near)
  #new_alleles_df$new.allele.near <- ifelse(is.na(new_alleles_df$new.allele.near),'NA',new_alleles_df$new.allele.near)
  new_alleles_df <- as.data.frame(new_alleles_df[,c("sample","new.allele.name","new.allele.near","geno","categ","variable.pos.bases")])
  
  # Merge new allele genotype info with SNP position info
  new_alleles_df.merged <- merge(new_alleles_df,new_alleles_pos_var_df, by=0,all = T)
  new_alleles_df.merged$variable.pos.bases <- NULL
  new_alleles_df.merged$Row.names          <- NULL
  
  new_alleles_df.merged
}



##==================================================================
## PING Allele Caller Functions ------------------------------------
##==================================================================

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

## CHECKING VCF DATA FRAME - check.vcf, format.vcf, em.format

# Format VCF data frame for allele calling functions
format.vcf <- function(x, DPthresh = 6) {
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
  KIR_sample_snps <- KIR_sample_snps[!KIR_sample_snps$depth < DPthresh,]
  
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
emformat.vcf <- function(x,DPthresh = 6){
  KIR_sample_snps <- format.vcf(x, DPthresh)
  
  if(is.null(KIR_sample_snps)){
    return(NULL)
  }
  
  KIR_sample_snps$badcall <- KIR_sample_snps$var_called+KIR_sample_snps$fwrev
  # Removing DP4 threshold of 3
  #KIR_sample_snps <- KIR_sample_snps[!KIR_sample_snps$badcall == 1, ]
  
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
check.vcf <- function(x, DPthresh = 6){
  
  if(length(x[,1]) == 383 && current.locus == "2DL23"){
    # Generating position information for the current locus
    SOS_locus_lookup <<- snps_gen.vcf("2DL2")
    
    positions <<- gsub("X", "", names(SOS_locus_lookup)[-1])
  }
  
  KIR_sample_snps <- emformat.vcf(x,DPthresh)
  
  if(length(x[,1]) == 383 && current.locus == "2DL23"){
    # Generating position information for the current locus
    SOS_locus_lookup <<- snps_gen.vcf(current.locus)
    
    positions <<- gsub("X", "", names(SOS_locus_lookup)[-1])
  }
  
  good_vcf <- ifelse(dim(KIR_sample_snps)[1] > 2, "yes", "no")
  
  if(is.null(KIR_sample_snps)){
    good_vcf <- "no"
  }
  
  return(good_vcf)
}



## ENUMERATING ALL ALLELE COMBINATIONS

# ERROR Log function for haplo.enum() function
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



## ALLELE CALLING
#  - Call genotypes from VCF data frame
allele_call.vcf <- function(x, sample,DPthresh = 6){
  
  if(length(x[,1]) == 383 && current.locus == "2DL23"){
    # Generating position information for the current locus
    SOS_locus_lookup <<- snps_gen.vcf("2DL2")
    
    positions <<- gsub("X", "", names(SOS_locus_lookup)[-1])
  }
  
  KIR_sample_snps <- emformat.vcf(x,DPthresh)
  
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
    
    ###allele calling
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
    } else{
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
  #if( copy_number == 1 && length(poss_genos[,1]) !=0 && poss_genos$X1 == poss_genos$X2 && current.locus != "2DL23"){
  #  poss_genos$X2 <- paste0("KIR", current.locus, "_null")
  #}
  
  if(length(x[,1]) == 383 && current.locus == "2DL23"){
    # Generating position information for the current locus
    SOS_locus_lookup <<- snps_gen.vcf(current.locus)
    
    positions <<- gsub("X", "", names(SOS_locus_lookup)[-1])
  }
  
  genos <- poss_genos
  genos
}





## NEW ALLELES: Call possible new alleles from VCF data frame
new_alleles.vcf <- function(x, DPthresh=6){
  KIR_sample_snps <- emformat.vcf(x,DPthresh)
  
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



# Call possible new snps from VCF data frame
new_snps_call.vcf <- function(x, DPthresh = 6){
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
  # DP cutt off
  KIR_sample<-KIR_sample[!KIR_sample$depth < DPthresh,]
  
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


# Converts calls to Ucode
get_calls.vcf <- function(x, DPthresh = 6){
  KIR_sample_snps <- emformat.vcf(x, DPthresh)
  
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



## PRIMARY ALLELE CALLER FUNCTION
#  - Genotype caller (multi sample)
ping.caller <- function (sample, current.locus, DPthresh = 6) {
  
  # Format Sample ID for Genotype output
  sample_id <- unlist(strsplit(sample,"_"))[1]
  if(grepl("_", copy_number)){
    sample_id <- paste0(sample_id, "*")
  }
  
  
  # Read in vcf files and check if they are good ----------------------------
  
  vcf.location <- paste0(results.directory, "Vcf/")
  
  current_wd <- getwd()
  setwd(vcf.location)
  
  vcf_files <- list.files(pattern = paste0(sample, "_", current.locus, "nuc.vcf"))
  
  ## 2DL23 specific vcf_file finding
  if(current.locus == "2DL23"){
    if(presence2DL2or3$two){
      vcf_files <- c(vcf_files, list.files(pattern = paste0(sample, "_", "2DL2", "nuc.vcf")))
    }
    if(presence2DL2or3$three){
      vcf_files <- c(vcf_files, list.files(pattern = paste0(sample, "_", "2DL3", "nuc.vcf")))
    }
  }
  
  setwd(current_wd)
  
  # Return NULL if VCF files are empty (VCF files with header info only)
  for (i in vcf_files) {
    vcf_test <- tryCatch(read.table(paste0(vcf.location, i), header = F), error=function(e) NULL)
    if (is.null(vcf_test)){
      write.table(paste0(sample_id," empty_vcf"), paste0(results.directory, "KIRcaller/genos_", current.locus, ".txt"), append = TRUE, row.names=F, col.names=F, quote=F)
      return(NULL)
    }
  }
  
  vcf_list = lapply(paste0(vcf.location, vcf_files), read.table, header = FALSE)
  names(vcf_list) <- vcf_files
  
  vcf_list_indels <- lapply(vcf_list, indels.vcf)
  
  vcf_list <- lapply(vcf_list, nodels.vcf)
  
  check_vcf.out <- lapply(vcf_list, check.vcf,DPthresh = DPthresh)
  vcf_check.df <- data.frame(as.matrix(check_vcf.out))
  
  good <- data.frame(unlist(vcf_check.df))[,1] == "yes"
  vcf_files.good = vcf_files[good]
  
  # Insufficient Depth in VCF file for this sample
  if(length(vcf_files.good) == 0){
    write.table(paste0(sample_id," insufficient_depth"), paste0(results.directory, "KIRcaller/genos_", current.locus, ".txt"), append = TRUE, row.names=F, col.names=F, quote=F)
    return(NULL)
  }
  
  vcf_list.good = lapply(paste0(vcf.location, vcf_files.good), read.table, header=F)
  vcf_list.good <- lapply(vcf_list.good, nodels.vcf)
  names(vcf_list.good) <- vcf_files.good
  
  vcf_files.bad <- data.frame(vcf_files[!good])
  
  write.table(vcf_files.bad, paste0(results.directory, "KIRcaller/bad_", current.locus, "_files.txt"), append = TRUE, row.names = F, col.names = F, quote = F)
  
  indel_presence <- any(lapply(vcf_list_indels, nrow) > 0)
  
  if(indel_presence) {
    if(length(vcf_list_indels) != 10){
      vcf_list_indels <- do.call(rbind, vcf_list_indels)
    }
    write.table(vcf_list_indels, paste0(results.directory, "KIRcaller/indels_", current.locus, ".txt"), append = TRUE, row.names = F, col.names = F)
  }
  
  # Run genotype caller -----------------------------------------------------
  
  genos.out <- lapply(vcf_list.good, allele_call.vcf, sample, DPthresh)
  
  if(all(unlist(lapply(genos.out, is.null))) || all(unlist(lapply(genos.out, nrow)) == 0)){
    cat("\nNo genotype found, moving on to look for new alleles.\n")
    genos.df <- as.data.frame(do.call(rbind, genos.out))
    genos.df <- genos.df[!genos.df$X1 == "new",]
    genos.df <- genos.df[!genos.df$X2 == "new",]
    genos.df$X1 <- NULL
    genos.df$X2 <- NULL
  }else{
    
    if(current.locus == "2DL23"){
      
      genos.out <- genos.out[!sapply(genos.out, is.null)]
      genos.out <- genos.out[!sapply(genos.out,nrow) == 0]
      
      # To make sure new alleles don't break anything
      if(!any(lapply(genos.out, nrow) == 0) && length(genos.out) > 1){
        
        # Intersecting results
        genos.out <- genos_2DL23(genos.out)
      }else if(all(lapply(genos.out, nrow)) == 0){
        return(NULL)
      }
    }
    
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
    
    ## NEW ALLELES altered function: New alleles output converted into a genos string ------------
    newalleles.output.df.list <- lapply(no_types.list, new_alleles_genos.vcf, current.locus=current.locus,DPthresh=DPthresh,mismatchThresh = 1)
    
    # unlist and combine all new allele info for this sample and locus
    newalleles.output.colnames  <- names(newalleles.output.df.list[[1]])
    newalleles.output.df        <- data.frame(matrix(nrow = 0, ncol = length(newalleles.output.colnames)))
    names(newalleles.output.df) <- newalleles.output.colnames
    if (!is.null(unlist(newalleles.output.df))){
      newalleles.output.df <- do.call("rbind", newalleles.output.df.list)
    }
    
    # Output new alleles results along with genos
    if(file.exists(paste0(results.directory, "KIRcaller/snps_newalleles_", current.locus, ".txt"))) {
      write.table(newalleles.output.df, paste0(results.directory, "KIRcaller/snps_newalleles_", current.locus, ".txt"), append = T, row.names=F, col.names=F, quote=F)
    }else{
      write.table(newalleles.output.df, paste0(results.directory, "KIRcaller/snps_newalleles_", current.locus, ".txt"), append = F,row.names=F, col.names=T, quote=F)
    }
    
    # Write 'new' into genos output
    write.table(paste0(sample_id," new"), paste0(results.directory, "KIRcaller/genos_", current.locus, ".txt"), append = T, row.names=F, col.names=F, quote=F)
    
    ##
    
    newalleles.out <- lapply(no_types.list, new_alleles.vcf, DPthresh=DPthresh)
    
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
  
  newsnps.out <- lapply(vcf_list.good, new_snps_call.vcf,DPthresh = DPthresh)
  newsnps.df <- as.data.frame(do.call(rbind, newsnps.out))
  
  if(file.exists(paste0(results.directory, "KIRcaller/newsnps_", current.locus, ".txt"))) {
    write.table(newsnps.df, paste0(results.directory, "KIRcaller/newsnps_", current.locus, ".txt"), append = TRUE, row.names=T, col.names=F, quote=F)
  }else{
    write.table(newsnps.df, paste0(results.directory, "KIRcaller/newsnps_", current.locus, ".txt"), row.names=T, col.names=T, quote=F)
  }
  
  
  snps.out <- lapply(vcf_list.good, get_calls.vcf, DPthresh=DPthresh)
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

