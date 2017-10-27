# Copyright 2016 Wesley Marin, Paul Norman, Jill Hollenbach
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

ping_gc_kff_version <- '1.2'
cat(paste0('PING_gc_kff version: ', ping_gc_kff_version))

ping_gc_kff <- function(
  sample.location = "PING_sequences/",
  fastq.pattern.1 = "_1.fastq.gz",
  fastq.pattern.2 = "_2.fastq.gz",
  threshold.KFF = 0.2,
  read.cap = 120000,
  results.directory = "",
  probelist.file = 'Resources/gc_resources/probelist_2017_09_19.csv'
  ) {
  
  library(data.table)

  
  # Build sequence list ------------------------------------------------------
  get_sample_list <- function() {
    
    # This is the non-recursive version
    sample_list = list.files(file.path(sample.location), pattern = fastq.pattern.1)
    
    if (is.na(sample_list[1])) {
      string <- paste("No sequences found in", sample.location, "using fastq pattern", fastq.pattern.1)
      stop(string)
    } else {
      sample_list <- gsub(fastq.pattern.1, "", sample_list)
      cat(paste("Found sequences: ", paste(sample_list, collapse = "\n"), sep = "\n"))
      cat("\n")
      return(sample_list)
    }
  }
  
  # Creates results directory, defaults to GC_results
  create_results_directory <- function() {
    cat("----- Getting PING ready -----\n\n")
    cat(paste0("Current working directory: ", getwd(), '\n\n'))
    
    if(results.directory != ""){
      save_to <- results.directory
    }else{
      save_to <- paste0("GC_KFF_results_", format(Sys.time(), "%Y_%m_%d_%H_%M"), "/")
      
      count <- 1
      while(file.exists(save_to)) {
        save_to <- paste0("GC_KFF_results_", format(Sys.time(), "%Y_%m_%d_%H_%M"), "_", count, "/")
        count <- count + 1
      }
    }
    dir.create(file.path(save_to), showWarnings = F)
    save_to <- normalizePath(save_to)
    
    cat(paste("Results being saved to", save_to, "\n\n"))
    
    return(save_to)
  }
  
  # Build primer table from forward and reverse primer files
  make_primer_table <- function(primerlist.file = probelist.file) {
    primer_table = read.csv(primerlist.file)
    
    primer_table$Sequence <- gsub(" ", "", primer_table[["Sequence"]])
    
    return(primer_table)
  }
  
  # Counting primer matches to a table
  count_primer_matches <- function(primer.table, sample, is.gz) {
    primer_count_table <- data.frame(matrix(0, nrow=length(primer.table[,1]), ncol=1))
    rownames(primer_count_table) <- primer.table$Name
    colnames(primer_count_table) <- sample
    
    
    # Counting primer matches in first fastq file
    cat(paste0(sample, fastq.pattern.1, "\n"))
        
    if(is.gz){
      inFile <- file.path(sample.location, paste0(sample, fastq.pattern.1))
      file_contents <- fread(paste("zcat", inFile), sep="\n", nrows = read.cap, header=FALSE)
    }else{
      inFile <- file.path(sample.location, paste0(sample, fastq.pattern.1))
      file_contents <- fread(inFile, sep="\n", nrows = read.cap, header=FALSE)
    }
        
    file_contents <- file_contents[seq(2, length(file_contents[[1]]), 4)]
        
    for(j in 1:length(primer.table[,"Sequence"])) {
      primer_count_table[j, sample] <- sum(length(grep(primer.table[j, "Sequence"], file_contents[[1]], fixed=TRUE)), primer_count_table[j, sample])
    }
        
    remove(file_contents)
    
    # Counting primer matches in second fastq file
    cat(paste0(sample, fastq.pattern.2, "\n"))
    
    if(is.gz){
      inFile <- file.path(sample.location, paste0(sample, fastq.pattern.2))
      file_contents <- fread(paste("zcat", inFile), sep="\n", nrows = read.cap, header=FALSE)
    }else{
      inFile <- file.path(sample.location, paste0(sample, fastq.pattern.2))
      file_contents <- fread(inFile, sep="\n", nrows = read.cap, header=FALSE)
    }
    
    file_contents <- file_contents[seq(2, length(file_contents[[1]]), 4)]
    
    for(j in 1:length(primer.table[,"Sequence"])) {
      primer_count_table[j, sample] <- sum(length(grep(primer.table[j, "Sequence"], file_contents[[1]], fixed=TRUE)), primer_count_table[j, sample])
    }
    
    remove(file_contents)
    
    return(primer_count_table)
  }
  
  # Normalize and check against threshold
  reduce_and_normalize_primer_counts <- function(all.primer.count.table, sample) {
    
    # Copy only the locus presence probes into a new table (the rest are for allele calling)
    locus_probes_index <- grep(">", rownames(all.primer.count.table))
    locus_probes_names <- rownames(all.primer.count.table)[locus_probes_index]
    locus_primer_count_table <- data.frame(all.primer.count.table[locus_probes_names,],
                                           row.names=locus_probes_names)
    colnames(locus_primer_count_table) <- sample
    
    # Split the forward probes from the reverse probes
    non_rc_names <- unlist(strsplit(rownames(locus_primer_count_table), '_rc'))
    non_rc_names <- unique(non_rc_names)
    
    # Create a new table with only the forward probes
    locus_primer_sum_table <- data.frame(matrix(0, nrow=length(non_rc_names), ncol=1), 
                                         row.names=non_rc_names,
                                         stringsAsFactors=FALSE)
    colnames(locus_primer_sum_table) <- sample
    
    # Add the sum of the forward and reverse probes to the sum table
    for(rowname in rownames(locus_primer_sum_table)){
      rowname_rc <- paste0(rowname,'_rc')
      rowname_list <- c(rowname, rowname_rc)
      
      locus_primer_sum_table[rowname,sample] <- sum(locus_primer_count_table[rowname_list,sample])
    }
    
    # Split the probe names into their base locus names, assemble a locus list based off this
    locus_names <- tstrsplit(rownames(locus_primer_sum_table), '>')[[2]]
    locus_names <- unique(locus_names)
    
    # Create a new table with only the locus list
    normalized_table <- data.frame(matrix(0, nrow=length(locus_names), ncol=1),
                                   row.names=locus_names,
                                   stringsAsFactors=FALSE)
    colnames(normalized_table) <- sample
    
    # Normalized the mean of the probes of each loci to the mean of 3DL3
    # Also check if the ratio of means is above the KFF threshold, save the boolean
    for(locus in locus_names){
      locus_list <- grep(paste0('>',locus,'>'), rownames(locus_primer_sum_table))
      locus_3DL3 <- grep('>3DL3>', rownames(locus_primer_sum_table))

      locus_mean <- sum(locus_primer_sum_table[locus_list,sample]) / length(locus_list)
      locus_3DL3_mean <- sum(locus_primer_sum_table[locus_3DL3,sample]) / length(locus_3DL3)
      normalized_table[locus,sample] <- ((locus_mean / locus_3DL3_mean)>=threshold.KFF)*1
    }
    
    # Return this boolean table
    return(normalized_table)
  }
  
  # Setting results directory and trimming reads
  # Find out if the first file is gzipped or not
  find_gz <- function(fastq.pattern){
    if(is.na(fastq.pattern)){
      return(FALSE)
    }else{
      gz <- last(unlist(strsplit(fastq.pattern, ".", fixed = TRUE))) == "gz"
      return(gz)
    }
  }
  
  # Setting up run environment -------------------------------
  
  results_directory <- create_results_directory()
  
  sample_list <- get_sample_list()
  
  # Run KFF -------------------------------------------------

  cat("\n")
  cat("----- Running KFF -----\n")
    
  cat("\n")
  cat("Building probe table. \n")
  primer_table <- make_primer_table()
  is_gz <- find_gz(fastq.pattern.1)
  for(sample in sample_list){
    
    
    cat("\n")
    cat("Counting probes in: \n")
    primer_match_table <- count_primer_matches(primer_table, sample, is_gz)
    
    file_name <- file.path(results_directory,'raw_kff_counts.csv')
    
    if(file.exists(file_name)){
      cat("\nAppending probe matches to raw_kff_counts.csv\n")
      write.table(t(primer_match_table), file = file_name, sep=',', append=TRUE, col.names = FALSE)
    }else{
      cat("\nWriting probe matches to raw_kff_counts.csv\n")
      write.table(t(primer_match_table), file = file_name, sep=',')
    }
  
      
    cat("\n")
    cat("Reducing and normalizing probe counts. \n")
    normalized_results <- reduce_and_normalize_primer_counts(primer_match_table, sample)
    
    file_name <- file.path(results_directory,'kff_results.csv')
    
    if(file.exists(file_name)){
      cat("Appending presence / absence information to kff_results.csv\n")
      write.table(t(normalized_results), file = file_name, sep=',', append=TRUE, col.names = FALSE)
    }else{
      cat("Writing presence / absence information to kff_results.csv\n")
      write.table(t(normalized_results), file = file_name, sep=',')
    }
    cat("\n----------------------------------------------------------------\n")
  }
  cat("\n")
  cat("Finished! Please look at kff_results.csv.")
}
