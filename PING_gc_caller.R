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

ping_gc <- function(
  run.MIRA = TRUE, 
  run.KFF = TRUE, 
  make.graphs = TRUE, 
  sample.location = "PING_sequences/", 
  threshold.file = "Resources/gc_resources/defaultThresholds.txt", 
  threshold.KFF = 0.2, 
  read.cap = 120000,
  results.directory = ""
  ) {
  
  library(data.table)
  
  
  # Build sequence list ------------------------------------------------------
  
  get_sequence_list <- function(folder.name = sample.location, file.pattern = "_1.fastq") {
    sequence_list = list.files(file.path(folder.name), pattern = file.pattern)
    if (is.na(sequence_list[1])) {
      stop("No sequences found, please place fastq files in the PING_sequences folder.")
    } else {
      sequence_list <- gsub(file.pattern, "", sequence_list)
      cat(paste("Found sequences: ", paste(sequence_list, collapse = "\n"), sep = "\n"))
      return(sequence_list)
    }
  }
  
  # Finds files that are smaller than the recommended size
  files_too_small <- function(sequence.list){
    small_samples <- NA
    
    for(sample in sequence.list){
      line_count <- as.numeric(system2("wc", c("-l", "<", paste0(sample.location, sample, "_1.fastq")), stdout = T))
      if(line_count < read.cap){
        small_samples <- c(small_samples, paste0(sample.location, sample, "_1.fastq"))
      }
      
      line_count <- as.numeric(system2("wc", c("-l", "<", paste0(sample.location, sample, "_2.fastq")), stdout = T))
      if(line_count < read.cap){
        small_samples <- c(small_samples, paste0(sample.location, sample, "_2.fastq"))
      }
    }
    
    small_samples <- small_samples[!is.na(small_samples)]
    
    if(length(small_samples) > 0){
      error.string <- "\n\nThese files have fewer than the recommended lines, excluding these sequences:\n"
      cat(error.string)
      
      cat(paste(small_samples), sep='\n')
      
      cat("\nContinuing PING_gc_caller analysis...\n\n")
    }
    
    return(small_samples)
  }
  
  # Filters sequence_list, taking out samples that are too small
  filter_sequences <- function(too.small, sequence.list){
    bad_sequences <- unlist(strsplit(too.small, sample.location))
    bad_sequences <- bad_sequences[bad_sequences != ""]
    bad_sequences <- unlist(strsplit(bad_sequences, "_1.fastq"))
    bad_sequences <- unlist(strsplit(bad_sequences, "_2.fastq"))
    
    bad_sequences <- unique(bad_sequences)
    
    sequence.list <- sequence.list[-match(bad_sequences, sequence.list)]
    
    return(sequence.list)
  }
  
  sequence_list <- get_sequence_list()
  
  too_small <- files_too_small(sequence_list)
  
  if(length(too_small) > 0){
    sequence_list <- filter_sequences(too_small, sequence_list)
  }
  
  results_directory <- function() {
    
    if(results.directory != ""){
      save_to <- results.directory
    }else{
    
      save_to <- paste0("GC_results_", format(Sys.time(), "%Y_%m_%d_%H:%M"), "/")
      
      count <- 1
      while(file.exists(save_to)) {
        save_to <- paste0("GC_results_", format(Sys.time(), "%Y_%m_%d_%H:%M"), "_", count, "/")
        count <- count + 1
      }
    }
    
    dir.create(save_to)
    return(save_to)
  }
  
  # Caps sequences at read.cap lines
  cut_fastq <- function(file.name, read.cap, post.file.name) {
    file_contents <- fread(file.name, sep="\n", nrows = read.cap, header = F)
    write.table(file_contents, file = post.file.name, quote = F, row.names = F, col.names = F)
  }
  
  # KFF functions -----------------------------------------------------------
  
  # Build primer table from forward and reverse primer files
  
  make_primer_table <- function(primerlist.file = "Resources/gc_resources/Primerlist.txt") {
    primer_table = read.delim(primerlist.file, header = FALSE)
    colnames(primer_table) <- c("Locus", "Primer")
    
    primer_table$Primer <- gsub(" ", "", primer_table[["Primer"]])
    
    return(primer_table)
  }
  
  
  # Counting primer matches to a table
  
  count_primer_matches <- function(primer.table, sequence.list) {
    primer_count_table <- data.frame(matrix(0, nrow=length(primer.table[,1]), ncol=length(sequence.list)))
    rownames(primer_count_table) <- primer.table$Locus
    colnames(primer_count_table) <- sequence.list
    
    counter <- 1
    while(counter <= 2) {
      
      for(i in 1:length(sequence.list)) {
        cat(paste0(sequence.list[i], "_", counter, "\n"))
        inFile <- paste0(sample.location, sequence.list[i], "_", counter, ".fastq")
        file_contents <- fread(inFile, sep="\n", nrows = read.cap, header=FALSE)
        file_contents <- file_contents[seq(2, length(file_contents[[1]]), 4)]
        
        for(j in 1:length(primer.table[,"Primer"])) {
          primer_count_table[j, sequence.list[i]] <- sum(length(grep(primer.table[j, "Primer"], file_contents[[1]], fixed=TRUE)), primer_count_table[j, sequence.list[i]])
        }
        
        remove(file_contents)
      }
      
      counter = counter + 1
    }
    
    return(primer_count_table)
  }
  
  
  # Normalize and check against threshold
  
  reduce_and_normalize_primer_counts <- function(primer.count.table, sequence.list) {
    
    
    # Adding forward and reverse primer counts
    
    for(i in 1:length(primer.count.table[,1])) {
      primer.count.table[rownames(primer.count.table)[i],] <- primer.count.table[rownames(primer.count.table)[i],] + primer.count.table[paste0(rownames(primer.count.table)[i],"rc"),]
    }
    
    primer_rows <- rownames(primer.count.table)
    primer_rows <- primer_rows[1:(length(primer_rows)/2)]
    
    
    # Making a new frame from added values
    
    reduced_table <- as.data.frame(primer.count.table[1:(length(primer.count.table[,1])/2),], row.names = primer_rows)
    colnames(reduced_table) <- sequence.list
    
    loci_list <- unlist(strsplit(rownames(reduced_table), "_"))
    loci_list <- unique(loci_list[grep(">", loci_list)])
    
    normalized_table <- data.frame(matrix(0, nrow=length(loci_list), ncol=length(sequence.list)))
    rownames(normalized_table) <- loci_list
    colnames(normalized_table) <- sequence.list
    
    
    # Build further reduced table by combining counts for each locus
    
    for(j in 1:length(loci_list)) {
      for(i in 1:length(sequence.list)) {
        normalized_table[loci_list[j], sequence.list[i]] <- sum(reduced_table[grep(paste0(loci_list[j],"_"), rownames(reduced_table)), sequence.list[i]])
      }
      normalized_table[loci_list[j],] <- normalized_table[loci_list[j],] / (length(grep(loci_list[j], rownames(reduced_table)))*2)
    }
    
    for(j in 1:length(sequence.list)) {
      
      
      # Here we are normalizing against 3DL3 count
      
      normalized_table[,sequence.list[j]] <- normalized_table[,sequence.list[j]]/normalized_table[">3DL3", sequence.list[j]]
      
      
      # And then checking to see if those results are above the set threshold value
      
      normalized_table[,sequence.list[j]] <- (normalized_table[,sequence.list[j]] > threshold.KFF) * 1
    }
    
    
    # Taking out all the arrows from the loci list
    
    new_list <- unique(unlist(strsplit(loci_list, ">")))
    new_list <- new_list[2:length(new_list)]
    rownames(normalized_table) <- sapply(strsplit(rownames(normalized_table), ">"), '[', 2)
    
    return(normalized_table)
  }
  
  
  # Setting results directory and trimming reads
  
  results <- results_directory()
  cat("\n\n")
  cat("Setting results directory to", results)
  
  # Run KFF -------------------------------------------------
  
  if(run.KFF) {

    cat("\n\n")
    cat("----- Running KFF -----\n")
    
    cat("\n")
    cat("Building primer table. \n")
    primer_table <- make_primer_table()
    
    cat("\n")
    cat("Counting primers in: \n")
    primer_match_table <- count_primer_matches(primer_table, sequence_list)
    
    cat("\n")
    cat("Reducing and normalizing primer counts. \n")
    normalized_results <- reduce_and_normalize_primer_counts(primer_match_table, sequence_list)
    
    write.csv(normalized_results, file = paste0(results, "kff_results.csv"))
    
    cat("\n")
    cat("Finished! Please look at kff_results.csv.")
    
  }
  
  
  # MIRA functions ------------------------------------------------------------
  
  ## This is an experimental section meant to differentiate kff neg results in MIRA graphs
  check_kff_results <- function(){
    if(file.exists(paste0(results, "kff_results.csv"))){
      kff_results <- read.csv(paste0(results, "kff_results.csv"))
    }
    
    return(kff_results)
  }
  
  # This part runs MIRA on all sequences in the Sequence folder, based on preferences set in the conf file
  # Returns a table of MIRA read counts
  
  build_mira_counts <- function(sequence.list, mira.reference = "Resources/gc_resources/MIRA_configuration.conf") {
    dir.create(".miraReads", showWarnings = F)
    file.remove(file.path(".miraReads", list.files(".miraReads")))
    
    for (i in 1:length(sequence.list)) {
      
      cut_fastq(paste0(sample.location, sequence.list[i], "_1.fastq"), read.cap, paste0(".miraReads/", sequence.list[i], "_1.fastq"))
      cut_fastq(paste0(sample.location, sequence.list[i], "_2.fastq"), read.cap, paste0(".miraReads/", sequence.list[i], "_2.fastq"))
      
      
      cat(paste0(sequence.list[i], "\n"))
      
      system2(system2("which", "mira", stdout = TRUE), mira.reference, stdout = "MIRA_log.txt")
      
      file.remove(file.path(".miraReads", list.files(".miraReads")))
      
      mira_filepath <- "KIR_allele_assembly/KIR_allele_d_info/KIR_allele_info_contigstats.txt"
      
      if (!file.exists(file.path(mira_filepath))) {
        unlink(file.path("miraReads", recursive = TRUE))
        unlink(file.path("KIR_allele_assembly"), recursive = TRUE)
        stop("Program stopped, please check MIRA_log.txt for more information")
      }
      
      mira_count_file <- read.delim(file.path(mira_filepath))
      
      if (i == 1) {
        mira_count_table <- data.frame(matrix(NA, nrow=length(mira_count_file[,1]), ncol=length(sequence.list)))
        rownames(mira_count_table) <- mira_count_file$X..name
        colnames(mira_count_table) <- sequence.list
        mira_count_table[sequence.list[i]] <- mira_count_file$X..reads
      } else {
        mira_count_table[sequence.list[i]] <- mira_count_file$X..reads
      }
    }
    
    unlink(file.path(".miraReads"), recursive = TRUE)
    unlink(file.path("KIR_allele_assembly"), recursive = TRUE)
    
    return(mira_count_table)
  }
  
  
  # Normalizing MIRA read counts to 3DL3
  
  normalize_mira_reads <- function(sequence.list, mira.count.table) {
    for(i in 1:length(mira.count.table[1,])) {
      mira.count.table[sequence.list[i]] <- mira.count.table[sequence.list[i]]/mira.count.table[grep("3DL3", rownames(mira.count.table)),sequence.list[i]]
    }
    
    return(mira.count.table)
  }

  
  # Creating graphs and allowing user to manually set threshold values
  # Returns table of threshold values
  
  create_graphs <- function(normalized.mira.counts, threshold.file) {
    
    threshold_table <- data.frame(matrix(NA, nrow=15, ncol=11))
    default_table <- read.delim(threshold.file)
    threshold_table[,1:length(default_table)] <- default_table
    
    if(file.exists(paste0(results, "kff_results.csv"))){
      # Loading KFF results
      kff_results <- check_kff_results()
    }
    
    for(i in 1:length(normalized.mira.counts[,1])) {
      flag <- 'n'
      
      while(flag == 'n') {
  
        x_lim <- length(normalized.mira.counts[i,])
        y_lim <- max(normalized.mira.counts[i,]) + 0.2
        
        par(bg = "white")
        
        plot(as.numeric(normalized.mira.counts[i, order(normalized.mira.counts[i,])]), ylab = "Locus Ratio", xlab = "Sample", xlim = c(1, x_lim), ylim = c(0, y_lim))
        
        ## Making KFF neg results show as red X's
        if(file.exists(paste0(results, "kff_results.csv"))){

          # Finding the ordered colnames of mira results
          kff_sample_order <- colnames(normalized.mira.counts)[order(normalized.mira.counts[i,])]
          
          # Matching mira and kff locus
          kff_locus <- rownames(normalized.mira.counts)[i]
          kff_locus <- unlist(strsplit(kff_locus, "_"))[1]
          
          # Putting kff results in the same order
          kff_locus_results <- kff_results[grep(kff_locus, kff_results[,1]), kff_sample_order]
          
          # Finding all KFF neg results
          kff_neg_samples <- colnames(kff_locus_results)[grepl(0, kff_locus_results)]
          
          if(length(kff_neg_samples) != 0){
            # Converting these back into MIRA points
            mira_kff_neg_points <- normalized.mira.counts[i, kff_neg_samples]
            
            # Graphing
            par(new = T)
            plot((1:length(mira_kff_neg_points)), mira_kff_neg_points, ylab = NA, xlab = NA, axes = F, col = "red", pch = 4, xlim = c(1, x_lim), ylim = c(0, y_lim))
          }
        }
        
        title(paste0(rownames(normalized.mira.counts)[i],"/ 3DL3"))
        abline(h = threshold_table[i,2:sum((!is.na(threshold_table[i,]))*1)], col = "blue", lty=2)
        
        threshold_number <- 0
        
        while(!(threshold_number >= 1 & threshold_number <= 10) && !is.na(as.numeric(threshold_number))) {
          threshold_number <- readline("How many threshold values for this locus? Will take values 1 - 10 (or Enter to keep defaults):  ")
          threshold_number <- as.numeric(threshold_number)
        }
        
        if(is.na(threshold_number)) {
          thresholds <- threshold_table[i,2:sum((!is.na(threshold_table[i,]))*1)]
        } else if(threshold_number < 1 | threshold_number > 10) {
          cat("Please select a value between 1 and 10")
          return    
        } else {
          cat(paste0("\nPlease click ", threshold_number, " points on the graph.\n"))
          thresholds <- round(locator(n = threshold_number, type = "n")$y, digits = 3)
          thresholds <- sort(thresholds)
        }
        
        abline(h = thresholds, col = "red")
        
        thresholds_legend <- thresholds
        
        for(j in 1:length(thresholds)) {
          thresholds_legend[j] <- paste0("Copy number ", j - 1, " is < ", thresholds[j])
        }
        
        legend("bottomright", as.character(thresholds_legend), title="Thresholds")
        
        cat("\n")
        flag <- readline("Would you like to keep these threshold values (y/n)? ")
        flag <- tolower(flag)
        
        cat("\n--------------------------------\n\n")
        if(flag == 'n') {
          dev.off()
        } else {
          threshold_table[i,2:length(threshold_table[i,])] <- NA
          threshold_table[i, 2:(length(thresholds) + 1)] <- thresholds
          
          save_and_remove_graphs(results, threshold_table[i, 1])
          
        }
      }
    }
    
    return(threshold_table)
  }
  
  save_and_remove_graphs <- function(results.directory, graph.name) {
    count <- 1
    
    #dev.set(which = dev.next())
    
    while(dev.interactive()){
      dev.copy(pdf, paste0(results.directory, graph.name[count], ".pdf"))
      dev.off()
      dev.off()
      count <- count + 1
    }
  }
  
  
  # Determines locus counts based on threshold values
  # Returns table with locus counts
  
  determine_locus_counts <- function(threshold.table, normalized.mira.counts) {
    mira_locus_counts <- normalized.mira.counts
    
    for(j in 1:length(normalized.mira.counts[1,])) {
      
      for(i in 1:length(normalized.mira.counts[,1])) {
        in_count <- 0
        
        for(k in 2:length(threshold.table[1,])) {
          row_number <- grep(threshold.table[i,1], rownames(normalized.mira.counts))
          
          if((!is.na(threshold.table[i,k]) && !is.null(threshold.table[i,k]) && (normalized.mira.counts[row_number, j] >= threshold.table[i,k]))) {
            in_count <- in_count + 1
          } else {
            mira_locus_counts[row_number, j] <- in_count
          }
        }
      }
    }
    
    return(mira_locus_counts)
  }
  
  
  # Run MIRA ----------------------------------------------------------------
  
  if(run.MIRA) {

    cat("\n\n\n")
    cat("----- Running MIRA -----\n")
    prelimResults <- ""
    
    
    # Running MIRA and reading in results -------------------------------------
    
    cat("\nDetermining read counts for:\n")
    mira_count_table <- build_mira_counts(sequence_list)
    
    if(!exists("mira_count_table")) {
      stop("Ran into problems running MIRA. Please check MIRA_log.txt")
    }
    cat("Done.\n")
    
    ## Adding the results folder name to the MIRA_count_table file
    mira_count_table_recalc <- mira_count_table
    
    add_result_row <- data.frame(matrix(data = NA, nrow = 1, ncol = ncol(mira_count_table_recalc)))
    add_result_row[1,1] <- results
    
    colnames(add_result_row) <- colnames(mira_count_table_recalc)
    mira_count_table_recalc <- rbind(mira_count_table_recalc, add_result_row)
    
    
    write.csv(mira_count_table_recalc, file = "MIRA_count_table.csv")
    write.csv(mira_count_table, file = paste0(results, "MIRA_count_table.csv"))
    
    cat("\nNormalizing MIRA reads to 3DL3\n")
    normalized_mira_table <- normalize_mira_reads(sequence_list, mira_count_table)

    cat("Done.\n")
    
    
    # Creating graphs and setting thresholds ---------------------------------- 
    
    if(make.graphs){
      cat("\n\nATTENTION!--------------------------------------------------------------\n")
      cat("\nGenerating graphs, default threshold values will be marked with blue lines. \n")
      cat("\nYou will be able to mark threshold values for each locus by entering how many threshold values you would like\nand clicking where the threshold values should be on the graph\n")
      cat("\nDefault threshold values will be marked with a blue dotted line,\nenter 0 for the number of threshold values if you would like to keep these default values\n\n")
      
      while(!exists("threshold_table")) {
        threshold_table <- create_graphs(normalized_mira_table, threshold.file)
      }
    } else {
      threshold_table <- data.frame(matrix(NA, nrow=15, ncol=11))
      default_table <- read.delim(threshold.file)
      threshold_table[,1:length(default_table)] <- default_table
    }
    
    locus_count_table <- determine_locus_counts(threshold_table, normalized_mira_table)
    rownames(locus_count_table) <- gsub("_.*$", "", rownames(locus_count_table))
    
    # Deleting graphs and saving results --------------------------------------
    
    write.csv(locus_count_table, file=paste0(results,"MIRA_results.csv"))
    cat("\nFinished! Please look at MIRA_results.csv.\n") 
    
  } # End of MIRA
  
  
  if(run.MIRA && run.KFF) {
    mira_results <- read.csv(paste0(results, "MIRA_results.csv"), row.names = 1, check.names = FALSE)
    kff_results <- read.csv(paste0(results, "kff_results.csv"), row.names = 1, check.names = FALSE)
    
    com_names <- unique(c(rownames(mira_results), rownames(kff_results)))
    com_table <- merge(kff_results, mira_results, by=0, suffixes = c(".kff",".mira"))
    final_table <- data.frame(matrix(0, nrow=length(com_names), ncol=length(sequence_list)))
    colnames(final_table) <- sequence_list
    rownames(final_table) <- com_names
    
    for(j in 1:length(sequence_list)) {
      for(i in 1:length(com_names)) {
        kff_value <- com_table[i, paste0(sequence_list[j], ".kff")]
        mira_value <- com_table[i, paste0(sequence_list[j], ".mira")]
        if( (kff_value && mira_value) || kff_value == mira_value ) {
          final_table[com_table[i,"Row.names"],sequence_list[j]] <- mira_value
        } else {
          final_table[com_table[i,"Row.names"],sequence_list[j]] <- paste0(com_table[i,paste0(sequence_list[j],".kff")],"_",com_table[i,paste0(sequence_list[j],".mira")])
        }
      }
    }
    write.csv(final_table, file = paste0(results, "Combined_results.csv"))
  }

}

ping_recalc <- function(mira.csv = "MIRA_count_table.csv", threshold.file = "Resources/gc_resources/defaultThresholds.txt", make.graphs = TRUE) {
  
  results_directory <- function() {
    save_to <- paste0("GC_results_", format(Sys.time(), "%Y_%m_%d_%H:%M"), "/")
    
    count <- 1
    while(file.exists(save_to)) {
      save_to <- paste0("GC_results_", format(Sys.time(), "%Y_%m_%d_%H:%M"), "_", count, "/")
      count <- count + 1
    }
    
    dir.create(save_to)
    return(save_to)
  }
  
  ## This is an experimental section meant to differentiate kff neg results in MIRA graphs
  check_kff_results <- function(){
    if(file.exists(paste0(old.results, "kff_results.csv"))){
      kff_results <- read.csv(paste0(old.results, "kff_results.csv"))
    }
    
    return(kff_results)
  }

  # Normalizing MIRA read counts to 3DL3
  
  normalize_mira_reads <- function(sequence.list, mira.count.table) {
    for(i in 1:length(mira.count.table[1,])) {
      
      count_3DL3 <- mira.count.table[grep("3DL3", rownames(mira.count.table)), sequence.list[i]]
      count_3DL3 <- as.numeric(as.character(count_3DL3))
      
      mira.count.table[sequence.list[i]] <- mira.count.table[sequence.list[i]]/count_3DL3
    }
    
    return(mira.count.table)
  }
  
  
  # Creating graphs and allowing user to manually set threshold values
  # Returns table of threshold values
  
  create_graphs <- function(normalized.mira.counts, threshold.file) {
    
    threshold_table <- data.frame(matrix(NA, nrow=15, ncol=11))
    default_table <- read.delim(threshold.file)
    threshold_table[,1:length(default_table)] <- default_table
    
    # Loading KFF results
    if(file.exists(paste0(old.results, "kff_results.csv"))){
      kff_results <- check_kff_results()
    }
    
    for(i in 1:length(normalized.mira.counts[,1])) {
      flag <- 'n'
        
      while(flag == 'n') {
        
        x_lim <- length(normalized.mira.counts[i,])
        y_lim <- max(normalized.mira.counts[i,]) + 0.2
          
        par(bg = "white")
        
        plot(as.numeric(normalized.mira.counts[i, order(normalized.mira.counts[i,])]), ylab = "Locus Ratio", xlab = "Sample", xlim = c(1, x_lim), ylim = c(0, y_lim))
        
        ## Making KFF neg results show as red X's
        if(file.exists(paste0(old.results, "kff_results.csv"))){
          # Finding the ordered colnames of mira results
          kff_sample_order <- colnames(normalized.mira.counts)[order(normalized.mira.counts[i,])]
          
          # Matching mira and kff locus
          kff_locus <- rownames(normalized.mira.counts)[i]
          kff_locus <- unlist(strsplit(kff_locus, "_"))[1]
          
          # Putting kff results in the same order
          kff_locus_results <- kff_results[grep(kff_locus, kff_results[,1]), kff_sample_order]
          
          # Finding all KFF neg results
          kff_neg_samples <- colnames(kff_locus_results)[grepl(0, kff_locus_results)]
          
          if(length(kff_neg_samples) != 0){
            # Converting these back into MIRA points
            mira_kff_neg_points <- normalized.mira.counts[i, kff_neg_samples]
            
            # Graphing
            par(new = T)
            plot((1:length(mira_kff_neg_points)), mira_kff_neg_points, ylab = NA, xlab = NA, axes = F, col = "red", pch = 4, xlim = c(1, x_lim), ylim = c(0, y_lim))
          }
        }
        
        title(paste0(rownames(normalized.mira.counts)[i],"/ 3DL3"))
        abline(h = threshold_table[i,2:sum((!is.na(threshold_table[i,]))*1)], col = "blue", lty=2)
          
        threshold_number <- 0
          
        while(!(threshold_number >= 1 & threshold_number <= 10) && !is.na(as.numeric(threshold_number))) {
          threshold_number <- readline("How many threshold values for this locus? Will take values 1 - 10 (or Enter to keep defaults):  ")
          threshold_number <- as.numeric(threshold_number)
        }
          
        if(is.na(threshold_number)) {
          thresholds <- threshold_table[i,2:sum((!is.na(threshold_table[i,]))*1)]
        } else if(threshold_number < 1 | threshold_number > 10) {
          cat("Please select a value between 1 and 10")
          return    
        } else {
          cat(paste0("\nPlease click ", threshold_number, " points on the graph.\n"))
          thresholds <- round(locator(n = threshold_number, type = "n")$y, digits = 3)
          thresholds <- sort(thresholds)
        }
          
        abline(h = thresholds, col = "red")
        
        thresholds_legend <- thresholds
        
        for(j in 1:length(thresholds)) {
          thresholds_legend[j] <- paste0("Copy number ", j - 1, " is < ", thresholds[j])
        }
        
        legend("bottomright", as.character(thresholds_legend), title="Thresholds")
          
        cat("\n")
        flag <- readline("Would you like to keep these threshold values (y/n)? ")
        flag <- tolower(flag)
          
        cat("\n--------------------------------\n\n")
        if(flag == 'n') {
          dev.off()
        } else {
          threshold_table[i,2:length(threshold_table[i,])] <- NA
          threshold_table[i, 2:(length(thresholds) + 1)] <- thresholds
          
          save_and_remove_graphs(results, threshold_table[i, 1])
          
        }
      }
    }

    return(threshold_table)
  }
  
  save_and_remove_graphs <- function(results.directory, graph.name) {
    count <- 1
    
    #dev.set(which = dev.next())
    
    while(dev.interactive()){
      dev.copy(pdf, paste0(results.directory, graph.name[count], ".pdf"))
      dev.off()
      dev.off()
      count <- count + 1
    }
  }
  
  
  # Determines locus counts based on threshold values
  # Returns table with locus counts
  
  determine_locus_counts <- function(threshold.table, normalized.mira.counts) {
    mira_locus_counts <- normalized.mira.counts
    
    for(j in 1:length(normalized.mira.counts[1,])) {
      
      for(i in 1:length(normalized.mira.counts[,1])) {
        in_count <- 0
        
        for(k in 2:length(threshold.table[1,])) {
          row_number <- grep(threshold.table[i,1], rownames(normalized.mira.counts))
          
          if((!is.na(threshold.table[i,k]) && !is.null(threshold.table[i, k]) && (normalized.mira.counts[row_number, j] >= threshold.table[i,k]))) {
            in_count <- in_count + 1
          } else {
            mira_locus_counts[row_number, j] <- in_count
          }
        }
      }
    }
    
    return(mira_locus_counts)
  }
  
  cat("\n\n\n")
  cat("----- Running MIRA -----\n")
  
  results <- results_directory()
  cat("\n\n")
  cat("Setting results directory to", results)
  
  prelimResults <- ""
  
  mira_count_table <- read.csv(mira.csv, check.names = F)
  
  old.results <- mira_count_table[nrow(mira_count_table), 2]
  mira_count_table <- mira_count_table[1: (nrow(mira_count_table) - 1), ]
  
  row.names(mira_count_table) <- mira_count_table[,1]
  mira_count_table <- mira_count_table[,2:length(mira_count_table[1,])]
  mira_count_table[,1] <- as.numeric(as.character(mira_count_table[,1]))
    
  sequence_list <- colnames(mira_count_table)
  
  # Running MIRA and reading in results -------------------------------------
    
    
  cat("\n\nNormalizing MIRA reads to 3DL3\n")
  normalized_mira_table <- normalize_mira_reads(sequence_list, mira_count_table)
    
  cat("Done.\n")
    
    
  # Creating graphs and setting thresholds ---------------------------------- 
  
  if(make.graphs){
    cat("\n\nATTENTION!--------------------------------------------------------------\n")
    cat("\nGenerating graphs, default threshold values will be marked with blue lines. \n")
    cat("\nYou will be able to mark threshold values for each locus by entering how many threshold values you would like\nand clicking where the threshold values should be on the graph\n")
    cat("\nDefault threshold values will be marked with a blue dotted line,\npress Enter with no input for the number of threshold values if you would like to keep these default values\n\n")
    
    while(!exists("threshold_table")) {
      threshold_table <- create_graphs(normalized_mira_table, threshold.file)
    }
    
    thresh_columns <- c("KIR", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    cat("Writing set thresholds to New_thresholds.csv\n\n")
    write.table(threshold_table, file = "New_thresholds.txt", quote = F, sep = "\t", na = "", row.names = F, col.names = thresh_columns)
    write.table(threshold_table, file = paste0(results, "New_thresholds.txt"), quote = F, sep = "\t", na = "", row.names = F, col.names = thresh_columns)
    
    } else {
    
    threshold_table <- data.frame(matrix(NA, nrow=15, ncol=11))
    default_table <- read.delim(threshold.file)
    threshold_table[,1:length(default_table)] <- default_table
  }
  
  locus_count_table <- determine_locus_counts(threshold_table, normalized_mira_table)
  rownames(locus_count_table) <- gsub("_.*$", "", rownames(locus_count_table))
  
  # Saving results --------------------------------------
  
  write.csv(locus_count_table, file=paste0(results,"MIRA_results.csv"))
  cat("\nFinished! Please look at MIRA_results.csv.") 
  
  # End of MIRA

}
