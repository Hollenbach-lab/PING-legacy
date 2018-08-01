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

ping_gc_caller <- function(
  run.MIRA = TRUE,
  make.graphs = TRUE,
  sample.location = "PING_sequences/",
  threshold.file = "Resources/gc_resources/defaultThresholds.txt",
  threshold.KFF = 0.2,
  read.cap = 120000,
  results.directory = "",
  probelist.file = 'Resources/gc_resources/probelist_2017_09_19.csv',
  fastq.pattern.1 = "_1.fastq.gz",
  fastq.pattern.2 = "_2.fastq.gz"
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
      full_kff_table <- read.csv(file_name,row.names = 1)
      merged_kff_table <- merge(full_kff_table, normalized_results, by = 0,all = T)
      row.names(merged_kff_table) <- merged_kff_table[,1]
      merged_kff_table[,1] <- NULL
      write.csv(merged_kff_table, file = file_name)
    }else{
      cat("Writing presence / absence information to kff_results.csv\n")
      write.csv(normalized_results, file = file_name)
    }
    cat("\n----------------------------------------------------------------\n")
  }
  cat("\n")
  cat("Finished! Please look at kff_results.csv.")
  
  
  # MIRA functions ------------------------------------------------------------

  
  ## This is an experimental section meant to differentiate kff neg results in MIRA graphs
  check_kff_results <- function(){
    if(file.exists(file.path(results_directory, "kff_results.csv"))){
      kff_results <- read.csv(file.path(results_directory, "kff_results.csv"), check.names = F)
    }
    
    return(kff_results)
  }
  
  # Caps sequences at read.cap lines
  cut_fastq <- function(file.name, read.cap, post.file.name, is.gz) {
    
    # Uncompressing if is.gz and moving to .miraReads, this happens one at a time, so space should not be a concern
    if(is.gz){
      file_contents <- fread(paste("zcat", file.name), sep="\n", nrows = read.cap, header = F)
    }else{
      file_contents <- fread(file.name, sep="\n", nrows = read.cap, header = F)
    }
    
    write.table(file_contents, file = post.file.name, quote = F, row.names = F, col.names = F)
  }
  
  # This part runs MIRA on all sequences in the Sequence folder, based on preferences set in the conf file
  # Returns a table of MIRA read counts
  
  build_mira_counts <- function(sequence.list, mira.reference = "Resources/gc_resources/MIRA_configuration.conf") {
    dir.create(".miraReads", showWarnings = F)
    file.remove(file.path(".miraReads", list.files(".miraReads")))
    
    for (i in 1:length(sequence.list)) {
      
      cut_fastq(file.path(sample.location, paste0(sequence.list[i], fastq.pattern.1)), read.cap, paste0(".miraReads/", sequence.list[i], "_1.fastq"), is_gz)
      cut_fastq(file.path(sample.location, paste0(sequence.list[i], fastq.pattern.2)), read.cap, paste0(".miraReads/", sequence.list[i], "_2.fastq"), is_gz)
      
      
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
    
    # Loading KFF results
    if(file.exists(file.path(results_directory, "kff_results.csv"))){
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
        if(file.exists(file.path(results_directory, "kff_results.csv"))){
          # Finding the ordered colnames of mira results
          kff_sample_order <- colnames(normalized.mira.counts)[order(normalized.mira.counts[i,])]
          
          # Matching mira and kff locus
          kff_locus <- rownames(normalized.mira.counts)[i]
          kff_locus <- unlist(strsplit(kff_locus, "_"))[1]
          
          # Putting kff results in the same order
          kff_locus_results <- kff_results[grep(kff_locus, kff_results[,1]), kff_sample_order]
          
          # Finding all KFF neg results
          kff_neg_samples <- colnames(kff_locus_results)[grepl(0, kff_locus_results)]
          
          # Mira values
          mira_values <- as.numeric(normalized.mira.counts[i, order(normalized.mira.counts[i,])])
          
          
          if(length(kff_neg_samples) != 0){
            # Converting these back into MIRA points
            mira_kff_neg_points <- as.numeric(normalized.mira.counts[i, kff_neg_samples])
            
            # X axis values
            x_vals <- match(mira_kff_neg_points, mira_values)
            
            # Graphing
            par(new = T)
            plot(x_vals, mira_kff_neg_points, ylab = NA, xlab = NA, axes = F, col = "red", pch = 4, xlim = c(1, x_lim), ylim = c(0, y_lim))
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
        
        legend("topleft", as.character(thresholds_legend), title="Thresholds")
        
        cat("\n")
        flag <- readline("Would you like to keep these threshold values (y/n)? ")
        flag <- tolower(flag)
        
        cat("\n--------------------------------\n\n")
        if(flag == 'n') {
          dev.off()
        } else {
          threshold_table[i,2:length(threshold_table[i,])] <- NA
          threshold_table[i, 2:(length(thresholds) + 1)] <- thresholds
          
          save_and_remove_graphs(results_directory, threshold_table[i, 1])
          
        }
      }
    }
    
    return(threshold_table)
  }
  
  save_and_remove_graphs <- function(results.directory, graph.name) {
    count <- 1
    
    #dev.set(which = dev.next())
    
    while(dev.interactive()){
      dev.copy(pdf, file.path(results.directory, paste0(graph.name[count], ".pdf")))
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
    mira_count_table <- build_mira_counts(sample_list)
    
    if(!exists("mira_count_table")) {
      stop("Ran into problems running MIRA. Please check MIRA_log.txt")
    }
    cat("Done.\n")
    
    ## Adding the results folder name to the MIRA_count_table file
    mira_count_table_recalc <- mira_count_table
    
    add_result_row <- data.frame(matrix(data = NA, nrow = 1, ncol = ncol(mira_count_table_recalc)))
    add_result_row[1,1] <- results_directory
    
    colnames(add_result_row) <- colnames(mira_count_table_recalc)
    mira_count_table_recalc <- rbind(mira_count_table_recalc, add_result_row)
    
    
    write.csv(mira_count_table_recalc, file = "MIRA_count_table.csv")
    write.csv(mira_count_table, file = file.path(results_directory, "MIRA_count_table.csv"))
    
    cat("\nNormalizing MIRA reads to 3DL3\n")
    normalized_mira_table <- normalize_mira_reads(sample_list, mira_count_table)

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
      
      thresh_columns <- c("KIR", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
      cat("Writing set thresholds to New_thresholds.csv\n\n")
      write.table(threshold_table, file = "New_thresholds.txt", quote = F, sep = "\t", na = "", row.names = F, col.names = thresh_columns)
      write.table(threshold_table, file = file.path(results_directory, "New_thresholds.txt"), quote = F, sep = "\t", na = "", row.names = F, col.names = thresh_columns)
      
    } else {
      threshold_table <- data.frame(matrix(NA, nrow=15, ncol=11))
      default_table <- read.delim(threshold.file)
      threshold_table[,1:length(default_table)] <- default_table
    }
    
    locus_count_table <- determine_locus_counts(threshold_table, normalized_mira_table)
    rownames(locus_count_table) <- gsub("_.*$", "", rownames(locus_count_table))
    
    # Deleting graphs and saving results --------------------------------------
    
    write.csv(locus_count_table, file=file.path(results_directory,"MIRA_results.csv"))
    cat("\nFinished! Please look at MIRA_results.csv.\n") 
    
  } # End of MIRA
  
  
  if(run.MIRA) {
    mira_results <- read.csv(file.path(results_directory, "MIRA_results.csv"), row.names = 1, check.names = FALSE)
    kff_results <- read.csv(file.path(results_directory, "kff_results.csv"), row.names = 1, check.names = FALSE)
    
    com_names <- unique(c(rownames(mira_results), rownames(kff_results)))
    com_table <- merge(kff_results, mira_results, by=0, suffixes = c(".kff",".mira"))
    final_table <- data.frame(matrix(0, nrow=length(com_names), ncol=length(sample_list)))
    colnames(final_table) <- sample_list
    rownames(final_table) <- com_names
    
    for(j in 1:length(sample_list)) {
      for(i in 1:length(com_names)) {
        kff_value <- com_table[i, paste0(sample_list[j], ".kff")]
        mira_value <- com_table[i, paste0(sample_list[j], ".mira")]
        if( (kff_value && mira_value) || kff_value == mira_value ) {
          final_table[com_table[i,"Row.names"],sample_list[j]] <- mira_value
        } else {
          final_table[com_table[i,"Row.names"],sample_list[j]] <- paste0(com_table[i,paste0(sample_list[j],".kff")],"_",com_table[i,paste0(sample_list[j],".mira")])
        }
      }
    }
    write.csv(final_table, file = file.path(results_directory, "Combined_results.csv"))
  }

}

ping_recalc <- function(
  mira.csv = "MIRA_count_table.csv",
  kff.results.directory = "",
  threshold.file = "Resources/gc_resources/defaultThresholds.txt", 
  make.graphs = TRUE,
  results.directory = ""
  ) {
  
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
  
  ## This is an experimental section meant to differentiate kff neg results in MIRA graphs
  check_kff_results <- function(){
    if(file.exists(paste0(old.results, "kff_results.csv"))){
      kff_results <- read.csv(paste0(old.results, "kff_results.csv"), check.names = F)
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
          
          # Mira values
          mira_values <- as.numeric(normalized.mira.counts[i, order(normalized.mira.counts[i,])])

          
          if(length(kff_neg_samples) != 0){
            # Converting these back into MIRA points
            mira_kff_neg_points <- as.numeric(normalized.mira.counts[i, kff_neg_samples])
            
            # X axis values
            x_vals <- match(mira_kff_neg_points, mira_values)
            
            # Graphing
            par(new = T)
            plot(x_vals, mira_kff_neg_points, ylab = NA, xlab = NA, axes = F, col = "red", pch = 4, xlim = c(1, x_lim), ylim = c(0, y_lim))
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
        
        legend("topleft", as.character(thresholds_legend), title="Thresholds")
          
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
  
  #old.results <- mira_count_table[nrow(mira_count_table), 2]
  old.results <- kff.results.directory
  mira_count_table <- mira_count_table[1: nrow(mira_count_table), ]
  
  row.names(mira_count_table) <- mira_count_table[,1]
  mira_count_table <- mira_count_table[,2:length(mira_count_table[1,]), drop=FALSE]
  mira_count_table[,1] <- as.numeric(as.character(mira_count_table[,1]))
    
  sample_list <- colnames(mira_count_table)
  
  # Running MIRA and reading in results -------------------------------------
    
    
  cat("\n\nNormalizing MIRA reads to 3DL3\n")
  normalized_mira_table <- normalize_mira_reads(sample_list, mira_count_table)
    
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
