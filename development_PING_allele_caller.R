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


ping_allele_caller_dev <- function(
  sample.location = "PING_sequences/",
  fastq.pattern.1 = "_1.fastq.gz",
  fastq.pattern.2 = "_2.fastq.gz",
  bowtie.threads = 4,
  supported.loci = c("2DL1", "2DL23", "2DL4", "2DL5", "2DS3", "2DS4", "2DS5", "2DP1", "3DL1", "3DS1", "3DL2", "3DL3"),
  ping.gc.output = "Combined_results.csv",
  results.directory = "",
  resource_location = 'Resources/caller_resources/',
  allele_resources  = 'Resources/caller_resources/Formatted_caller_resources/',
  DPthresh          = 6
){
  
  library(data.table) ## Used for fread
  library(ape) ## Used for read.dna an seg.sites
  library(stringr)
  library(digest) ## use for md5 hash algorithm
  source("Resources/locus_functions.R", local = TRUE)
  source("Resources/haplo_functions.R", local = TRUE)
  source("Resources/caller_functions.R",local = TRUE)
  
  
  
  ###################
  # MAIN CODE: Allele caller -----------------------------------------
  ###################
  
  # EXCLUDE this approach for new ------------------------------------
  # INITIALIZE empty Data Frame to stor all New allele Variants
  #  - This needs to be done for ALL KIR loci, since one sample will have genos results for multiple KIR loci
  #newVariantsDataFrame.list <<- initialize_new_alleles_df(allele_resources)
  
  ### TEMP test code for adding new variant --------------------------
  # New variant Data Frame (first list element is 2DL1)
  #  - SOS_locus_lookup must already be generated for 2DL1 for this test code to work!
  #locus_from_list_frame     <- newVariantsDataFrame.list[[1]][[1]] 
  #newVariantsDataFrame.2DL1 <- newVariantsDataFrame.list[[1]][[2]] 
  
  #new_variant_df <- data.frame(matrix(nrow=1, ncol=length(names(SOS_locus_lookup))))
  #names(new_variant_df) <- names(SOS_locus_lookup)
  #new_variant_df$allele <- NULL
  #test_id_iter          <- 1
  #test_new_id           <- sprintf("new%04d",test_id_iter)
  #row.names(new_variant_df) <- test_new_id
  #test_vec_new              <- c("A","T","G")
  #test_pos_new              <- c("X1279","X4788","X4866")
  #new_variant_df[,test_pos_new] <- test_vec_new
  
  #newVariantsDataFrame.2DL1           <- rbind(newVariantsDataFrame.2DL1,new_variant_df)
  #newVariantsDataFrame.list[[1]][[2]] <- newVariantsDataFrame.2DL1
  
  
  # START PING -------------------------------------------------------
  
  # Create a results directory
  results.directory <- results_directory(results.directory)
  
  # Create results subfolders
  dir.create(results.directory, showWarnings = F)
  dir.create(paste0(results.directory, "Vcf"), showWarnings = F)
  dir.create(paste0(results.directory, "Vcf_Genomic"), showWarnings = F)
  dir.create(paste0(results.directory, "KIRcaller"), showWarnings = F)
  dir.create(paste0(results.directory, "Fastq"), showWarnings = F)
  cat("Results directories created.\n\n")
  
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
  
  # Pull in the GC results (original PING function)
  #gc_results <- gc.results()
  
  ## Read in Custom GC table for all samples
  gc_results <- read.csv(ping.gc.output,header = TRUE, row.names = 1, check.names = F)
  gc_results <- as.data.frame(t(gc_results))
  # Line needed to get results from Newest GC caller to work
  row.names(gc_results) <- tstrsplit(row.names(gc_results),"KIR")[[2]]
  
  # The master list is the intersection of sample name matches in sequence_list and gc_results
  master_list <- master.creator(sequence_list, gc_results)
  
  # Initialize Vector of reference KIR alleles for formatting new allele names
  #  - Alleles used as reference for IPD-KIR multiple alignments for each KIR locus as of June 2018
  ref_allele_names <- c(
    "KIR2DL1_001",
    "KIR2DL2_00101",
    "KIR2DL3_00101",
    "KIR2DL4_00101",
    "KIR2DL5A_00101",
    "KIR2DP1_00102",
    "KIR2DS1_001",
    "KIR2DS2_00101",
    "KIR2DS3_00101",
    "KIR2DS4_00101",
    "KIR2DS5_001",
    "KIR3DL1_00101",
    "KIR3DL2_00101",
    "KIR3DL3_00101",
    "KIR3DP1_001",
    "KIR3DS1_010"
  )
  
  # Output the sample list used by the allele caller
  write.table(as.data.frame(master_list), file = paste0(results.directory,"sample_list.txt"), row.names = F, col.names = F,quote = F)
  
  ### START: Generating PING results per sample
  for(sample in master_list){
    # Checking to see what loci are present for each sample based on PING_gc results
    locus_presence <- gc_results[,pmatch(sample, colnames(gc_results))]
    
    # As long as the PING_gc locus result isn't 0, it will be run, this includes discrepencies
    loci.list <- rownames(gc_results)[locus_presence != 0]
    
    # Cutting down the loci.list to only include supported loci
    loci.list <- loci.list[loci.list %in% supported.loci]
    
    # Mutating 2DS5 to 2DS35 for allele calling
    if("2DS5" %in% loci.list){loci.list[loci.list == "2DS5"] <- "2DS35"}
    
    # Going backwards from master_list sample name to sequence_list sample name
    sample <- sequence_list[pmatch(sample, sequence_list)]
    
    for(current.locus in loci.list){
      if(current.locus == "2DS3" && any("2DS35" %in% loci.list)){next}
      
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
      if (current.locus == "2DL23"){
        for (locus in c("2DL2","2DL3","2DL23")){
          alleles_gen.vcf(locus, resource.location = allele_resources)  # Re-formatted Resource files
          snps_gen.vcf(locus)  # Generate necessary files to create a lookup table for 2DL23,2DL3 or 2DL2
        }
      }else{
        alleles_gen.vcf(current.locus, resource.location = allele_resources)  # Re-formatted Resource files
      }
      
      
      # Generating position information for the current locus
      SOS_locus_lookup <- snps_gen.vcf(current.locus)
      positions        <- gsub("X", "", names(SOS_locus_lookup)[-1])
      
      # Running the bowtie2 scripts
      ping.locus(sample, current.locus, is_gz)
      
      # 2DL23 specific logic
      if("2DL23" == current.locus){
        ## Determine if 2DL2+ 2DL3+
        # 2DL2+ test
        pos2DL2 <- FALSE
        vcf_table <- tryCatch(read.table(paste0(results.directory, "Vcf/", sample, "_2DL2nuc.vcf"), header = F), error=function(e) NULL)
        
        if(!is.null(vcf_table)){
          vcf_table <- nodels.vcf(vcf_table)
        }
        
        if(length(vcf_table[,1]) == 383) {
          pos2DL2 <- TRUE
          alleles_gen.vcf("2DL2", resource.location = allele_resources)
        }
        
        # 2DL3+ test
        pos2DL3 <- FALSE
        vcf_table <- tryCatch(read.table(paste0(results.directory, "Vcf/", sample, "_2DL3nuc.vcf"), header = F), error=function(e) NULL)
        
        if(!is.null(vcf_table)){
          vcf_table <- nodels.vcf(vcf_table)
        }
        
        if(length(vcf_table[,1]) == 386) {
          pos2DL3 <- TRUE
        }
        
        
        presence2DL2or3 <- list("two" = pos2DL2, "three" = pos2DL3)
      }
      
      ## ALLELE CALLING: Running the allele calling scripts
      has_genotype <- tryCatch(ping.caller(sample, current.locus,DPthresh), error=function(e) NULL)
      
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

