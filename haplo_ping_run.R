library(gtools)
library(data.table)
library(stringr)
source('Resources/haplo_functions.R')

ping_haplo <- function(
  sample.location = '',
  fastq.pattern.1 = '_1.fastq.gz',
  fastq.pattern.2 = '_2.fastq.gz',
  bowtie.threads  = 4,
  results.directory = '',
  combined.csv.file = ''
){
  ipdkir_allele_df <- get_ipdkir_allele_df(kir_gen_path)
  ipdkir_nuc_df <- get_ipdkir_allele_df(kir_nuc_path)
  
  results.directory = file.path(results.directory)
  combined.csv.file = normalizePath(combined.csv.file, mustWork=T)

  master_haplo_path <- normalizePath(file.path(haplo_resources_directory, 'master_haplo_iteration_testing_v4.csv'), mustWork=T)
  
  ## Making sure we got the right output from ping_kff
  kff_results_path = normalizePath(file.path(results.directory, 'raw_kff_counts.csv'), mustWork=T)
  
  cat('\nCreating gc_input.csv\n')
  ## figuring out 2DL2/3 copy number, creates gc_input.csv file
  prepare_gc_input(kff_results_path,combined.csv.file,results.directory)
  
  ## Looking for the prepare_gc_input output file
  gc_input_path = normalizePath(file.path(results.directory,'gc_input.csv'), mustWork=T)
  
  cat('\nReading in gc_input.csv\n')
  ## Read in thegc_input.csv file
  master_gc_table <- read_master_gc(gc_input_path)
  
  ##  ADD WAY TO INPUT MASTER GC TABLE FOR ALL SAMPLES ##
  
  
  
  cat('\nReading in allele references\n')
  ## Reading in the reference allele name table
  master_haplo_table <- read_master_haplo(master_haplo_path)
  
  cat('\nReading in blacklisted alleles\n')
  ## Read in the blacklisted alleles (will not be called)
  allele_blacklist_table <- read_blacklist(file.path(caller_resources_directory, 'allele_blacklist.csv'))
  
  ## Normalizing loci names between the GC input and the reference table
  master_loci_intersection <- intersect(rownames(master_haplo_table), colnames(master_gc_table))
  
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
  
  sample_list = get_sample_list()
  
  ## Intersecting the samples found in sample.location with the GC results
  sample_list <- intersect(sample_list, rownames(master_gc_table))
  
  ## Cutting down the GC table to only include samples that will run
  gc_frame <- master_gc_table[sample_list,master_loci_intersection]
  
  ## Initializing data frames to store results
  type_frame <- data.frame(matrix(nrow=length(sample_list),ncol=length(master_loci_intersection)),row.names=sample_list)
  colnames(type_frame) <- master_loci_intersection
  distance_frame <- type_frame
  
  ## Initialize list for storing what samples are skipped
  skipped_samples <- c()
  for(sample_id in sample_list){
    cat('\n\n')
    cat(sample_id, '\n')
    
    #distance_frame <- read.table(file.path(results.directory, 'distance_frame.csv'), stringsAsFactors = F, sep = ',',check.names = F)
    #type_frame <- read.table(file.path(results.directory, 'type_frame.csv'), stringsAsFactors = F, sep = ',',check.names = F)
    #gc_frame <- read.table(file.path(results.directory, 'gc_frame.csv'),stringsAsFactors = F,sep=',',check.names=F)
    
    ## These two lines are unnecessary unless reading in the tables from file (previous three commented lines)
    distance_frame[sample_id,] <- NA
    type_frame[sample_id,] <- NA
    
    ## Getting the gc results for this sample
    sample_haplo_raw <- gc_frame[sample_id,]
    
    ## Figuring out what loci to build a reference for
    sample_haplo_loci <- names(sample_haplo_raw[,sample_haplo_raw > 0,drop=F])
    
    ## Going over 5 iterations of haplotype alignment
    for(i in 1:5){
      haplo <- master_haplo_table[sample_haplo_loci,i]
      ping_haplo_aligner(sample.location=sample.location,fastq.pattern.1=fastq.pattern.1,fastq.pattern.2=fastq.pattern.2,
                         results.directory=results.directory,bowtie.threads = bowtie.threads,sample.haplotype=haplo,sample.name=sample_id,
                         ipdkir_allele_df = ipdkir_allele_df, ipdkir_nuc_df = ipdkir_nuc_df, haplo.iteration=i)
    }
      
    result_list <- ping_haplo_caller(sample.name=sample_id,results.directory=results.directory,
                                      ipdkir_allele_df = ipdkir_allele_df, ipdkir_nuc_df = ipdkir_nuc_df)
      
    if(all(result_list == 'DID NOT PASS')){
      cat('\n\nThis sample did not align correctly. Moving on.\n\n')
      skipped_samples <- c(skipped_samples, sample_id)
      next
    }
      
    cat('\n\nMoving to genotype calling.\n')
    for(current_locus in sample_haplo_loci){
      copy_number <- sample_haplo_raw[,current_locus]
      
      ## Depending on experimental validation results, for certain loci we
      ## may add conditions to keep the KFF result over the MIRA result or vice vera
      if(is.character(copy_number)){ ## Hack for conflicting gc
        copy_number = 1
      }
        
      allele_caller_input <- allele_caller_input_formatting(current_locus, result_list, msf_directory, allele_blacklist_table)
      allele_type_list <- allele_caller(allele_caller_input$allele_calling_frame_list, allele_caller_input$possible_allele_frame, copy_number)
      
      ### -------- NEW ALLELES CAN GO ANYWHERE BELOW HERE --------- ###
      ### Variable position frame: allele_caller_input$possible_allele_frame
      ### Vcf position frame: allele_caller_input$
      
      if(all(allele_type_list == 'DID NOT PASS')){
        cat('\n\nThis sample did not have the correct GC input. Moving on.\n\n')
        skipped_samples <- c(skipped_samples, sample_id)
        type_frame[sample_id,] <- "nocall"
        distance_frame[sample_id,] <- "nocall"
        ## Add a 'nocall' flag into genotype output
        break
      }
      
      type_frame[sample_id,current_locus] <- paste0(allele_type_list$allele_names, collapse=' ')
      write.table(type_frame,file=file.path(results.directory,'type_frame.csv'),quote=F,sep=',')
      
      distance_frame[sample_id,current_locus] <- paste0(allele_type_list$distance, collapse=' ')
      write.table(distance_frame,file=file.path(results.directory,'distance_frame.csv'),quote=F,sep=',')
      cat('\n')
    }
    cat('\n\n-------')
  }
  
  cat('\n\nThese samples were skipped:\n')
  cat(paste0(skipped_samples, collapse='\n'))
}
