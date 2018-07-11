# Format Caller Resources
#  - Locus alleles FASTA files and BED position conversion tables (convert between CDS and Genomic positions)
update_caller_resources <- function(
  ping.dir.path = "./PING/",
  supported.loci = c("2DL1", "2DL23", "2DL4", "2DL5", "2DS3", "2DS4", "2DS5", "2DP1", "3DL1", "3DS1", "3DL2", "3DL3"),
  results.directory = "./Formatted_caller_resources",
  resource_location = ""
) {
  
  ping.dir.path <- normalizePath(ping.dir.path, mustWork = T)
  setwd(paste0(ping.dir.path,"/"))
  library(data.table) ## Used for fread
  library(ape) ## Used for read.dna an seg.sites
  library(stringr)
  # Source HAPLOcaller functions for:
  #  - msf_to_allele_frame, is_nuc, num_unique_nuc, iterate 
  source(paste0(ping.dir.path,"/Resources/haplo_functions.R"), local = TRUE)
  
  
  ###################
  # INPUT FORMATS ----------------------------------------------------
  ###################
  
  # Create Results Directory if it does not exist
  dir.create(results.directory, showWarnings = F)
  cat(paste0(">> RESULTS DIR:\t", results.directory,"\n"))
  
  # Format all PATHs the same
  results.directory <- normalizePath(results.directory, mustWork = T)
  resource_location <- normalizePath(resource_location, mustWork = T)
  
  
  ###################
  # FUNCTIONS --------------------------------------------------------
  ###################
  
  # Process all unknown sequnce in the MSF data frame to determine whether it should be filled in or not
  #  - 'sep' represents the char that defines unknown sequence, by default it is '*'
  process_unknown_seq_in_allele_frame <- function(allele_frame,ref_allele,sep='*'){
    column_vec <- names(allele_frame)
    row_vec    <- row.names(allele_frame)
    
    # Create data frame of number of unique nucleotides per position
    nuc_unique_frame        <- data.frame(matrix(nrow=1, ncol=length(column_vec)), row.names = "Unique")
    names(nuc_unique_frame) <- column_vec
    for(pos in column_vec){
      nuc_vec                <- allele_frame[,pos]
      num_unique             <- num_unique_nuc(nuc_vec)
      nuc_unique_frame[,pos] <- num_unique
    }
    
    # Determine whether to fill in the gaps for alleles with unknown sequence
    for(allele in row_vec){
      allele_seq_vec <- allele_frame[allele,]
      if(any(allele_seq_vec == sep)){
        inv_allele_seq_vec     <- as.data.frame(t(allele_seq_vec))
        unknown_seq_pos        <- names(inv_allele_seq_vec[inv_allele_seq_vec[,1] == sep,])
        unknown_seq_unique_vec <- nuc_unique_frame[,unknown_seq_pos]
        
        # If any of the positions with unknown sequence have variability (more than one unique nucleotide)
        # do not replace it with sequence from the reference allele
        replace_unknown_seq     <- TRUE
        if(any(unknown_seq_unique_vec > 1)){
          replace_unknown_seq <- FALSE
        }
        
        if(replace_unknown_seq){
          allele_seq_vec[,unknown_seq_pos] <- ref_allele[,unknown_seq_pos] 
        }
        
        # Add Changes (or no changes) back into the allele dataframe
        allele_frame[allele,] <- allele_seq_vec
        
      }
    }
    
    return(allele_frame)
  }
  
  # Format allele string to fit the output format for PING's caller resource files
  format_allele_string <- function(allele, formatLength = TRUE){
    name_comp <- tstrsplit(allele,'*',fixed=T)
    allele_id <- name_comp[2]
    
    # If length of the alleles numeric ID is greater than 5, trim of the excess numbers (code for noncoding variation)
    if (formatLength){
      if (nchar(allele_id) > 5) {name_comp[2] <- substr(allele_id,1,5)}
    }
    
    allele_formatted <- paste0(">",paste(name_comp,collapse = '_'))
    
    return(allele_formatted)
  }
  
  
  ###################
  # MAIN CODE: -------------------------------------------------------
  ###################
  
  
  # Alleles used as reference for IPD-KIR multiple alignments for each KIR locus
  ref_allele_names <- c(
    "KIR2DL1*001",
    "KIR2DL2*0010101",
    "KIR2DL3*0010101",
    "KIR2DL4*00101",
    "KIR2DL5A*0010101",
    "KIR2DP1*0010201",
    "KIR2DS1*001",
    "KIR2DS2*0010101",
    "KIR2DS3*00101",
    "KIR2DS4*0010101",
    "KIR2DS5*001",
    "KIR3DL1*0010101",
    "KIR3DL2*0010101",
    "KIR3DL3*00101",
    "KIR3DP1*001",
    "KIR3DS1*010"
  )
  
  ## Iterate through each locus MSF file
  msf.files <- dir(resource_location, pattern = ".msf")
  
  cat(paste0("\n>> CREATING FORMATTED ALLELE RESOURCES:\n"))
  for(file in msf.files){
    current_locus <- tstrsplit(file,"_", fixed=T)[[1]]
    current_locus <- tstrsplit(current_locus,"KIR",fixed=T)[[2]]
    output_fasta_fh <- paste0(results.directory, "/All_", current_locus, "_preKFF.fas")
    
    cat(paste0(current_locus,"\t",output_fasta_fh,"\n"))
    
    # Generate Data frame from MSF resource file
    allele_frame <- msf_to_allele_frame(resource_location,current_locus)
    ref_allele   <- allele_frame[row.names(allele_frame) %in% ref_allele_names[grep(current_locus,ref_allele_names)],]
    allele_frame.processed <- process_unknown_seq_in_allele_frame(allele_frame,ref_allele)
    
    # Format MSF as FASTA entries 
    nucleotides <- c('A','T','G','C')
    row_names_vec <- row.names(allele_frame.processed)
    fasta_entries <- vector(mode="character", length = length(row_names_vec))
    
    # Format all non nucleotide BPs and fill in fasta_entries
    for(i in 1:length(row_names_vec)){
      allele <- row_names_vec[i]
      non_nuc_seq <- setdiff(allele_frame.processed[allele,],nucleotides)
      # Replace special characters ('*','.', etc.) with '-' for all positions in the Allele data frame
      if(length(non_nuc_seq) > 0){
        non_nuc_seq[1:length(non_nuc_seq)] <- '-'
        allele_frame.processed[allele,names(non_nuc_seq)] <- non_nuc_seq
      }
      
      # Collapse all positions to form full CDS sequence
      # - Make an exception to formatting the allele name for for KIR2DL4*00801
      fasta_allele <- allele
      if (length(grep("KIR2DL4*00801",fasta_allele,fixed = T)) != 0){
        fasta_allele <- format_allele_string(allele, formatLength = F)
      } else{
        fasta_allele <- format_allele_string(allele)
      }
      
      fasta_seq    <- paste(allele_frame.processed[allele,],collapse = '')
      entry        <- paste0(fasta_allele,"\t",fasta_seq)
      fasta_entries[i] <- entry
      
    }
    
    
    ## Make an exception for KIR2DL4*0080101/02/03/04/05
    if (current_locus == "2DL4"){
      subset_entries      <- c()
      subset_uniq_alleles <- c()
      # Isolate FASTA entries for KIR2DL4*00801
      for(i in 1:length(fasta_entries)){
        entry      <- unlist(fasta_entries[[i]])
        allele_str <- strsplit(entry,"\t")[[1]][1]
        allele_seq <- strsplit(entry,"\t")[[1]][2]
        if (length(grep("KIR2DL4_00801",allele_str,fixed = T)) != 0){
          subset_entries <- c(subset_entries, entry)
          subset_uniq_alleles <- c(subset_uniq_alleles,allele_seq)
        }
      }
      
      # Collapse unique allele variants of KIR2DL4*00801
      unique_allele00801_seqs <- as.matrix(unique(as.data.frame(subset_uniq_alleles)))
      alleles_to_keep         <- c()
      for (aseq in unique_allele00801_seqs){
        allele_str <- ""
        for(i in 1:length(subset_entries)){
          entry      <- unlist(subset_entries[[i]])
          allele_str <- strsplit(entry,"\t")[[1]][1]
          allele_seq <- strsplit(entry,"\t")[[1]][2]
          if (aseq == allele_seq){break}
        }
        alleles_to_keep <- c(alleles_to_keep,allele_str)
      }
      
      # Remove redundant KIR2DL4*00801 entries from list of FASTA allele entries
      fasta_entries_filtered = c()
      for (i in 1:length(fasta_entries)){
        entry      <- unlist(fasta_entries[[i]])
        allele_str <- strsplit(entry,"\t")[[1]][1]
        allele_seq <- strsplit(entry,"\t")[[1]][2]
        if (length(grep("KIR2DL4_00801",allele_str,fixed = T)) != 0){
          if (allele_str %in% alleles_to_keep){
            fasta_entries_filtered <- c(fasta_entries_filtered,entry)
          }
        }else{
          fasta_entries_filtered <- c(fasta_entries_filtered,entry)
        }
      }
      
      fasta_entries <- fasta_entries_filtered
    }
    
    ## Collapse all allele entries so only unique entries are left
    fasta_entries <- unique(as.data.frame(fasta_entries))
    ## Output Updated resource FASTA to Resource Directory
    write.table(fasta_entries,file = output_fasta_fh,quote = F,row.names = F,col.names = F)
    
  }
  
  ## Create Resource Files for 2DL23 and 2DS35
  #  - Concatenate the 2DL2 and 2DL3 resource files
  output_fasta_2DL2_fh  <- paste0(results.directory, "/All_2DL2_preKFF.fas")
  output_fasta_2DL3_fh  <- paste0(results.directory, "/All_2DL3_preKFF.fas")
  output_fasta_2DL23_fh <- paste0(results.directory, "/All_2DL23_preKFF.fas")
  
  system2("cat", c(output_fasta_2DL2_fh, ">",  output_fasta_2DL23_fh))
  system2("cat", c(output_fasta_2DL3_fh, ">>", output_fasta_2DL23_fh))
  
  #  - Concatenate the 2DS3 and 2DS5 resource files
  output_fasta_2DS3_fh  <- paste0(results.directory, "/All_2DS3_preKFF.fas")
  output_fasta_2DS5_fh  <- paste0(results.directory, "/All_2DS5_preKFF.fas")
  output_fasta_2DS35_fh <- paste0(results.directory, "/All_2DS35_preKFF.fas")
  
  system2("cat", c(output_fasta_2DS3_fh, ">",  output_fasta_2DS35_fh))
  system2("cat", c(output_fasta_2DS5_fh, ">>", output_fasta_2DS35_fh))
  
  cat(paste0("\n>> FINISHED <<\n"))
}


