## Generate a haplotype bed file (Non-coding) -------------------
make_bed_full <- function(haplo_directory, allele_list, bed_file_path, ipdkir_gen_df){
  ## 2DL3 and 3DS1 were guestimates, they might cause problems (11/16/18 WMM)
  # These are hardcoded exon lengths
  
  first=T
  for(allele in allele_list){
    current_locus <- unlist(strsplit(allele, '*', fixed=T))[1]
    current_locus <- unlist(strsplit(current_locus, 'KIR', fixed=T))[2]
    if(!(current_locus %in% names(exon_lengths))){
      stop(paste0(allele, ': locus not found in exon_lengths'))
    }
    gen_chr_string <- ipdkir_gen_df[allele,]
    current_locus_allele_frame <- msf_to_allele_frame(haplo_directory, current_locus)
    
    exon_list <- c()
    start <- 0
    end <- nchar(gen_chr_string)
    # for(exon_length in exon_lengths[[current_locus]]){
    #   start <- total_length
    #   end <- (exon_length+total_length)-1
    #   current_exon <- current_locus_allele_frame[allele,start:end]
    #   current_exon <- paste0(current_exon[,names(which(lapply(current_exon, is_nuc) == TRUE))], collapse='')
    #   exon_list <- c(exon_list, current_exon)
    #   total_length <- total_length + exon_length
    #   
    #   if(str_count(gen_chr_string, current_exon) != 1){
    #     cat('\n\nDuplication! ', allele, ' ', nchar(current_exon), ' ', current_exon)
    #     
    #     known_duplications <- c('KIR3DL3*00101', 'KIR3DL3*00201', 'KIR3DL3*0030101', 'KIR3DS1*0130101', 'KIR3DL3*0090101', 'KIR3DL3*00801', 'KIR3DL3*007',
    #                             'KIR3DL3*00202', 'KIR3DL3*00206', 'KIR3DL3*00601', 'KIR3DS1*055', 'KIR3DS1*078', 'KIR3DL3*00802', 'KIR3DL3*036',
    #                             'KIR3DS1*0130108', 'KIR3DL3*005')
    #     
    #     if(allele %in% known_duplications & grepl(current_exon, gen_chr_string)){
    #       cat('\nThis is a known duplication')
    #     }else{
    #       stop()
    #     }
    #   }
    # }
    
    #exon_position_frame <- data.frame(matrix(nrow=length(exon_lengths[[current_locus]]), ncol=3))
    
    allele_position_frame <- data.frame(matrix(nrow = 1,ncol = 3))
    allele_position_frame[1,] <- c(allele, start,end)
    #attribute_list <- lapply(exon_list, regexpr, gen_chr_string)
    #allele_name <- sub(current_locus, paste0('KIR', current_locus), allele)
    #allele_name <- sub('*', '_', allele_name, fixed=T)
    
    ## This section matches exon strings to their position in the allele. In the case of a duplication, the first match after the previous exon is taken.
    # gen_chr_substr <- gen_chr_string
    # substr_start <- 1
    # for(i in 1:length(exon_list)){
    #   exon <- exon_list[i]
    # 
    #   if(nchar(exon) == 0){
    #     next()
    #   }
    #   exon_attributes <- regexpr(exon, gen_chr_substr)
    #   exon_start <- exon_attributes[1] + substr_start-1
    #   exon_end <- exon_start + as.integer(attr(exon_attributes, 'match.length'))-1
    # 
    #   substr_start <- exon_end
    #   gen_chr_substr <- substr(gen_chr_string, substr_start, nchar(gen_chr_string))
    #   exon_position_frame[i,] <- c(allele, exon_start, exon_end)
    # }
    # exon_position_frame <- exon_position_frame[!is.na(exon_position_frame[,1]),,drop=F]
    # exon_position_frame[,'X2'] <- as.integer(exon_position_frame[,'X2']) - 1

    write.table(allele_position_frame, bed_file_path, quote=F, sep=' ', row.names=F, col.names=F, append=!first)
    first=F
  }
}

## --------------------------------------------------------------



moved_sample_list <- c('IND00008','IND00093','IND00012','IND00014','IND00020','IND00033','IND00045',
                       'IND00057','IND00069','IND00081','IND00094','IND00026','IND00194')

moved_sample_conversion <- list('IND00012'='IND00014',
                                'IND00014'='IND00012',
                                'IND00020'='IND00008',
                                'IND00033'='IND00020',
                                'IND00045'='IND00033',
                                'IND00057'='IND00045',
                                'IND00069'='IND00057',
                                'IND00081'='IND00069',
                                'IND00094'='IND00081')

## 2DL3 and 3DS1 were guestimates, they might cause problems (11/16/18 WMM)
# These are hardcoded exon lengths
exon_lengths <- list('2DL1'=c(34, 36, 300, 294, 51, 102, 53, 177),
                     '2DL2'=c(34, 36, 300, 294, 51, 102, 53, 177),
                     '2DL3'=c(34, 36, 300, 294, 51, 105, 51, 155),
                     '2DL4'=c(40, 36, 285, 294, 51, 105, 53, 270),
                     '2DL5A'=c(34, 36, 285, 294, 51, 105, 53, 270),
                     '2DL5B'=c(34, 36, 285, 294, 51, 105, 53, 270),
                     '2DS1'=c(34, 36, 300, 294, 51, 105, 53, 42),
                     '2DS2'=c(34, 36, 300, 294, 51, 105, 53, 42),
                     '2DS3'=c(34, 36, 300, 294, 51, 105, 53, 42),
                     '2DS4'=c(34, 36, 300, 294, 51, 105, 53, 42),
                     '2DS5'=c(34, 36, 300, 294, 51, 105, 53, 42),
                     '3DL1'=c(34, 36, 285, 300, 294, 51, 105, 53, 177),
                     '3DL2'=c(34, 36, 285, 300, 294, 51, 105, 53, 210),
                     '3DL3'=c(34, 36, 285, 300, 294, 105, 53, 126),
                     '3DS1'=c(34, 36, 285, 300, 294, 51, 105, 52, 8),
                     '2DP1'=c(34, 36, 299, 294, 51, 105, 53, 177))



### Creating IPDKIR Reference sequence map: allele sequence data fram from IPDKIR Reference file (CORRECT approach)
get_ipdkir_allele_df <- function(REF_IPDKIR) {
  # Reference FASTA = KIR_gen.fasta available on IPDKIR
  ipdkir_ref <- read.table(file = REF_IPDKIR, sep = '\n', header = F)
  
  ref_fasta_ipdkir_df <- data.frame()
  allele_iter <- 1
  
  allele_id  <- ipdkir_ref[1,]
  allele_seq <- ""
  for(i in 2:length(row.names(ipdkir_ref))) {
    header_status <- ifelse(grepl(">", ipdkir_ref[i,]), TRUE, FALSE)
    if(header_status == TRUE) {
      allele_id <- as.character(allele_id)
      id_pieces <- unlist(strsplit(allele_id, "\\s+"))
      allele_name <- id_pieces[2]
      allele_name <- gsub("KIR", "", allele_name, fixed = T)
      #allele_name_mod <- gsub("*", "_", allele_name, fixed=T)
      formatted_allele_name <- paste("KIR", allele_name, sep = "")
      #cat(paste0(allele_iter, ">->->", allele_id, "\n", formatted_allele_name, "\n"))
      ref_fasta_ipdkir_df[allele_iter,1] <- formatted_allele_name
      ref_fasta_ipdkir_df[allele_iter,2] <- allele_seq 
      
      allele_iter <- allele_iter + 1
      allele_id <- ipdkir_ref[i,]
      allele_seq <- ""
    }else{
      allele_seq <- paste(allele_seq, ipdkir_ref[i,], sep = "")
    }
  }
  # last Allele and Seq
  last_allele_id <- as.character(allele_id)
  id_pieces <- unlist(strsplit(last_allele_id, "\\s+"))
  allele_name <- id_pieces[2]
  allele_name <- gsub("KIR", "", allele_name, fixed = T)
  #allele_name_mod <- gsub("*", "_", allele_name, fixed=T)
  formatted_allele_name <- paste("KIR", allele_name, sep = "")
  ref_fasta_ipdkir_df[allele_iter,1] <- formatted_allele_name
  ref_fasta_ipdkir_df[allele_iter,2] <- allele_seq 
  
  # Final Formats for reference map
  names(ref_fasta_ipdkir_df) <- c("allele", "seq")
  row.names(ref_fasta_ipdkir_df) <- ref_fasta_ipdkir_df$allele
  ref_fasta_ipdkir_df$allele <- NULL
  
  return(ref_fasta_ipdkir_df)
}

## Extract Data Frame of KIR alleles that are full length (not just coding sequence)
extract_full_alleles <- function(gen_df, nuc_df, thresh, locus) {
  gen <- gen_df
  nuc <- nuc_df
  locus <- paste0('KIR', locus)
  
  length_df <- data.frame()
  allele_vec = row.names(gen)
  for(i in 1:length(allele_vec)) {
    allele          <- allele_vec[i]
    full_allele_seq <- gen[allele,1]
    cds_allele_seq  <- nuc[allele,1]
    full_allele_len <- nchar(full_allele_seq)
    cds_allele_len  <- nchar(cds_allele_seq)
    length_df[i,1] <- allele
    length_df[i,2] <- full_allele_len
    length_df[i,3] <- cds_allele_len
    length_df[i,4] <- full_allele_len / cds_allele_len
  }
  names(length_df) <- c("allele", "full_length", "CDS_length", "norm_full_length")
  
  
  if(locus == 'KIR2DL5'){
    length_df <- subset(length_df, grepl(locus, length_df$allele, fixed=T))
  }else{
    length_df <- subset(length_df, grepl(paste0(locus, "*"), length_df$allele, fixed=T))
  }
  norm_length_list <- as.numeric(length_df$norm_full_length)
  max_norm_length <- max(norm_length_list)
  threshold <- thresh * max_norm_length
  
  
  length_df_subset <- subset(length_df, length_df$norm_full_length >= threshold)
  full_allele_vec <- length_df_subset$allele
  
  #full_allele_vec <- gsub('_', '*', good_alleles, fixed=T)
  
  return(full_allele_vec)
}

### Generates a Haplotype FASTA file
generate_haplo_ref <- function(allele_list, ipdkir_df, file_path) {
  
  # iterate through each allele in vector and assemble haplo FASTA file
  sink(file_path)
  
  for(i in 1:length(allele_list)) {
    allele <- allele_list[i]
    
    allele_ref_seq <- as.character(ipdkir_df[row.names(ipdkir_df) == allele,])
    if(allele != "neg"){
      if(length(allele_ref_seq) == 0){
        allele_pieces <- unlist(strsplit(allele, split = '_'))
        a_number <- allele_pieces[2]
        if (nchar(a_number) == 5) {
          a_number_new <- paste(a_number, "01", sep = "") # always choose first non-coding variant (for the time being)
          a_new <- paste(allele_pieces[1], a_number_new, sep = "_")
          allele_ref_seq <- as.character(ipdkir_df[row.names(ipdkir_df) == a_new,])
        }
      }
    }
    
    if(allele != "neg" && length(allele_ref_seq) > 0){cat(paste0(">", allele, "\n", allele_ref_seq, "\n"))}
    
  }
  
  sink()
}

## Generate a haplotype bed file
make_bed <- function(haplo_directory, allele_list, bed_file_path, ipdkir_gen_df){
  ## 2DL3 and 3DS1 were guestimates, they might cause problems (11/16/18 WMM)
  # These are hardcoded exon lengths
  exon_lengths <- list('2DL1'=c(34, 36, 300, 294, 51, 102, 53, 177),
                       '2DL2'=c(34, 36, 300, 294, 51, 102, 53, 177),
                       '2DL3'=c(34, 36, 300, 294, 51, 105, 51, 155),
                       '2DL4'=c(40, 36, 285, 294, 51, 105, 53, 270),
                       '2DL5A'=c(34, 36, 285, 294, 51, 105, 53, 270),
                       '2DL5B'=c(34, 36, 285, 294, 51, 105, 53, 270),
                       '2DS1'=c(34, 36, 300, 294, 51, 105, 53, 42),
                       '2DS2'=c(34, 36, 300, 294, 51, 105, 53, 42),
                       '2DS3'=c(34, 36, 300, 294, 51, 105, 53, 42),
                       '2DS4'=c(34, 36, 300, 294, 51, 105, 53, 42),
                       '2DS5'=c(34, 36, 300, 294, 51, 105, 53, 42),
                       '3DL1'=c(34, 36, 285, 300, 294, 51, 105, 53, 177),
                       '3DL2'=c(34, 36, 285, 300, 294, 51, 105, 53, 210),
                       '3DL3'=c(34, 36, 285, 300, 294, 105, 53, 126),
                       '3DS1'=c(34, 36, 285, 300, 294, 51, 105, 52, 8),
                       '2DP1'=c(34, 36, 299, 294, 51, 105, 53, 177))
  first=T
  for(allele in allele_list){
    current_locus <- unlist(strsplit(allele, '*', fixed=T))[1]
    current_locus <- unlist(strsplit(current_locus, 'KIR', fixed=T))[2]
    if(!(current_locus %in% names(exon_lengths))){
      stop(paste0(allele, ': locus not found in exon_lengths'))
    }
    gen_chr_string <- ipdkir_gen_df[allele,]
    current_locus_allele_frame <- msf_to_allele_frame(haplo_directory, current_locus)
    
    exon_list <- c()
    total_length <- 1
    for(exon_length in exon_lengths[[current_locus]]){
      start <- total_length
      end <- (exon_length+total_length)-1
      current_exon <- current_locus_allele_frame[allele,start:end]
      current_exon <- paste0(current_exon[,names(which(lapply(current_exon, is_nuc) == TRUE))], collapse='')
      exon_list <- c(exon_list, current_exon)
      total_length <- total_length + exon_length
      
      if(str_count(gen_chr_string, current_exon) != 1){
        cat('\n\nDuplication! ', allele, ' ', nchar(current_exon), ' ', current_exon)
        
        known_duplications <- c('KIR3DL3*00101', 'KIR3DL3*00201', 'KIR3DL3*0030101', 'KIR3DS1*0130101', 'KIR3DL3*0090101', 'KIR3DL3*00801', 'KIR3DL3*007',
                                'KIR3DL3*00202', 'KIR3DL3*00206', 'KIR3DL3*00601', 'KIR3DS1*055', 'KIR3DS1*078', 'KIR3DL3*00802', 'KIR3DL3*036',
                                'KIR3DS1*0130108', 'KIR3DL3*005')
        
        if(allele %in% known_duplications & grepl(current_exon, gen_chr_string)){
          cat('\nThis is a known duplication')
        }else{
          stop()
        }
      }
    }
    
    exon_position_frame <- data.frame(matrix(nrow=length(exon_lengths[[current_locus]]), ncol=3))
    
    #attribute_list <- lapply(exon_list, regexpr, gen_chr_string)
    #allele_name <- sub(current_locus, paste0('KIR', current_locus), allele)
    #allele_name <- sub('*', '_', allele_name, fixed=T)
    
    ## This section matches exon strings to their position in the allele. In the case of a duplication, the first match after the previous exon is taken.
    gen_chr_substr <- gen_chr_string
    substr_start <- 1
    for(i in 1:length(exon_list)){
      exon <- exon_list[i]
      
      if(nchar(exon) == 0){
        next()
      }
      exon_attributes <- regexpr(exon, gen_chr_substr)
      exon_start <- exon_attributes[1] + substr_start-1
      exon_end <- exon_start + as.integer(attr(exon_attributes, 'match.length'))-1
      
      substr_start <- exon_end
      gen_chr_substr <- substr(gen_chr_string, substr_start, nchar(gen_chr_string))
      exon_position_frame[i,] <- c(allele, exon_start, exon_end)
    }
    exon_position_frame <- exon_position_frame[!is.na(exon_position_frame[,1]),,drop=F]
    exon_position_frame[,'X2'] <- as.integer(exon_position_frame[,'X2']) - 1
    
    write.table(exon_position_frame, bed_file_path, quote=F, sep=' ', row.names=F, col.names=F, append=!first)
    first=F
  }
}

## Support function for make_bed, this function turns an allele into a character vector
gen_to_chr_vect <- function(haplo_directory, allele){
  gen_file_path <- file.path(haplo_directory, 'KIR_gen.fasta')
  
  gen_frame <- read.table(gen_file_path, sep='\n', stringsAsFactors = F)
  all_allele_start_pos <- grep('>', gen_frame[,1], fixed=T)
  allele_start_pos <- grep(paste0(allele,' '), gen_frame[,1], fixed=T)
  allele_end_pos <- all_allele_start_pos[grep(allele_start_pos, all_allele_start_pos)+1][1]
  allele_gen_frame <- gen_frame[(allele_start_pos+1):(allele_end_pos-1),]
  allele_gen_string <- paste0(allele_gen_frame, collapse='')
  
  allele_gen_chr_vect <- strsplit(allele_gen_string, '')[[1]]
  return(allele_gen_chr_vect)
}

## Support funciton for make_bed, this function turns a msf file into a dataframe of aligned alleles
msf_to_allele_frame <- function(haplo_directory, current_locus){
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

## Returns boolean
is_nuc <- function(chr){
  nuc_list <- c('A', 'T', 'C', 'G')
  return(as.character(chr) %in% nuc_list)
}

## Counts how many unique nucleotides are in a vector
num_unique_nuc <- function(vector_to_check){
  return(sum(is_nuc(names(table(vector_to_check)))))
}

## Gotta iterate
iterate <- function(first_value, last_value){
  return(as.integer(first_value):as.integer(last_value))
}

## Condenses an allele frame into only positions that have multiple nucleotides
allele_frame_to_unique_pos_frame <- function(allele_frame, ref_allele){
  unique_pos_frame <- make_unique_pos_frame(allele_frame)
  return(unique_pos_frame)
}

## Subsets an allele frame to only have positions from ref_allele
make_relative_allele_frame <- function(allele_frame, ref_allele, list_allele_pos_conv){
  allele_frame <- allele_frame[,is_nuc(allele_frame[ref_allele,])]
  
  if(length(list_allele_pos_conv) != ncol(allele_frame)){
    cat('\n\nBed allele and relative allele frame do not match up!')
    stop()
  }
  
  colnames(allele_frame) <- list_allele_pos_conv
  
  return(allele_frame)
}

## Condenses allele fram into only positions that have multiple nucleotides
make_unique_pos_frame <- function(pos_frame){
  num_unique <- lapply(pos_frame, num_unique_nuc)
  unique_pos_frame <- pos_frame[,names(num_unique[num_unique > 1]), drop=F]
  return(unique_pos_frame)
}

get_haplo_allele_list <- function(bed_file_path){
  bed_file_frame <- data.frame(read.table(bed_file_path, header=F, stringsAsFactors = F))
  bed_allele_list <- unique(bed_file_frame[,'V1'])
  return(bed_allele_list)
}

alternate_make_bed_to_pos_conv <- function(bed_file_path, allele_frame, ref_allele){
  bed_file_frame <- data.frame(read.table(bed_file_path, header=F, stringsAsFactors = F))
  bed_allele_list <- unique(bed_file_frame[,'V1'])
  list_allele_bed_to_pos <- list()
  
  ## Sanity checks
  if(!(ref_allele %in% bed_allele_list)){
    cat('\n\n', ref_allele, ' reference allele is not found in the bed file!')
    stop()
  }else if(!(ref_allele %in% rownames(allele_frame))){
    cat('\n\n', ref_allele, ' reference allele is not found in the allele frame!')
    stop()
  }
  total_sum_check <- 0
  
  bed_allele_frame <- bed_file_frame[bed_file_frame[,'V1'] == ref_allele,]
  allele_frame <- allele_frame[,is_nuc(allele_frame[ref_allele,])]
  allele_frame_col <- 1
  
  for(i in 1:nrow(bed_allele_frame)){
    first_value <- bed_allele_frame[i,2]+1
    last_value <- bed_allele_frame[i,3]
    iter_values <- iterate(first_value, last_value)
    total_sum_check <- total_sum_check + length(iter_values)

    for(iter in iter_values){
      list_allele_bed_to_pos[[paste0('X',iter)]] <- paste0(colnames(allele_frame)[allele_frame_col], 'X')
      allele_frame_col <- allele_frame_col + 1
    }
  }
  
  ## Another sanity check
  if(total_sum_check != ncol(allele_frame)){
    cat('\n\n', ref_allele, ' bed coordinates do not match up with the allele frame!')
    stop()
  }
  
  return(list_allele_bed_to_pos)
}

## Returns a list that converts between bed positions and allele positions
make_bed_to_pos_conv <- function(bed_file_path){
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

## Find distances between alleles
make_allele_distance_frame <- function(frame_allele_unique){
  frame_allele_distance <- data.frame(matrix(nrow=nrow(frame_allele_unique), ncol=nrow(frame_allele_unique)),
                                      row.names = rownames(frame_allele_unique))
  
  colnames(frame_allele_distance) <- rownames(frame_allele_distance)
  
  #base_num <- ncol(frame_allele_unique)
  
  for(row_allele in rownames(frame_allele_distance)){
    #base_num <- sum(is_nuc(frame_allele_unique[row_allele,]))
    
    for(col_allele in colnames(frame_allele_distance)){
      base_num <- min(c(sum(is_nuc(frame_allele_unique[row_allele,])), sum(is_nuc(frame_allele_unique[col_allele,]))))
      distance <- sum(unlist(lapply(X = frame_allele_unique[c(row_allele, col_allele),], FUN = num_unique_nuc)))
      distance <- distance - base_num
      frame_allele_distance[row_allele, col_allele] <- distance
    }
  }
  return(frame_allele_distance)
}

make_allele_distance_list_from_ref <- function(frame_allele_unique, ref_allele){
  list_allele_distance <- list()
  for(allele_name in rownames(frame_allele_unique)){
    base_num <- sum(is_nuc(frame_allele_unique[ref_allele,]))
    
    distance <- sum(unlist(lapply(frame_allele_unique[c(allele_name, ref_allele),], num_unique_nuc)))
    distance <- distance - base_num
    list_allele_distance[[allele_name]] <- distance
  }
  return(list_allele_distance)
}

## Returns closest alleles to ref_allele (excludes any with distance of 0)
find_closest_alleles_from_allele <- function(ref_allele, list_allele_distance){
  
  if(!(ref_allele %in% names(list_allele_distance))){
    cat('\n\n', ref_allele, ': Reference allele not found!')
    stop()
  }
  
  non_ref_list <- list_allele_distance[which(list_allele_distance > 0)]
  min_non_ref <- min(unlist(non_ref_list))
  allele_names <- names(which(list_allele_distance == min_non_ref))
  
  return(list(distance=min_non_ref, alleles=allele_names))
}

## Returns furthest alleles to ref_allele (excludes any with distance of 0)
find_furthest_alleles_from_allele <- function(ref_allele, list_allele_distance){
  
  if(!(ref_allele %in% names(list_allele_distance))){
    cat('\n\n', ref_allele, ': Reference allele not found!')
    stop()
  }
  
  non_ref_list <- list_allele_distance[which(list_allele_distance > 0)]
  max_non_ref <- max(unlist(non_ref_list))
  allele_names <- names(which(list_allele_distance == max_non_ref))
  
  return(list(distance=max_non_ref, alleles=allele_names))
}

## Returns furthest alleles to all alleles in ref_allele_list (needs to be greater than 0)
find_furthest_alleles_from_ref_list <- function(ref_allele_list, frame_allele_unique){
  frame_ref_allele <- frame_allele_unique[ref_allele_list,,drop=F]
  frame_distance <- data.frame(matrix(nrow=nrow(frame_allele_unique), ncol=1),
                              row.names = rownames(frame_allele_unique))
  colnames(frame_distance) <- 'distance'
  frame_distance[,'distance'] <- 0

  ## Any positions with . should be cut out here
  only_nucs_pos_list <- names(which(unlist(lapply(lapply(frame_ref_allele, is_nuc), all))))
  frame_ref_allele <- frame_ref_allele[,only_nucs_pos_list,drop=F]
  frame_allele_unique <- frame_allele_unique[,only_nucs_pos_list,drop=F]
  
  ## Counting mismatches
  for(snp_pos in only_nucs_pos_list){
    mismatch_rows <- !(frame_allele_unique[,snp_pos] %in% frame_ref_allele[,snp_pos])
    mismatch_rownames <- rownames(frame_allele_unique[mismatch_rows,])
  
    if(length(mismatch_rownames) > 0){
      frame_distance[mismatch_rownames,'distance'] <- frame_distance[mismatch_rownames,'distance'] + 1
    }
  }
  
  max_distance <- max(frame_distance[,'distance'])
  max_indices <- which(frame_distance[,'distance'] == max_distance)
  max_alleles <- rownames(frame_distance[max_indices,,drop=F])
  return(list(distance=max_distance, alleles=max_alleles))
}

## Returns closest alleles to snp_frame
find_allele_distances_from_snps <- function(snp_frame, frame_allele_unique_updated){
  # Takes in a SNP position dataframe with SNP info, matches it to the closest alleles in the
  # unique frame. Returns the names of the closest alleles
  
  cols_in_both <- intersect(colnames(frame_allele_unique_updated), colnames(snp_frame))
  frame_allele_unique_updated['SNP',cols_in_both] <- snp_frame['SNP',cols_in_both]
  cols_only_in_allele_frame <- setdiff(colnames(frame_allele_unique_updated), colnames(snp_frame))
  frame_allele_unique_updated['SNP',cols_only_in_allele_frame] <- '.'
  
  base_num <- ncol(frame_allele_unique_updated)
  distance_list <- list()
  for(row_allele in rownames(frame_allele_unique_updated)){
    distance <- sum(unlist(lapply(X=frame_allele_unique_updated[c('SNP', row_allele),,drop=F], num_unique_nuc)))
    distance <- distance - base_num
    distance_list[[row_allele]] <- distance
  }
  
  distance_list <- unlist(distance_list)
  distance_list <- distance_list[!(names(distance_list) == 'SNP')]
  distance_min <- min(distance_list)
  closest_alleles <- names(which(distance_list == distance_min))

  return(list(distance_list=distance_list, min_distance=distance_min, closest_alleles=closest_alleles))
}

## Vcf quality control
filter_vcf_frame <- function(vcf_frame, vcf_threshold){
  
  if(nrow(vcf_frame) > 0){ ## Empty vcf_frames should pass right through
    dp_list <- as.integer(tstrsplit(tstrsplit(vcf_frame$V8, ';')[[1]], 'DP=')[[2]])
    cut_off <- vcf_threshold*mean(dp_list)
    vcf_frame <- vcf_frame[((dp_list >= cut_off) & (dp_list >= 20)),,drop=F]
  }
  
  if(nrow(vcf_frame) > 0){
    ## Removing low quality positions
    vcf_frame <- vcf_frame[(vcf_frame$V10 != './.'),,drop=F]
  }
  
  if(nrow(vcf_frame) > 0){  
    ## Removing low DP4 positions
    dp4_thresh <- 3
    
    attribute_list <- strsplit(vcf_frame$V8, ';')
    dp4_pos <- unlist(lapply(attribute_list, FUN=(function(x) grep('DP4',x))))
    
    dp4_list <- lapply(1:length(attribute_list), FUN=(function(x) attribute_list[[x]][dp4_pos[x]]))
    dp4_list <- tstrsplit(dp4_list,'DP4=')[[2]]
    dp4_split_list <- strsplit(dp4_list, ',')
    
    dp4_cutoff <- unlist(lapply(dp4_split_list, (function(x) 
      ((as.integer(x[1])>dp4_thresh | as.integer(x[2])>dp4_thresh) || (as.integer(x[3])>dp4_thresh | as.integer(x[4])>dp4_thresh)))
      ))
    
    vcf_frame <- vcf_frame[dp4_cutoff,]
  }
    
    #cat('\n\nCutoff: ', cut_off, '\nMean: ', mean(dp_list), '\nSD: ', sd(dp_list))
  if(nrow(vcf_frame) < 1){
    cat('\nRemoved all rows of vcf_frame in filter_vcf_frame')
      #top()
  }
  return(vcf_frame)
}

mean_dp_vcf_frame <- function(vcf_frame){
  dp_list <- as.integer(tstrsplit(tstrsplit(vcf_frame$V8, ';')[[1]], 'DP=')[[2]])
  return(mean(dp_list))
}

## Call alleles from locus specific vcf result frames
allele_caller <- function(allele_calling_frame_list, possible_allele_frame, n_alleles){
  
  ## Saving a copy of the possible allele frame, since it will be cut down during this function
  original_possible_allele_frame <- possible_allele_frame
  variant_positions <- colnames(possible_allele_frame)
  
  ## Filtering any calls in the VCF files that do not pass the threshold
  allele_calling_frame_list_filtered <- lapply(allele_calling_frame_list, filter_vcf_frame, 0.1)
  
  ## Finding all positions that are in any VCF frame after filtration
  filtered_rownames <- lapply(allele_calling_frame_list_filtered, rownames)
  unique_filtered_rownames <- unique(unlist(filtered_rownames))
  
  ## If there are no filtered rownames, then say this locus is negative for this sample
  if(length(unique_filtered_rownames) == 0){
    return(list(distance='(0/1)',allele_names=c('neg'))) ## Changed from 'DID NOT PASS'
  }
  
  ## Ordering the rownames
  rowname_order <- order(as.integer(tstrsplit(unique_filtered_rownames, 'X')[[2]]))
  unique_filtered_rownames <- unique_filtered_rownames[rowname_order]
  
  ## Setting up a data.frame for holding all called snps across all VCF iterations
  called_snp_frame <- data.frame(matrix(nrow=length(unique_filtered_rownames),ncol=5), row.names = unique_filtered_rownames)
  not_calls <- list()
  
  ## Filling out called_snp_frame, as well as recording any bad calls
  for(snp_pos in unique_filtered_rownames){
    not_calls[[snp_pos]] <- c()
    for(calling_frame in allele_calling_frame_list_filtered){
      if(snp_pos %in% rownames(calling_frame)){
        geno_call <- calling_frame[snp_pos,10]
        geno_call <- tstrsplit(geno_call, ':')[[1]]
        snp_calls <- calling_frame[snp_pos,4:5]
        
        combined_snp_calls <- paste0(snp_calls, collapse=',')
        snp_call_vector <- strsplit(combined_snp_calls,',')[[1]]
        
        geno_call_vector <- strsplit(geno_call, '/')[[1]]
        geno_call_vector <- as.integer(unique(geno_call_vector))
        geno_call_vector <- geno_call_vector + 1
        
        if(!(1 %in% geno_call_vector)){
          not_calls[[snp_pos]] <- unique(c(not_calls[[snp_pos]], snp_call_vector[1]))
        }
        
        snp_geno_calls <- snp_call_vector[geno_call_vector]
        
        for(snp_call in snp_geno_calls){
          if(snp_call %in% unlist(called_snp_frame[snp_pos,])){
            next
          }else{
            col_pos_in_frame <- sum(!is.na(called_snp_frame[snp_pos,]))+1
            called_snp_frame[snp_pos,col_pos_in_frame] <- snp_call
          }
        }
      }
    }
  }
  
  ## Iterating through the bad calls to cut down the possible allele frame
  for(snp_pos in names(not_calls)){
    not_snp_vector <- not_calls[[snp_pos]]
    
    ## If a bad call is in a non-variant position, assume it is from a contaminating read and move on to the next bad call
    ## This could be a bad assumption, it will be something to look back on
    if(!(snp_pos %in% colnames(possible_allele_frame))){
      cat('\nContaminating bad call.\n')
      next
    }
    
    
    for(not_snp in not_snp_vector){
      
      if(not_snp %in% called_snp_frame[snp_pos,]){
        cat('\n\nSkipping bad call that was also a good call.\n\n')
        next
      }
      
      good_alleles <- rownames(possible_allele_frame)[possible_allele_frame[,snp_pos,drop=F] != not_snp]
      possible_allele_frame <- possible_allele_frame[good_alleles,,drop=F]
      
      if(nrow(possible_allele_frame) == 0){
        return(list(distance='(0/0)',allele_names=c('newR')))
      }
    }
  }
  
  het_positions <- rownames(called_snp_frame)[is_nuc(called_snp_frame[,'X2'])]
  
  
  if(!all(het_positions %in% variant_positions)){
    ## This will probably be changed to catch new alleles
    cat('\nNot all het calls are at variant positions. This likely indicates contaminating reads. Check out.\n')
    cat('For now these het calls are being excluded.\n')
    cat('POS: ', paste0(het_positions, collapse=', '))
    het_positions <- intersect(het_positions, variant_positions)
  }
  
  called_snp_frame <- called_snp_frame[colnames(original_possible_allele_frame),]
  
  ## If there are no het positions, then only call a single allele, no matter what n_alleles is
  if(all(is.na(called_snp_frame[,2]))){
    n_alleles <- 1
  }
  
  possible_allele_frame <- make_unique_pos_frame(possible_allele_frame)
  
  if(ncol(possible_allele_frame) == 0){
    possible_allele_frame <- original_possible_allele_frame[rownames(possible_allele_frame),]
  }
  
  ## Het positions, this is in response to a bad_call that dropped the correct allele
  ## the solutions I am going to use is to make sure the het positions for the correct allele
  ## are not dropped so the called alleles get the correct penalty

  
  if(length(het_positions) > 0){
    
    ## If the possible allele frame is missing these positions, add them back in
    if(!all(het_positions %in% colnames(possible_allele_frame))){
      cat('\nSome het positions were dropped for this sample due to bad_call.\n')
      to_be_added_to_score <- length(setdiff(het_positions, colnames(possible_allele_frame)))
    }else{
      to_be_added_to_score <- 0
    }
  }else{
    to_be_added_to_score <- 0
  }
  
  important_snp_pos <- intersect(rownames(called_snp_frame), colnames(possible_allele_frame))
  
  ## excluded SNP positions now count toward allele mismatch score
  base_allele_mismatch_score <- ncol(possible_allele_frame) - length(important_snp_pos)
  allele_mismatch_scale <- ncol(possible_allele_frame)
  
  called_snp_frame <- called_snp_frame[important_snp_pos,]
  possible_allele_frame <- possible_allele_frame[,important_snp_pos,drop=F]
  
  if(nrow(possible_allele_frame) > 1){
    possible_allele_frame <- make_unique_pos_frame(possible_allele_frame)
  }else if(nrow(possible_allele_frame) == 1){
    n_alleles = 1
  }else{
    cat('\nNo rows in possible allele frame!!')
    stop()
  }
  
  ## Single allele calling
  if(n_alleles == 1){
    single_allele_distance_frame <- data.frame(matrix(nrow=nrow(possible_allele_frame),ncol=1), row.names = rownames(possible_allele_frame))
    colnames(single_allele_distance_frame) <- 'distance'
    
    for(allele_name in rownames(single_allele_distance_frame)){
      allele_score <- 0
      allele_sequence <- possible_allele_frame[allele_name,,drop=F]
      
      for(snp_pos in colnames(allele_sequence)){
        allele_snp <- allele_sequence[,snp_pos]
        
        called_snp_vector <- unlist(called_snp_frame[snp_pos,])
        called_snp_vector <- called_snp_vector[!is.na(called_snp_vector)]
        
        if(!(allele_snp %in% called_snp_vector)){
          allele_score <- allele_score+1
        }
        
        extra_snps <- sum(allele_snp != called_snp_vector)
        
        allele_score <- allele_score+extra_snps
      }
      single_allele_distance_frame[allele_name,1] <- allele_score + to_be_added_to_score
    }
    min_allele_distance <- min(single_allele_distance_frame)
    min_allele_names <- rownames(single_allele_distance_frame[single_allele_distance_frame == min_allele_distance,,drop=F])
    cat('\n')
    cat(min_allele_distance)
    cat('\n')
    cat(min_allele_names)
  }
  
  ## Multi allele calling
  if(n_alleles > 1){
    ## Hardcoding the number of returned alleles to be 2
    n_alleles <- 2
    
    multi_rownames <- combinations(n=nrow(possible_allele_frame), r=n_alleles, v=rownames(possible_allele_frame))
    multi_rownames <- apply(multi_rownames,1,paste0,collapse='+')
    multi_allele_distance_frame <- data.frame(matrix(nrow=length(multi_rownames),ncol=1),row.names=multi_rownames)
    colnames(multi_allele_distance_frame) <- 'distance'
    
    for(multi_allele_name in rownames(multi_allele_distance_frame)){
      allele_score <- 0
      allele_vector <- strsplit(multi_allele_name, '+', fixed=T)[[1]]
      
      allele_sequence <- possible_allele_frame[allele_vector,,drop=F]
      for(snp_pos in colnames(allele_sequence)){
        allele_snps <- allele_sequence[,snp_pos]
        allele_snps <- unique(allele_snps)
        called_snp_vector <- unlist(called_snp_frame[snp_pos,])
        called_snp_vector <- called_snp_vector[!is.na(called_snp_vector)]
        
        missed_snps <- length(setdiff(called_snp_vector, allele_snps))
        extra_snps <- length(setdiff(allele_snps,called_snp_vector))
        allele_score <- allele_score + missed_snps + extra_snps
      }
      multi_allele_distance_frame[multi_allele_name,1] <- allele_score + to_be_added_to_score
    }
    
    min_allele_distance <- min(multi_allele_distance_frame)
    min_allele_names <- rownames(multi_allele_distance_frame[multi_allele_distance_frame == min_allele_distance,,drop=F])
    cat('\n')
    cat(min_allele_distance)
    cat('\n')
    cat(min_allele_names)
  }
  #min_allele_distance <- round((min_allele_distance / length(min_allele_names)), 4)
  base_allele_mismatch_scale <- round((base_allele_mismatch_score / allele_mismatch_scale), 4)
  min_allele_distance <- paste0('(', min_allele_distance, '/', base_allele_mismatch_scale, ')')
  return(list(distance=min_allele_distance,allele_names=min_allele_names))
}

## Make possible_allele_frame and locus specific vcf result frames to use with the allele_caller function
allele_caller_input_formatting <- function(current_locus, all_locus_result_frame_list, haplo_directory, allele_blacklist_table){
  current_locus_allele_frame <- msf_to_allele_frame(haplo_directory, current_locus)
  
  current_locus_allele_frame <- unique(current_locus_allele_frame) ## This eliminates alleles with the same coding sequence
  
  if(current_locus %in% colnames(allele_blacklist_table)){
    current_locus_allele_frame <- current_locus_allele_frame[!(rownames(current_locus_allele_frame) %in% allele_blacklist_table[,current_locus]),] ## Remove alleles on blacklist
  }
  
  unique_pos_frame <- make_unique_pos_frame(current_locus_allele_frame)
  colnames(unique_pos_frame) <- paste0(colnames(unique_pos_frame), 'X')
  
  allele_calling_frame_list <- all_locus_result_frame_list[[current_locus]]
  
  possible_allele_frame <- unique_pos_frame
  
  # Returning the two important inputs for the allele_caller function
  return(list(possible_allele_frame=possible_allele_frame,allele_calling_frame_list=allele_calling_frame_list))
}

## Takes in a '(x,x)' error score and returns the missing error (the second x).
distance_missing <- function(distance_value){
  missing_value <- strsplit(distance_value, '/')[[1]][2]
  missing_value <- strsplit(missing_value, ')')[[1]][1]
  missing_value <- as.double(missing_value)
  return(missing_value)
}

## Takes in a row of the distance_frame, checks to see if the 2DL2 or 2DL3 missing error value is above the threshold, indicating that the locus is not present
## and GC should be adjusted
gc_corrector_2DL23 <- function(distance_frame_row, missing_threshold = 0.5){
  return_list <- list()
  vector_2DL23 <- distance_frame_row[,c('2DL2','2DL3'),drop=F]
  
  missing_data <- any(is.na(vector_2DL23))
  if(missing_data){
    return_list[['2DL2']] <- 1
    return_list[['2DL3']] <- 1
  }else{
    missing_2DL2 <- distance_missing(vector_2DL23$`2DL2`) >= missing_threshold
    missing_2DL3 <- distance_missing(vector_2DL23$`2DL3`) >= missing_threshold
    
    if(missing_2DL2 & missing_2DL3){
      cat('\n\nSomething went wrong with 2DL23 GC. Skipping this sample.')
      return('DID NOT PASS')
    }else if(missing_2DL2){
      return_list[['2DL2']] <- 0
      return_list[['2DL3']] <- 2
    }else if(missing_2DL3){
      return_list[['2DL2']] <- 2
      return_list[['2DL3']] <- 0
    }else{
      return_list[['2DL2']] <- 1
      return_list[['2DL3']] <- 1
    }
  }
  
  return(return_list)
}

remove_sample_result_directory <- function(results.directory, sample_id){
  file_path_list <- list.files(results.directory, recursive = T, pattern=sample_id, full.names = T)
  file.remove(file_path_list)
  
  dir_path_list <- list.files(results.directory, recursive = T, pattern=sample_id, full.names = T, include.dirs = T)
  
  for(dir_path in dir_path_list){
    file.remove(list.files(dir_path, recursive=T, include.dirs=T, full.names=T))
    file.remove(dir_path)
  }
}

## Speed-up stuff
haplo_resources_directory <- 'Resources/haplo_resources'
kir_gen_path <- file.path(haplo_resources_directory, 'KIR_gen.fasta')
kir_nuc_path <- file.path(haplo_resources_directory, 'KIR_nuc.fasta')


