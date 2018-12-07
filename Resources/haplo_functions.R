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
                     '3DS1'=c(34, 36, 285, 300, 294, 51, 106, 51, 8),
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
                       '3DS1'=c(34, 36, 285, 300, 294, 51, 106, 51, 8),
                       '2DP1'=c(34, 36, 299, 294, 51, 105, 53, 177))
  first=T
  for(allele in allele_list){
    current_locus <- unlist(strsplit(allele, '*', fixed=T))[1]
    current_locus <- unlist(strsplit(current_locus, 'KIR', fixed=T))[2]
    if(!(current_locus %in% names(exon_lengths))){
      stop(paste0(allele, ': locus not found in exon_lengths'))
    }
    gen_chr_string <- ipdkir_gen_df[allele,]
    current_locus_allele_frame <- msf_to_allele_frame(msf_directory, current_locus)
    
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
          stop("\n>>> Stopped because of a duplication\n")
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
msf_to_allele_frame <- function(msf_directory, current_locus){
  if(current_locus == '2DL5A' | current_locus == '2DL5B'){
    old_locus <- current_locus
    current_locus <- '2DL5'
  }
  ## Read in MSF file for current_locus in msf_directory
  msf_file_path <- file.path(msf_directory,paste0('KIR',current_locus,'_nuc.msf'))
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
    #cut_off <- vcf_threshold*mean(dp_list)
    vcf_frame <- vcf_frame[(dp_list >= vcf_threshold),,drop=F]
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

## Call new alleles from integrated VCF output from all rounds of haplotye alignment
new_allele_caller <- function(SOS_locus_lookup,called_snp_frame,current_locus,DPthresh=6, mismatchThresh = 1){
  # Initialize Vector of reference KIR alleles for formatting new allele names
  #  - Alleles used as reference for IPD-KIR multiple alignments for each KIR locus as of June 2018
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
  
  # Store read position info
  all_reads_pos <- row.names(called_snp_frame)
  reads_pos <- row.names(called_snp_frame)
  reads_pos <- reads_pos[reads_pos %in% names(SOS_locus_lookup)]
  
  
  ##enumerate all possible genotypes
  hmat <- as.matrix(rbind(called_snp_frame$X1, called_snp_frame$X2))
  # Fill out second row beased on the current combined VCF organization, 
  # any position that is NA is the same as the the ref call in row 1
  for (i in 1:length(hmat[1,])){
    ref1 <- hmat[1,i]
    ref2 <- hmat[2,i]
    if(is.na(ref2)){
      hmat[2,i] <- ref1
    }
  }
  
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
  lookitup$string <- gsub('\\*','-',lookitup$string)   # Substitute '*' for '-'
  ##find allelic ambiguities in current lookup table
  
  all_amb<-c(1:length(lookitup))
  for(i in seq_along(lookitup$string)) {
    query_string <- gsub('\\*','-',lookitup$string[i]) # Substitute '*' for '-'
    l <- str_detect(lookitup$string,query_string)
    ll <- lookitup$allele[l]
    
    all_amb[i] <- paste(ll, collapse="/")
  }
  
  lookitup <- cbind(lookitup,as.character(all_amb), stringsAsFactors=FALSE)
  names(lookitup) <- c("string", "allele", "amb")
  
  ###allle calling
  allele1 <- ifelse(poss_genos$X1 %in% lookitup$string, (lookitup[match(poss_genos$X1, lookitup$string), 3]),"new")
  allele2 <- ifelse(poss_genos$X2 %in% lookitup$string, (lookitup[match(poss_genos$X2, lookitup$string), 3]),"new")
  
  genos <- data.frame(cbind(as.character(allele1), as.character(allele2)),stringsAsFactors = F)
  genos$string1 <- poss_genos$X1
  genos$string2 <- poss_genos$X2
  
  #### CREATE NEW ALLELE DATA FRAME HERE!!!  ################################################
  pos_vec <- reads_pos  # potentially remove 'X's from position names
  # Conversion from Genomic to CDS positions
  #pos_vec_conv        <- reads_pos_cds
  #names(pos_vec_conv) <- reads_pos_genomic
  
  new_alleles_vec      <- c()
  new_alleles_seq_vec  <- c()
  new_alleles_geno_vec <- c()
  new_alleles_categ    <- c() # possible genotypes will always be either new+known(allele) or new+new
  # Create Data Frame of all variable positions for new allele variants
  new_alleles_pos_var_df <- data.frame(matrix(nrow=0, ncol=length(pos_vec)))
  names(new_alleles_pos_var_df) <- pos_vec
  
  for(i in 1:length(row.names(genos))){
    A1                <- genos[i,1]
    A2                <- genos[i,2]
    A1_seq            <- genos[i,3]
    A2_seq            <- genos[i,4]
    A1_seq_vec        <- unlist(strsplit(A1_seq,""))
    A2_seq_vec        <- unlist(strsplit(A2_seq,""))
    names(A1_seq_vec) <- all_reads_pos
    names(A2_seq_vec) <- all_reads_pos
    
    Aref <- ref_allele_names[grep(current_locus,ref_allele_names)]
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
        new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec)))
        names(new_allele_entry)              <- pos_vec
        new_allele_entry[1,]                 <- "N"
        new_allele_entry[,names(A1_seq_vec)] <- A1_seq_vec
        row.names(new_allele_entry)          <- new_allele_name
        names(new_allele_entry)              <- as.character(pos_vec)
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
        new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec)))
        names(new_allele_entry)              <- pos_vec
        new_allele_entry[1,]                 <- "N"
        new_allele_entry[,names(A1_seq_vec)] <- A1_seq_vec
        row.names(new_allele_entry)          <- a1_new_allele_name
        names(new_allele_entry)              <- as.character(pos_vec)
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
        a1_new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec)))
        names(a1_new_allele_entry)              <- pos_vec
        a1_new_allele_entry[1,]                 <- "N"
        a1_new_allele_entry[,names(A1_seq_vec)] <- A1_seq_vec
        row.names(a1_new_allele_entry)          <- a1_new_allele_name
        names(a1_new_allele_entry)              <- as.character(pos_vec)
        new_alleles_pos_var_df                  <- rbind(new_alleles_pos_var_df, a1_new_allele_entry)
        
        a2_new_allele_entry                     <- data.frame(matrix(nrow=1, ncol=length(pos_vec)))
        names(a2_new_allele_entry)              <- pos_vec
        a2_new_allele_entry[1,]                 <- "N"
        a2_new_allele_entry[,names(A2_seq_vec)] <- A2_seq_vec
        row.names(a2_new_allele_entry)          <- a2_new_allele_name
        names(a2_new_allele_entry)              <- as.character(pos_vec)
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
  new_alleles_df$sample <- sample_id
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

## Call alleles from locus specific vcf result frames
allele_caller <- function(allele_calling_frame_list, possible_allele_frame, n_alleles, DPthresh = 6){
  
  ## Saving a copy of the possible allele frame, since it will be cut down during this function
  original_possible_allele_frame <- possible_allele_frame
  variant_positions <- colnames(possible_allele_frame)
  
  ## Filtering any calls in the VCF files that do not pass the threshold
  allele_calling_frame_list_filtered <- lapply(allele_calling_frame_list, filter_vcf_frame, DPthresh)
  
  ## Finding all positions that are in any VCF frame after filtration
  filtered_rownames <- lapply(allele_calling_frame_list_filtered, rownames)
  unique_filtered_rownames <- unique(unlist(filtered_rownames))
  
  ## If there are no filtered rownames, then say this locus is negative for this sample
  if(length(unique_filtered_rownames) == 0){
    return(list(distance='(0/1)',allele_names=c('nocall_no_snps_left_after_filter'), called_snp_frame=NULL)) ## Changed from 'DID NOT PASS'
  }
  
  ## Ordering the rownames
  rowname_order <- order(as.integer(tstrsplit(unique_filtered_rownames, 'X')[[2]]))
  unique_filtered_rownames <- unique_filtered_rownames[rowname_order]
  
  ## Setting up a data.frame for holding all called snps across all VCF iterations
  called_snp_frame <- data.frame(matrix(nrow=length(unique_filtered_rownames),ncol=5), row.names = unique_filtered_rownames)
  not_calls <- list()
  
  ## Initialize a list for counting how many rounds a snp shows up for each position
  snp_pos_count <- list()

  ## Filling out called_snp_frame, as well as recording any bad calls
  for(snp_pos in unique_filtered_rownames){
    
    ## initalize a vector of bad calls for this position
    not_calls[[snp_pos]] <- c()
    
    ## initialize the count for this snp position
    snp_pos_count[[snp_pos]] <- 0

    for(calling_frame in allele_calling_frame_list_filtered){
      if(snp_pos %in% rownames(calling_frame)){
        
        ## Iterate the count for this position
        snp_pos_count[[snp_pos]] <- as.integer(snp_pos_count[[snp_pos]]) + 1
        
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
            col_pos_in_frame <- which(snp_call == unlist(called_snp_frame[snp_pos,]))[[1]]
            #called_snp_frame[snp_pos,col_pos_in_frame+1] <- called_snp_frame[snp_pos,col_pos_in_frame+1]+1
            next
          }else{
            col_pos_in_frame <- sum(!is.na(called_snp_frame[snp_pos,]))+1
            called_snp_frame[snp_pos,col_pos_in_frame] <- snp_call
            #called_snp_frame[snp_pos,col_pos_in_frame+1] <- 1
          }
        }
      }
    }
  }
  
  ## This section removes SNPs that are not found across a minimum of min_rounds or snp_round_count, whichever is lowest
  #snp_round_count_cols <- c('X2','X4','X6','X8','X10')
  #for(snp_pos in names(snp_pos_count)){
  #  current_pos_count <-snp_pos_count[[snp_pos]]
  #  current_cols <- snp_round_count_cols[!is.na(called_snp_frame[snp_pos,snp_round_count_cols])]
  #  
  #  current_cols <- current_cols[!(called_snp_frame[snp_pos,current_cols] >= min_rounds)]
  #  
  #  if(length(current_cols) == 0){
  #    next
  #  }else{
  #    bad_cols <- current_cols[called_snp_frame[snp_pos,current_cols] != current_pos_count]
  #    
  #    for(bad_col in bad_cols){
  #      bad_col_index <- which(bad_col == colnames(called_snp_frame))
  #      called_snp_frame[snp_pos,bad_col_index] <- NA
  #      called_snp_frame[snp_pos,bad_col_index-1] <- NA
  #    }
  #  }
  #}
  
  #called_snp_frame <- called_snp_frame[,c('X1','X3','X5','X7','X9')]
  #colnames(called_snp_frame) <- c('X1','X2','X3','X4','X5')
  
  
  ## Iterating through the bad calls to cut down the possible allele frame
  for(snp_pos in names(not_calls)){
    not_snp_vector <- not_calls[[snp_pos]]
    
    ## If a bad call is in a non-variant position, assume it is from a contaminating read and move on to the next bad call
    ## This could be a bad assumption, it will be something to look back on
    if(!(snp_pos %in% colnames(possible_allele_frame))){
      cat('\nContaminating bad call.\n')
      return(list(distance='(0/1)',allele_names=c('nocall_notSNP_found_at_invariant_position'),called_snp_frame=called_snp_frame))
    }
    
    
    for(not_snp in not_snp_vector){
      
      if(not_snp %in% called_snp_frame[snp_pos,]){
        cat('\n\nSkipping bad call that was also a good call.\n\n')
        next
      }
      
      good_alleles <- rownames(possible_allele_frame)[possible_allele_frame[,snp_pos,drop=F] != not_snp]
      possible_allele_frame <- possible_allele_frame[good_alleles,,drop=F]
      
      if(nrow(possible_allele_frame) == 0){
        return(list(distance='(0/1)',allele_names=c('nocall_all_possible_alleles_dropped_due_to_notSNPS'),called_snp_frame=called_snp_frame))
      }
    }
  }
  
  het_positions <- rownames(called_snp_frame)[is_nuc(called_snp_frame[,'X2'])]
  
  
  if(!all(het_positions %in% variant_positions)){
    ## This will probably be changed to catch new alleles
    cat('\nNot all het calls are at variant positions. This likely indicates contaminating reads or new allele.\n')
    cat('Returning a no call.\n')
    cat('POS: ', paste0(het_positions, collapse=', '))
    #het_positions <- intersect(het_positions, variant_positions)
    return(list(distance='(0/1)',allele_names=c('nocall_het_call_at_invariant_pos'),called_snp_frame=called_snp_frame))
  }
  
  called_snp_frame <- called_snp_frame[colnames(original_possible_allele_frame),]
  
  ## If there are no het positions, then only call a single allele, no matter what n_alleles is
  if(all(is.na(called_snp_frame[,2]))){
    n_alleles <- 1
  }else{
    n_alleles <- 2
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
    return(list(distance='(0/1)',allele_names=c('nocall_no_allele_differentiating_SNPS_left'),called_snp_frame=called_snp_frame))
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
    multi_allele_distance_frame <- data.frame(matrix(nrow=length(multi_rownames),ncol=3),row.names=multi_rownames)
    colnames(multi_allele_distance_frame) <- c('distance','extra','missing')
    
    for(multi_allele_name in rownames(multi_allele_distance_frame)){
      allele_score <- 0
      total_missed_snps <- 0
      total_extra_snps <- 0
      allele_vector <- strsplit(multi_allele_name, '+', fixed=T)[[1]]
      
      allele_sequence <- possible_allele_frame[allele_vector,,drop=F]
      for(snp_pos in colnames(allele_sequence)){
        allele_snps <- allele_sequence[,snp_pos]
        allele_snps <- unique(allele_snps)
        
        ## Excluding deletion reference symbols ('.') from allele calling
        allele_snps <- allele_snps[unlist(lapply(allele_snps, is_nuc))]
        called_snp_vector <- unlist(called_snp_frame[snp_pos,])
        called_snp_vector <- called_snp_vector[!is.na(called_snp_vector)]
        
        missed_snps <- length(setdiff(called_snp_vector, allele_snps))
        extra_snps <- length(setdiff(allele_snps,called_snp_vector))
        allele_score <- allele_score + missed_snps + extra_snps
        total_missed_snps <- total_missed_snps + missed_snps
        total_extra_snps <- total_extra_snps + extra_snps
      }
      multi_allele_distance_frame[multi_allele_name,1] <- allele_score + to_be_added_to_score
      multi_allele_distance_frame[multi_allele_name,2] <- total_extra_snps
      multi_allele_distance_frame[multi_allele_name,3] <- total_missed_snps
    }
    
    #min_missing_distance <- min(multi_allele_distance_frame$missing)
    #min_allele_names <- rownames(multi_allele_distance_frame[multi_allele_distance_frame$missing == min_missing_distance,,drop=F])
    
    #min_extra_distance <- min(multi_allele_distance_frame[min_allele_names,'extra'])
    #min_allele_names <- min_allele_names[multi_allele_distance_frame[min_allele_names,'extra'] == min_extra_distance]
    
    ## Subtracting the pentaly for having missing snps to makeup for taking out snps that do not pass the min_rounds threshold
    ## to clarify the naming, 'extra' means extra snps in the called alleles that were not found in the called_snp_frame
    #min_allele_distance <- min(multi_allele_distance_frame[min_allele_names,'distance'] - multi_allele_distance_frame[min_allele_names,'extra'])
    min_allele_distance <- min(multi_allele_distance_frame$distance)
    min_allele_names <- rownames(multi_allele_distance_frame[multi_allele_distance_frame$distance == min_allele_distance,,drop=F])
    
    cat('\n')
    cat(min_allele_distance)
    cat('\n')
    cat(min_allele_names)
  }
  #min_allele_distance <- round((min_allele_distance / length(min_allele_names)), 4)
  
  base_allele_mismatch_scale <- round((base_allele_mismatch_score / allele_mismatch_scale), 4)
  return_min_allele_distance <- paste0('(', min_allele_distance, '/', base_allele_mismatch_scale, ')')
  
  if(min_allele_distance > 0){
    return(list(distance=return_min_allele_distance,allele_names=c('new_does_not_perfectly_match_known_alleles'),called_snp_frame=called_snp_frame))
  }
  
  return(list(distance=return_min_allele_distance,allele_names=min_allele_names,called_snp_frame=NULL))
}

## Make possible_allele_frame and locus specific vcf result frames to use with the allele_caller function
allele_caller_input_formatting <- function(current_locus, all_locus_result_frame_list, msf_directory, allele_blacklist_table){
  current_locus_allele_frame <- msf_to_allele_frame(msf_directory, current_locus)
  
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

## Reading in the gc_input.csv file
read_master_gc <- function(master_gc_path){
  gc_path <- normalizePath(master_gc_path)
  gc_table <- read.table(gc_path, sep=',', check.names=F, stringsAsFactors=F, header=T, row.names=1)
  return(gc_table)
}

## Reading in the possible haplotype reference alleles
read_master_haplo <- function(master_haplo_path){
  haplo_path <- normalizePath(master_haplo_path)
  haplo_table <- read.table(haplo_path, sep=',', check.names=F, stringsAsFactors=F, header=T, row.names=1)
  return(haplo_table)
}

## Read in the list of alleles not to call
read_blacklist <- function(allele_blacklist_path){
  blacklist_path <- normalizePath(allele_blacklist_path)
  blacklist_table <- read.table(blacklist_path, sep=',', check.names=F, stringsAsFactors=F, header=T)
  return(blacklist_table)
}

## Determine 2DL2_3 copy number and write gc_input.csv
prepare_gc_input <- function(raw.kff.file='',combined.csv.file='',results.directory=''){
  kff_results_table <- read.csv(raw.kff.file, stringsAsFactors = F, check.names = F)
  combined_table <- read.csv(combined.csv.file, stringsAsFactors = F, check.names = F, row.names = 1)
  #combined_table <- t(combined_table)
  combined_table_append <- data.frame(matrix('0',2,ncol(combined_table)),stringsAsFactors=F,check.names=F,row.names=c('2DL2','2DL3'))
  colnames(combined_table_append) <- colnames(combined_table)
  combined_table<-rbind(combined_table, combined_table_append)
  
  ## only keep sample names that are found in both tables
  kff_samples <- rownames(kff_results_table)
  com_samples <- colnames(combined_table)
  in_both <- intersect(kff_samples,com_samples)
  
  kff_results_table <- kff_results_table[in_both,]
  combined_table <- combined_table[,in_both,drop=F]
  
  bool_index_2DL2 <- unlist(lapply(colnames(kff_results_table), FUN = (function(x) grepl('*2DL2*short', x, fixed=T))))
  bool_index_2DL3 <- unlist(lapply(colnames(kff_results_table), FUN = (function(x) grepl('*2DL3*2', x, fixed=T))))
  
  threshold <- 10
  weird_samples <- c()
  for(sample_id in rownames(kff_results_table)){
    positive_2DL2 <- sum(kff_results_table[sample_id,bool_index_2DL2]) >= threshold
    positive_2DL3 <- sum(kff_results_table[sample_id,bool_index_2DL3]) >= threshold
    
    if(positive_2DL2 & positive_2DL3){
      combined_table['2DL2',sample_id] <- 1
      combined_table['2DL3',sample_id] <- 1
    }else if(positive_2DL2){
      combined_table['2DL2',sample_id] <- 2
      combined_table['2DL3',sample_id] <- 0
    }else if(positive_2DL3){
      combined_table['2DL2',sample_id] <- 0
      combined_table['2DL3',sample_id] <- 2
    }else{
      combined_table['2DL2',sample_id] <- 0
      combined_table['2DL3',sample_id] <- 0
    }
  }
  
  write.csv(t(combined_table), file = file.path(results.directory,'gc_input.csv'), quote = F)
}

ping_haplo_aligner_version <- '1.1'
cat(paste0('PING_haplo_aligner version: ', ping_haplo_aligner_version, '\n'))

ping_haplo_aligner <- function(
  sample.location='',
  fastq.pattern.1 = "_1.fastq.gz",
  fastq.pattern.2 = "_2.fastq.gz",
  bowtie.threads = 4,
  results.directory = '',
  sample.name = '',
  sample.haplotype = haplo_1,
  ipdkir_allele_df = ipdkir_allele_df,
  ipdkir_nuc_df = ipdkir_nuc_df,
  haplo.iteration=1
){
  library(data.table)
  library(stringr)
  source("Resources/haplo_functions.R", local = TRUE)
  
  haplo_resources_directory <- normalizePath(haplo_resources_directory, mustWork=T)
  
  # Creates results directory, defaults to Aligner_results
  create_results_directory <- function() {
    cat("----- Getting PING ready -----\n\n")
    cat(paste0("Current working directory: ", getwd(), '\n\n'))
    
    if(results.directory != ""){
      save_to <- results.directory
    }else{
      save_to <- paste0("Haplo_aligner_results_", format(Sys.time(), "%Y_%m_%d_%H_%M"), "/")
      
      count <- 1
      while(file.exists(save_to)) {
        save_to <- paste0("Haplo_aligner_results_", format(Sys.time(), "%Y_%m_%d_%H_%M"), "_", count, "/")
        count <- count + 1
      }
    }
    dir.create(file.path(save_to), showWarnings = F)
    save_to <- normalizePath(save_to)
    
    cat(paste("Results being saved to", save_to, "\n\n"))
    
    return(save_to)
  }
  
  # Creates results folders
  ping_ready <- function(current_sample, results_directory) {
    
    dir.create(file.path(results_directory), showWarnings = F)
    dir.create(file.path(results_directory, "Vcf"), showWarnings = F)
    dir.create(file.path(results_directory, "Fastq"), showWarnings = F)
    dir.create(file.path(results_directory, "Haplo_reference"), showWarnings=F)
    
    sample <- current_sample
    dir.create(file.path(results_directory, 'Vcf', sample), showWarnings=F)
    dir.create(file.path(results_directory, 'Vcf', sample, 'locus'), showWarnings=F)
    dir.create(file.path(results_directory, 'Vcf', sample, 'haplo'), showWarnings=F)
    dir.create(file.path(results_directory, 'Fastq', sample), showWarnings=F)
    dir.create(file.path(results_directory, 'Fastq', sample, 'locus'), showWarnings=F)
    dir.create(file.path(results_directory, 'Fastq', sample, 'haplo'), showWarnings=F)
    dir.create(file.path(results_directory, 'Haplo_reference', sample), showWarnings=F)
    
    cat("\n\nResults subdirectories created.\n\n")
  }
  
  # Create next reference iteration directory
  create_new_iteration_directory <- function(sample, results_directory, haplo_iteration = 1) {
    ## Finding how many haplotype reference iterations have been performed for a specific sample
    ## and creating a new directory for the next iteration
    
    haplo_ref_path <- file.path(results_directory, 'Haplo_reference', sample)
    haplo_dirs <- list.dirs(haplo_ref_path, recursive=F, full.names=F)
    
    #if(length(haplo_dirs) == 0){
    #  haplo_iteration <- 1
    #}else{
    #  haplo_iteration_raw <- tstrsplit(haplo_dirs, '_')[[2]]
    #  haplo_iteration <- max(as.integer(haplo_iteration_raw)) + 1
    #}
    
    save_to <- file.path(results_directory, 'Haplo_reference', sample, paste0('reference_', haplo_iteration))
    if (dir.exists(save_to)){
      unlink(save_to,recursive = T)
    }
    dir.create(save_to, showWarnings=F)
    
    iteration_filepath <- normalizePath(save_to)
    return(list(filepath = iteration_filepath, iteration = haplo_iteration))
  }
  
  # Finds sequences
  get_sample <- function() {
    
    # This is the non-recursive version
    sample_list = list.files(file.path(sample.location), pattern = sample.name)
    
    if (is.na(sample_list[1])) {
      string <- paste("No sequences found in", sample.location, "using ", sample.name)
      stop(string)
    } else {
      
      if(!any(grepl(fastq.pattern.1, sample_list, fixed=T))){
        string <- paste0("Fastq pattern ", fastq.pattern.1, " not found in ", sample_list[1])
        stop(string)
      }
      if(!any(grepl(fastq.pattern.1, sample_list, fixed=T))){
        string <- paste0("Fastq pattern ", fastq.pattern.2, " not found in ", sample_list[2])
        stop(string)
      }
      
      sample_list <- gsub(fastq.pattern.1, "", sample_list)
      sample_list <- gsub(fastq.pattern.2, "", sample_list)
      sample_list <- unique(sample_list)
      cat(paste("Found sequences: ", paste(sample_list, collapse = "\n"), sep = "\n"))
      cat("\n")
      return(sample_list)
    }
  }
  
  ### Aligns Sample to Haplotype reference
  align_to_haplo_ref <- function(fasta_path, bed_file, sequence, iteration){
    ### 1. Align KIR extracted reads to haplo-reference
    bt2_p <- paste0("-p", bowtie.threads)
    bt2_5 <- "-5 3"
    bt2_3 <- "-3 7"
    #bt2_L <- "-L 20"
    bt2_i <- "-i S,1,0.5"
    bt2_min_score <- "--score-min L,0,-0.187"
    bt2_I <- "-I 75"
    bt2_X <- "-X 1000"
    
    bt_index_path <- substr(fasta_path, 1, (nchar(fasta_path)-6))
    bt2_x <- paste0("-x ", bt_index_path)
    
    bt2_1 <- paste0("-1 ", file.path(sample.location, paste0(sequence, fastq.pattern.1)))
    bt2_2 <- paste0("-2 ", file.path(sample.location, paste0(sequence, fastq.pattern.2)))
    
    sequence  <- last(unlist(strsplit(sequence, "/")))
    sam_file   <- paste0(sequence, ".haplo.sam")  # SAM output file path
    bt2_stream <- paste0("-S ", sam_file)    
    
    
    bt2_al_conc <- paste0("--al-conc-gz ", sequence, "_%.fastq.gz")
    
    bt2_un <- "--un dump.me"
    
    cat('\n\n', sequence,"\n\n")
    cat("bowtie2", bt2_p, bt2_5, bt2_3, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_stream, bt2_al_conc, bt2_un)
    system2("bowtie2", c(bt2_p, bt2_5, bt2_3, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_stream, bt2_al_conc, bt2_un))
    cat("\n\n")
    
    #file.remove(paste0(sequence, ".temp"))
    file.remove(paste0("dump.me"))
    
    
    ### 2. Convert SAM to BAM file
    samt_param <- "-b -q10 "
    bam_file <- paste0(sequence, ".haplo.bam")  # BAM output file
    samt_out <- paste0("-o ", bam_file)
    cat("samtools view", samt_param, sam_file, samt_out)
    system2("samtools", c("view",samt_param, sam_file, samt_out))
    cat("\n\n")
    
    ### 3. Sort BAM file
    samt_T <- "-T temp "
    sorted_bam_file <- paste0(sequence, ".haplo.sorted.bam")  # Sorted BAM output file
    samt_out_sorted <- paste0("-o ", sorted_bam_file)
    cat("samtools sort", samt_out_sorted, samt_T, bam_file)
    system2("samtools", c("sort", samt_out_sorted, samt_T, bam_file))
    cat("\n\n")
    
    ### 4. Generate VCF file
    mp_m <- "-m 3"
    mp_F <- "-F 0.0002"
    mp_u <- "-u "
    mp_f <- paste0("-f ", fasta_path)
    mp_l <- paste0("-l ", bed_file)
    vcf_filepath <- file.path(results_directory, 'Vcf', sequence, 'haplo', paste0(sequence, '.reference_', iteration, '.vcf'))
    mp_o <- paste0("-o ", vcf_filepath)
    mp_O <- "-O v"
    
    cat("samtools", "mpileup", mp_m, mp_u, mp_F, mp_f, sorted_bam_file, mp_l, "|", "bcftools", "call", "--multiallelic-caller", mp_O, mp_o)
    system2("samtools", c("mpileup", mp_m, mp_u, mp_F, mp_f, sorted_bam_file, mp_l, "|", "bcftools", "call", "--multiallelic-caller", mp_O, mp_o))
    cat("\n\n")
    
    
    file.remove(paste0(sequence,'_1.fastq.gz'))
    file.remove(paste0(sequence,'_2.fastq.gz'))
    file.remove(bam_file)
    file.remove(sam_file)
    file.remove(sorted_bam_file)
    
    
    return(vcf_filepath)
  }
  
  results_directory <- create_results_directory()
  current_sample <- get_sample()
  ping_ready(current_sample, results_directory)
  
  # Generate ipd_kir allele dataframe
  kir_gen_path <- file.path(ipd_kir_resources_directory, 'KIR_gen.fasta')
  kir_nuc_path <- file.path(ipd_kir_resources_directory, 'KIR_nuc.fasta')
  
  ## SPEED UP
  #ipdkir_allele_df <- get_ipdkir_allele_df(kir_gen_path)
  #ipdkir_nuc_df <- get_ipdkir_allele_df(kir_nuc_path)
  
  ## Generate iteration directory
  sample_haplo_iteration_list <- create_new_iteration_directory(current_sample, results_directory, haplo_iteration = haplo.iteration)
  haplo_ref_filepath <- file.path(sample_haplo_iteration_list$filepath, paste0(current_sample, '.haplo.fasta'))
  
  ## THIS FUNCTION INCORPORATES SINK
  generate_haplo_ref(sample.haplotype,ipdkir_allele_df, haplo_ref_filepath)
  
  ## Build Bowtie2 reference
  haplo_bowtie_index <- substr(haplo_ref_filepath, 1, nchar(haplo_ref_filepath)-6)
  system2("bowtie2-build", c(haplo_ref_filepath, haplo_bowtie_index))
  
  ## Make the bed file
  bed_filepath <- file.path(sample_haplo_iteration_list$filepath, paste0(current_sample, '.haplo.bed'))
  make_bed(haplo_resources_directory, sample.haplotype, bed_filepath, ipdkir_allele_df)
  
  ## Haplotype alignment
  vcf_filepath <- align_to_haplo_ref(haplo_ref_filepath, bed_filepath, current_sample, haplo.iteration)
  cat('\n\nAll finished!')
}

ping_haplo_caller_version <- '1.0'
cat(paste0('PING_haplo_caller version: ', ping_haplo_caller_version, '\n'))
ping_haplo_caller <- function(
  sample.name = '',
  results.directory = '',
  ipdkir_allele_df = ipdkir_allele_df,
  ipdkir_nuc_df = ipdkir_nuc_df
){
  library(data.table)
  library(stringr)
  #source("Resources/haplo_functions.R", local = TRUE)
  
  vcf_threshold = 0.1
  allele_threshold = 0.85
  
  results_directory <- normalizePath(results.directory, mustWork=T)
  
  # Finds sequence names
  get_full_sample_name <- function(sample_name, results_directory) {
    dir_list <- list.dirs(file.path(results_directory, 'Vcf'), full.names=F, recursive=F)
    if(length(dir_list) == 0){
      cat('\n\nNo sample directories were found in: ', results_directory)
      stop()
    }
    
    sample_index <- grep(sample_name, dir_list, fixed=T)
    
    if(length(sample_index) > 1){
      cat('\n\nMore than 1 directory matched: ', sample_name)
      stop()
    }else if(length(sample_index) == 0){
      cat('\n\nNo directories matched: ', sample_name)
      stop()
    }
    
    return(dir_list[sample_index])
  }
  
  call_alt_allele <- function(frame_allele_unique, vcf_frame){
    frame_allele_unique_updated <- frame_allele_unique
    vcf_frame <- filter_vcf_frame(vcf_frame, vcf_threshold)
    ref_allele <- vcf_frame[1,1]
    cat('\n\nRef allele: ', ref_allele)
    
    geno_bad_calls <- which(tstrsplit(vcf_frame$V10, ':')[[1]] == '1/1')
    
    if(length(geno_bad_calls)){
      cat('\n\nBad reference detected: ', ref_allele)
      #stop()
    }
    
    alt_pos_list <- which(vcf_frame$V5 != '.')
    alt_pos_list <- rownames(vcf_frame[alt_pos_list,])
    no_alt=F
    if(length(alt_pos_list) == 0){
      cat('\n\nNo alt calls found.')
      no_alt=T
    }else{
      cat('\n\n', length(alt_pos_list), ' alt calls found.')
    }
    
    ## Iterate through all alt calls
    for(current_pos in alt_pos_list){
      
      current_nuc <- vcf_frame[current_pos, 'V5']
      
      if(current_pos %in% colnames(frame_allele_unique_updated)){
        matched_alleles <- which(frame_allele_unique_updated[,current_pos] == current_nuc)
        frame_allele_unique_updated <- make_unique_pos_frame(frame_allele_unique_updated[matched_alleles,,drop=F])
      }else{
        ## Making sure cut out positions still match the nuc call
        sub_matched_alleles <- rownames(frame_allele_unique_updated)
        if(!all(frame_allele_unique[sub_matched_alleles, current_pos] == current_nuc)){
          frame_allele_unique_updated <- frame_allele_unique
          cat('\n\nCut position does not match SNP call.')
          break()
        }
      }
      
      if(nrow(frame_allele_unique_updated) == 0){
        cat('\n\nAll alleles dropped when trying to match.')
      }
    }
    
    ## First condition: No alt calls found
    
    if(no_alt){
      alleles_to_return <- list(primary=c(ref_allele), secondary=c())
    }else if(ncol(frame_allele_unique_updated) == 0){
      
      ## Second condition: Alt calls perfectly match an allele/s
      
      cat('\n\nPerfectly matched alt allele.')
      #alleles_to_return <- rownames(frame_allele_unique_updated)
      alleles_to_return <- list(primary=rownames(frame_allele_unique_updated), secondary=c())
    }else{
      
      ## Third condition: Ref and Alt calls will be matched to closest allele
      
      cat('\n\nFinding closest matching alt alleles.')
      snp_pos <- colnames(frame_allele_unique_updated)
      snp_nuc <- vcf_frame[snp_pos, 'V4']
      
      snp_frame <- data.frame(matrix(nrow=1, ncol=length(snp_pos)), row.names = c('SNP'))
      colnames(snp_frame) <- snp_pos
      snp_frame['SNP',] <- snp_nuc
      
      ## Bringing in ALT calls instead of only having REF calls
      alt_snp_list <- alt_pos_list[alt_pos_list %in% snp_pos]
      if(length(alt_snp_list) > 0){
        snp_frame['SNP',alt_snp_list] <- vcf_frame[alt_snp_list,'V5']
      }
      
      ## Drop na's
      snp_frame <- snp_frame[,!is.na(snp_frame), drop=F]
      frame_allele_unique_updated <- frame_allele_unique_updated[,colnames(snp_frame),drop=F]
      
      allele_distance_list <- find_allele_distances_from_snps(snp_frame, frame_allele_unique_updated)
      
      closest_alleles <- names(which(allele_distance_list$distance_list == allele_distance_list$min_distance))
      other_alleles <- names(which(allele_distance_list$distance_list > allele_distance_list$min_distance))
      
      min_distance <- allele_distance_list$min_distance
      
      cat('\n\nClosest allele distance: ', min_distance)
      cat('\nClosest alleles: ', closest_alleles)
      if(length(other_alleles) > 0){
        second_min_distance <- min(allele_distance_list$distance_list[other_alleles])
        second_closest_alleles <- names(which(allele_distance_list$distance_list == second_min_distance))
        cat('\n\nSecond closest allele distance: ', second_min_distance)
        cat('\nSecond closest alleles: ', second_closest_alleles)
        alleles_to_return <- list(primary=closest_alleles, secondary=second_closest_alleles)
      }else{
        cat('\n\nNo secondary alleles.')
        alleles_to_return <- list(primary=closest_alleles, secondary=c())
      }
    }
    
    return(alleles_to_return)
  }
  
  call_not_allele <- function(frame_allele_unique, vcf_frame){
    frame_allele_unique_updated <- frame_allele_unique
    ref_allele <- vcf_frame[1,1]
    cat('\n\nIncorrect ref allele: ', ref_allele)
    
    geno_bad_calls <- which(tstrsplit(vcf_frame$V10, ':')[[1]] == '1/1')
    bad_ref_pos_list <- rownames(vcf_frame[geno_bad_calls,])
    
    for(current_pos in bad_ref_pos_list){
      current_nuc <- vcf_frame[current_pos, 'V4']
      
      if(current_pos %in% colnames(frame_allele_unique_updated)){
        matched_alleles <- which(frame_allele_unique_updated[,current_pos] != current_nuc)
        frame_allele_unique_updated <- make_unique_pos_frame(frame_allele_unique_updated[matched_alleles,,drop=F])
      }else{
        ## Making sure cut out positions still match the nuc call
        sub_matched_alleles <- rownames(frame_allele_unique_updated)
        if(!all(frame_allele_unique[sub_matched_alleles, current_pos] != current_nuc)){
          stop('\n\nCut position does not match SNP call.')
        }else{
          cat('\n\nCut positions all match SNP call!')
        }
      }
    }
    
    geno_alt_calls <- which(tstrsplit(vcf_frame$V10, ':')[[1]] == '0/1')
    alt_call_list <- rownames(vcf_frame[geno_alt_calls,])
    
    for(current_pos in alt_call_list){
      ref_nuc <- vcf_frame[current_pos, 'V4']
      alt_nuc <- vcf_frame[current_pos, 'V5']
      
      current_nuc_list <- c(ref_nuc, alt_nuc, '.')
      
      if(current_pos %in% colnames(frame_allele_unique_updated)){
        matched_alleles <- which(frame_allele_unique_updated[,current_pos] %in% current_nuc_list)
        frame_allele_unique_updated <- make_unique_pos_frame(frame_allele_unique_updated[matched_alleles,,drop=F])
      }
    }
    
    
    if(ncol(frame_allele_unique_updated) == 0){
      cat('\n\nPerfectly matched alt reference allele.')
      alleles_to_return <- rownames(frame_allele_unique_updated)
    }else{
      cat('\n\nFinding closest matching alt reference alleles.')
      snp_pos <- colnames(frame_allele_unique_updated)
      snp_nuc <- vcf_frame[snp_pos, 'V4']
      
      snp_frame <- data.frame(matrix(nrow=1, ncol=length(snp_pos)), row.names = c('SNP'))
      colnames(snp_frame) <- snp_pos
      snp_frame['SNP',] <- snp_nuc
      
      ## Drop na's
      snp_frame <- snp_frame[,!is.na(snp_frame), drop=F]
      frame_allele_unique_updated <- frame_allele_unique_updated[,colnames(snp_frame),drop=F]
      
      allele_distance_list <- find_allele_distances_from_snps(snp_frame, frame_allele_unique_updated)
      
      closest_alleles <- names(which(allele_distance_list$distance_list == allele_distance_list$min_distance))
      min_distance <- allele_distance_list$min_distance
      
      cat('\n\nClosest allele distance: ', min_distance)
      alleles_to_return <- closest_alleles
    }
    cat('\n\nReturning possible alterative references.')
    return(alleles_to_return)
  }
  
  compare_Vcf_and_Haplo_reference <- function(vcf_directory, iteration_dir_list){
    both_are_same = TRUE
    
    num_of_vcf_files <- length(list.files(vcf_directory))
    
    if(length(iteration_dir_list) != num_of_vcf_files){
      both_are_same = FALSE
    }
    
    return(both_are_same)
  }
  
  ## Sample manipulation and setup
  sample_name <- get_full_sample_name(sample.name, results_directory)
  reference_directory <- normalizePath(file.path(results_directory, 'Haplo_reference', sample_name), mustWork=T)
  vcf_directory <- normalizePath(file.path(results_directory, 'Vcf', sample_name, 'haplo'), mustWork=T)
  iteration_dir_list <- list.dirs(file.path(results_directory, 'Haplo_reference', sample_name), recursive=F, full.names=F)
  
  ######
  # HARD RETURN IF NUMBER OF VCF FILES != NUMBER OF HAPLO REFERENCES
  sample_must_pass_this_checkpoint <- compare_Vcf_and_Haplo_reference(vcf_directory, iteration_dir_list)
  if(!sample_must_pass_this_checkpoint){
    return("DID NOT PASS")
  }
  #####
  # all good from here on out :)
  
  iteration_numbers <- as.integer(tstrsplit(iteration_dir_list, '_')[[2]])
  iteration_dir_list <- iteration_dir_list[order(iteration_numbers)]
  
  all_allele_calling_frame_list <- list()
  for(iteration_dir in iteration_dir_list){
    cat('\n\nReading in files from: ', iteration_dir)
    result_frame_list <- list()
    allele_frame_list <- list()
    allele_call_list <- list()
    allele_vcf_full_list <- list()
    
    bed_filepath <- file.path(reference_directory, iteration_dir, paste0(sample_name, '.haplo.bed'))
    bed_filepath <- normalizePath(bed_filepath, mustWork=T)
    
    ## Getting an allele list from the bed file
    haplo_allele_list <- get_haplo_allele_list(bed_filepath)
    
    ## Formatting the allele list to get a locus list
    haplo_locus_list <- tstrsplit(haplo_allele_list, '*', fixed=T)[[1]]
    haplo_locus_list <- tstrsplit(haplo_locus_list, 'KIR')[[2]]
    
    ## Find VCF filepath
    vcf_filepath <- file.path(vcf_directory, paste0(sample_name, '.', iteration_dir, '.vcf'))
    vcf_filepath <- normalizePath(vcf_filepath, mustWork=T)
    
    ## Read in VCF
    vcf_frame_raw <- read.table(vcf_filepath, header = F, stringsAsFactors = F)
    vcf_frame <- data.frame(vcf_frame_raw)
    
    ## Loop through every locus
    
    for(current_locus in haplo_locus_list){
      if(current_locus == '2DL5A' | current_locus == '2DL5B'){
        old_locus <- current_locus
        current_locus <- '2DL5'
      }
      
      if(iteration_dir == 'reference_1'){
        all_allele_calling_frame_list[[current_locus]] <- list()
      }
      
      ## Allele formatting
      allele_name <- haplo_allele_list[grep(current_locus, haplo_allele_list)]
      
      ## Creating current_locus allele frame (dataframe with all alleles)
      current_locus_allele_frame <- msf_to_allele_frame(msf_directory, current_locus)
      
      ## Create an allele specific bed conversion list to be used for VCF and allele frame naming
      list_allele_pos_conv <- alternate_make_bed_to_pos_conv(bed_filepath, current_locus_allele_frame, allele_name)
      
      ## Create a relative allele frame that only contains nucleotide positions for the given reference allele (deletion positions are removed)
      relative_allele_frame <- make_relative_allele_frame(current_locus_allele_frame, allele_name, list_allele_pos_conv)
      
      ## Create a unique position frame based on the relative allele frame
      unique_pos_frame <- allele_frame_to_unique_pos_frame(relative_allele_frame, allele_name)
      
      ## Subset VCF frame to only include current locus
      vcf_allele_frame <- vcf_frame[grep(allele_name, vcf_frame[,'V1'], fixed=T),]
      
      if(nrow(vcf_allele_frame) == 0){
        cat('\n\nNo aligned reads for this iteration and locus.\n\n')
        next
      }
      
      ## Removing Indels
      no_dels_vcf_allele_frame <- vcf_allele_frame[grep('INDEL', vcf_allele_frame[,'V8'], invert=T),]
      row.names(no_dels_vcf_allele_frame) <- paste0('X',no_dels_vcf_allele_frame[,'V2'])
      
      ## Match up row names of the vcf_allele_frame with conversion list
      row.names(no_dels_vcf_allele_frame) <- list_allele_pos_conv[row.names(no_dels_vcf_allele_frame)]
      
      ## Create a frame with only important positions
      #allele_calling_frame <- na.omit(no_dels_vcf_allele_frame[colnames(unique_pos_frame),])
      ## Instead we will return a frame with all positions
      allele_calling_frame <- no_dels_vcf_allele_frame
      
      ## Adding the allele_frame to the return list
      all_allele_calling_frame_list[[current_locus]][[iteration_dir]] <- allele_calling_frame
    }
  }
  return(all_allele_calling_frame_list)
}

## This function reads in the all KIR alleles aligned fasta file, outputs a list of kir alleles
read.all_kir_aligned_fasta <- function(fastaPath){
  ## Make sure the fasta path is actually a file
  fastaPath <- normalizePath(fastaPath, mustWork=T)
  
  ## Read in the fasta file
  fastaFile <- file(fastaPath, open='r')
  fileLines <- readLines(fastaFile)
  
  ## Initialize a list to store the alleles
  alleleList <- list()
  
  ## Iterate over each line in the file
  for(i in 1:length(fileLines)){
    ## split each line by 'IPD', this should isolate allele name lines
    currentLine <- strsplit(fileLines[i], 'IPD')[[1]]
    
    ## If the current line starts with '>', pull out the allele name and inialize a new list element
    if(currentLine[1] == '>'){
      alleleName <- strsplit(currentLine[2],' ')[[1]][2]
      alleleList[[alleleName]] <- ''
      
      ## Otherwise add the current line to the previous list element
    }else{
      #print(currentLine)
      alleleList[[alleleName]] <- paste0(alleleList[[alleleName]],currentLine)
    }
  }
  
  ## Close out the fasta file
  close(fastaFile)
  
  ## If the max allele length and min allele lenght are not the same throw an error
  if(max(unlist(lapply(alleleList, nchar))) != min(unlist(lapply(alleleList, nchar)))){
    stop('KIR allele lengths from the aligned fasta are not all the same.')
  }
  
  ## Return the completed list!
  return(alleleList)
}

## Speed-up stuff
haplo_resources_directory <- 'Resources/haplo_resources'
ipd_kir_resources_directory <- 'Resources/ipdkir_resources'
caller_resources_directory <- 'Resources/caller_resources'
msf_directory <- 'Resources/ipdkir_resources/IPD_KIR_MSF_Files/Formatted_MSF_Files'

kir_gen_path <- file.path(ipd_kir_resources_directory, 'genomic_sequence', 'KIR_gen.fasta')
kir_nuc_path <- file.path(ipd_kir_resources_directory, 'genomic_sequence', 'KIR_nuc.fasta')
all_aligned_kir_fasta_path <- file.path(ipd_kir_resources_directory, 'All_KIR_CDS_9-4-18.fas')


