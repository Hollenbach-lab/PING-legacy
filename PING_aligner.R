# Copyright 2017 Wesley Marin, Jill Hollenbach, Paul Norman, Ravi Dandekar
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

ping_aligner_version <- '1.0'
cat(paste0('PING_aligner version: ', ping_aligner_version))

ping_locus_aligner <- function(
  sample.location='PING_sequences/',
  fastq.pattern.1 = "_1.fastq.gz",
  fastq.pattern.2 = "_2.fastq.gz",
  bowtie.threads = 4,
  supported.loci = c("2DL1", "2DL23", "2DL4", "2DL5", "2DS3", "2DS4", "2DS5", "2DP1", "3DL1", "3DS1", "3DL2", "3DL3"),
  results.directory = '',
  kff.output = ''
){
  source("Resources/locus_functions.R", local = TRUE)
  library(data.table) ## Used for fread
  
  # Creates results directory, defaults to Aligner_results
  create_results_directory <- function() {
    cat("----- Getting PING ready -----\n\n")
    cat(paste0("Current working directory: ", getwd(), '\n\n'))
    
    if(results.directory != ""){
      save_to <- results.directory
    }else{
      save_to <- paste0("Aligner_results_", format(Sys.time(), "%Y_%m_%d_%H_%M"), "/")
      
      count <- 1
      while(file.exists(save_to)) {
        save_to <- paste0("Aligner_results_", format(Sys.time(), "%Y_%m_%d_%H_%M"), "_", count, "/")
        count <- count + 1
      }
    }
    dir.create(file.path(save_to), showWarnings = F)
    save_to <- normalizePath(save_to)
    
    cat(paste("Results being saved to", save_to, "\n\n"))
    
    return(save_to)
  }
  
  # Creates results folders
  ping_ready <- function(sample_list, results_directory) {
    
    dir.create(file.path(results_directory), showWarnings = F)
    dir.create(file.path(results_directory, "Vcf"), showWarnings = F)
    dir.create(file.path(results_directory, "Fastq"), showWarnings = F)
    
    for(sample in sample_list){
      dir.create(file.path(results_directory, 'Vcf', sample), showWarnings=F)
      dir.create(file.path(results_directory, 'Vcf', sample, 'locus'), showWarnings=F)
      dir.create(file.path(results_directory, 'Fastq', sample), showWarnings=F)
      dir.create(file.path(results_directory, 'Fastq', sample, 'locus'), showWarnings=F)
    }
    cat("\nResults subdirectories created.\n\n")
  }
  
  # Finds sequences
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
  
  # Reads in raw_kff_counts.csv table
  get_kff_probe_hits <- function() {
    if(kff.output == ''){
      kff_directory <- results_directory
    }else{
      kff_directory <- kff.output
    }
    
    cat('\nLooking for raw_kff_counts in: ', kff_directory)
    kff_directory <- normalizePath(kff_directory)
    
    kff_file <- file.path(kff_directory, 'raw_kff_counts.csv')
    
    if(file.exists(kff_file)){
      cat('\nFound: ', kff_file)
      kff_probe_table <- read.csv(kff_file, check.names = FALSE, stringsAsFactors = FALSE)
    }else{
      string <- paste('\nraw_kff_counts.csv was not found. Please make sure it exists.')
      stop(string)
    }
    return(kff_probe_table)
  }
  
  # Reads in kff_results.csv
  get_kff_locus_results <- function(){
    if(kff.output == ''){
      kff_directory <- results_directory
    }else{
      kff_directory <- kff.output
    }
    
    cat('\nLooking for kff_results in: ', kff_directory)
    kff_directory <- normalizePath(kff_directory)
    
    kff_file <- file.path(kff_directory, 'kff_results.csv')
    
    if(file.exists(kff_file)){
      cat('\nFound: ', kff_file)
      kff_probe_table <- read.csv(kff_file, check.names = FALSE, stringsAsFactors = FALSE)
    }else{
      string <- paste('\nkff_results.csv was not found. Please make sure it exists.')
      stop(string)
    }
    return(kff_probe_table)
  }
  
  # Checks if the matches of a specific prome meet the threshold
  kff_probe_check <- function(probe_table, probe_name, threshold=10){
    for_probe <- probe_name
    rev_probe <- paste0(probe_name,"_rc")
    
    probe_result <- ((probe_table[,for_probe] + probe_table[,rev_probe]) >= threshold)*1
    return(probe_result)
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
    
    if("2DS3" == current.locus) {
      KIR_2DS3(sample, is_gz)
    }
    
    if("2DS35" == current.locus){
      KIR_2DS35(sample, is_gz)
    }
    
    if("2DS4" == current.locus) {
      KIR_2DS4(sample, is_gz)
    }
    
    if("3DL1S1" == current.locus){
      KIR_3DL1S1(sample, is_gz)
    }
    
    if("3DL1" == current.locus){
      KIR_3DL1(sample, is_gz)
    }
    
    if("3DS1" == current.locus){
      KIR_3DS1(sample, is_gz)
    }
    
    if("3DL2" == current.locus) {
      KIR_3DL2(sample, is_gz)
    }
    
    if("3DL3" == current.locus) {
      KIR_3DL3(sample, is_gz)
    }
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
  
  results_directory <- create_results_directory()
  sample_list <- get_sample_list()
  ping_ready(sample_list, results_directory)
  kff_probe_table <- get_kff_probe_hits()
  kff_locus_table <- get_kff_locus_results()

  if('3DL1' %in% supported.loci || '3DS1' %in% supported.loci){
    supported.loci <- c(supported.loci, '3DL1S1')
  }
  
  if('2DS3' %in% supported.loci || '2DS5' %in% supported.loci){
    supported.loci <- c(supported.loci, '2DS35')
  }
  
  cat("\n\n---- Running PING locus alignment ----\n\n")
  for(sample in sample_list){
    loci_list <- kff_locus_table[sample,]
    probe_list <- kff_probe_table[sample,]
    for(current.locus in supported.loci){
      ping.locus(sample, current.locus, TRUE)
    }
    file.remove(list.files(pattern=sample, recursive = FALSE))
  }
  cat("\n\n---- Finished with PING locus alignment ----")
}