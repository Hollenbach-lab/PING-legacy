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

# VCF creation functions
KIR_2DL1 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is.gz = is_gz)
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  if(!"2DS1" %in% loci.list){
    
    # First pos filter
    bt2_fastq <- "-q"
    bt2_thread_parameter <- paste0("-p", bowtie.threads)
    bt2_5 <- "-5 3"
    bt2_3 <- "-3 7"
    bt2_L <- "-L 20"
    bt2_i <- "-i S,1,0.5"
    bt2_score_min <- "--score-min L,0,-0.187"
    bt2_I <- "-I 75"
    bt2_X <- "-X 1000"
    bt2_index <- "-x Resources/caller_resources/Filters/2DL1/All2DL1S1gen"
    bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
    bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
    bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL1S1in.fastq")
    bt2_sam <- paste0("-S ", sequence, ".temp")
    
    bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
    cat("\n", bowtie2_var, "\n")
    system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
    cat("\n")
    
    # First neg filter
    bt2_5 <- "-5 3"
    bt2_3 <- "-3 7"
    bt2_score_min <- "--score-min L,0,-0.155"
    bt2_index <- "-x Resources/caller_resources/Filters/2DL1/not2DL1S1"
    bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DL1S1in.1.fastq")
    bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DL1S1in.2.fastq")
    bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL1.fastq")
    bt2_sam <- paste0("-S ", sequence, ".temp")
    
    bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
    cat(bowtie2_var,"\n")
    system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
    cat("\n")
    
    # Final alignment
    bt2_score_min <- "--score-min L,0,-0.55"
    bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DL1AClong"
    bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1.1.fastq")
    bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1.2.fastq")
    bt2_sam <- paste0("-S ",sequence,"_2DL1.sam")
    bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL1.fastq")
    bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
    
    bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
    cat(bowtie2_var,"\n")
    system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
    cat("\n")
    
  }else{
    kff_results <- ping.kff(sequence, "2DL1", "Resources/caller_resources/KFF_2DL1.txt", 10)
    if((kff_results["2DL14and7"] == 1) || (kff_results["2DL14710b"] == 1)){
      
      # First pos filter
      bt2_fastq <- "-q"
      bt2_thread_parameter <- paste0("-p", bowtie.threads)
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_L <- "-L 20"
      bt2_i <- "-i S,1,0.5"
      bt2_score_min <- "--score-min L,0,-0.187"
      bt2_I <- "-I 75"
      bt2_X <- "-X 1000"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/All2DL1S1gen"
      bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
      bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
      bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL1S1in.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
      cat("\n", bowtie2_var, "\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
      cat("\n")
      
      # First neg filter
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_score_min <- "--score-min L,0,-0.155"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/not2DL1S1"
      bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DL1S1in.1.fastq")
      bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DL1S1in.2.fastq")
      bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL1.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
      cat("\n")
      
      # Second neg filter
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_score_min <- "--score-min L,0,-0.153"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DS1gen"
      bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1.1.fastq")
      bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1.2.fastq")
      bt2_un_conc <- paste0("--un-conc ", sequence,"_KIR2DL1b.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
      cat("\n")
      
      # Final alignment
      bt2_score_min <- "--score-min L,0,-0.55"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DL1AClong"
      bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1b.1.fastq")
      bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1b.2.fastq")
      bt2_sam <- paste0("-S ",sequence,"_2DL1.sam")
      bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL1.fastq")
      bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
      cat("\n")
      
      file.remove(paste0(sequence, "_KIR2DL1b.1.fastq"))
      file.remove(paste0(sequence, "_KIR2DL1b.2.fastq"))
    }else if((kff_results["2DL14and7"] == 0) && (kff_results["2DL14710b"] == 0)){
      
      # Pos filter
      bt2_fastq <- "-q"
      bt2_thread_parameter <- paste0("-p", bowtie.threads)
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_L <- "-L 20"
      bt2_i <- "-i S,1,0.5"
      bt2_score_min <- "--score-min L,0,-0.187"
      bt2_I <- "-I 75"
      bt2_X <- "-X 1000"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/All2DL1S1gen"
      bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
      bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
      bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL1S1in.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
      cat("\n", bowtie2_var, "\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
      cat("\n")
      
      # First neg filter
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_score_min <- "--score-min L,0,-0.155"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/not2DL1S1"
      bt2_sequence_1 <- paste0("-1 ","r", sequence,"_2DL1S1in.1.fastq")
      bt2_sequence_2 <- paste0("-2 ","r", sequence,"_2DL1S1in.2.fastq")
      bt2_un_conc <- paste0("--un-conc ", sequence,"_KIR2DL1.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
      cat("\n")
      
      # Second neg filter
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_score_min <- "--score-min L,0,-0.18"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DS1gen"
      bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1.1.fastq")
      bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1.2.fastq")
      bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL1b.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
      cat("\n")
      
      # Final alignment
      bt2_score_min <- "--score-min L,0,-0.55"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DL1AClong"
      bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1b.1.fastq")
      bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1b.2.fastq")
      bt2_sam <- paste0("-S ",sequence,"_2DL1.sam")
      bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL1.fastq")
      bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
      cat("\n")
      
      file.remove(paste0(sequence, "_KIR2DL1b.1.fastq"))
      file.remove(paste0(sequence, "_KIR2DL1b.2.fastq"))
    }
  }
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DL1.sam")
  st_out <- paste0("-o ",sequence,"_2DL1.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DL1.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DL1.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DL1/2DL1AClong.fas"
  st_in <- paste0(sequence,"_2DL1.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DL1/2DL1AClong.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DL1nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DL1cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DL1.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL1.1.fastq"))
  #file.copy(paste0(sequence, "_2DL1.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL1.2.fastq"))
  
  file.copy(paste0(sequence, "_2DL1nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DL1nuc.vcf"))
  
  file.remove(paste0(sequence, "_2DL1nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DL1S1in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DL1S1in.2.fastq"))
  file.remove(paste0(sequence,"_2DL1.1.fastq"))
  file.remove(paste0(sequence,"_2DL1.2.fastq"))
  file.remove(paste0(sequence,"_2DL1.sam"))
  file.remove(paste0(sequence,"_2DL1.bam"))
  file.remove(paste0(sequence,"_2DL1.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DL1.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DL1.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}


KIR_2DL1_Genomic <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is.gz = is_gz)
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  if(!"2DS1" %in% loci.list){
    
    # First pos filter
    bt2_fastq <- "-q"
    bt2_thread_parameter <- paste0("-p", bowtie.threads)
    bt2_5 <- "-5 3"
    bt2_3 <- "-3 7"
    bt2_L <- "-L 20"
    bt2_i <- "-i S,1,0.5"
    bt2_score_min <- "--score-min L,0,-0.187"
    bt2_I <- "-I 75"
    bt2_X <- "-X 1000"
    bt2_index <- "-x Resources/caller_resources/Filters/2DL1/All2DL1S1gen"
    bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
    bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
    bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL1S1in.fastq")
    bt2_sam <- paste0("-S ", sequence, ".temp")
    
    bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
    cat("\n", bowtie2_var, "\n")
    system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
    cat("\n")
    
    # First neg filter
    bt2_5 <- "-5 3"
    bt2_3 <- "-3 7"
    bt2_score_min <- "--score-min L,0,-0.155"
    bt2_index <- "-x Resources/caller_resources/Filters/2DL1/not2DL1S1"
    bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DL1S1in.1.fastq")
    bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DL1S1in.2.fastq")
    bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL1.fastq")
    bt2_sam <- paste0("-S ", sequence, ".temp")
    
    bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
    cat(bowtie2_var,"\n")
    system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
    cat("\n")
    
    # Final alignment
    bt2_score_min <- "--score-min L,0,-0.55"
    bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DL1AClong"
    bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1.1.fastq")
    bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1.2.fastq")
    bt2_sam <- paste0("-S ",sequence,"_2DL1.sam")
    bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL1.fastq")
    bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
    
    bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
    cat(bowtie2_var,"\n")
    system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
    cat("\n")
    
  }else{
    kff_results <- ping.kff(sequence, "2DL1", "Resources/caller_resources/KFF_2DL1.txt", 10)
    if((kff_results["2DL14and7"] == 1) || (kff_results["2DL14710b"] == 1)){
      
      # First pos filter
      bt2_fastq <- "-q"
      bt2_thread_parameter <- paste0("-p", bowtie.threads)
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_L <- "-L 20"
      bt2_i <- "-i S,1,0.5"
      bt2_score_min <- "--score-min L,0,-0.187"
      bt2_I <- "-I 75"
      bt2_X <- "-X 1000"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/All2DL1S1gen"
      bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
      bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
      bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL1S1in.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
      cat("\n", bowtie2_var, "\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
      cat("\n")
      
      # First neg filter
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_score_min <- "--score-min L,0,-0.155"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/not2DL1S1"
      bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DL1S1in.1.fastq")
      bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DL1S1in.2.fastq")
      bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL1.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
      cat("\n")
      
      # Second neg filter
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_score_min <- "--score-min L,0,-0.153"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DS1gen"
      bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1.1.fastq")
      bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1.2.fastq")
      bt2_un_conc <- paste0("--un-conc ", sequence,"_KIR2DL1b.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
      cat("\n")
      
      # Final alignment
      bt2_score_min <- "--score-min L,0,-0.55"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DL1AClong"
      bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1b.1.fastq")
      bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1b.2.fastq")
      bt2_sam <- paste0("-S ",sequence,"_2DL1.sam")
      bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL1.fastq")
      bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
      cat("\n")
      
      file.remove(paste0(sequence, "_KIR2DL1b.1.fastq"))
      file.remove(paste0(sequence, "_KIR2DL1b.2.fastq"))
    }else if((kff_results["2DL14and7"] == 0) && (kff_results["2DL14710b"] == 0)){
      
      # Pos filter
      bt2_fastq <- "-q"
      bt2_thread_parameter <- paste0("-p", bowtie.threads)
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_L <- "-L 20"
      bt2_i <- "-i S,1,0.5"
      bt2_score_min <- "--score-min L,0,-0.187"
      bt2_I <- "-I 75"
      bt2_X <- "-X 1000"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/All2DL1S1gen"
      bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
      bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
      bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL1S1in.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
      cat("\n", bowtie2_var, "\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
      cat("\n")
      
      # First neg filter
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_score_min <- "--score-min L,0,-0.155"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/not2DL1S1"
      bt2_sequence_1 <- paste0("-1 ","r", sequence,"_2DL1S1in.1.fastq")
      bt2_sequence_2 <- paste0("-2 ","r", sequence,"_2DL1S1in.2.fastq")
      bt2_un_conc <- paste0("--un-conc ", sequence,"_KIR2DL1.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
      cat("\n")
      
      # Second neg filter
      bt2_5 <- "-5 3"
      bt2_3 <- "-3 7"
      bt2_score_min <- "--score-min L,0,-0.18"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DS1gen"
      bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1.1.fastq")
      bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1.2.fastq")
      bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL1b.fastq")
      bt2_sam <- paste0("-S ", sequence, ".temp")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
      cat("\n")
      
      # Final alignment
      bt2_score_min <- "--score-min L,0,-0.55"
      bt2_index <- "-x Resources/caller_resources/Filters/2DL1/2DL1AClong"
      bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL1b.1.fastq")
      bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL1b.2.fastq")
      bt2_sam <- paste0("-S ",sequence,"_2DL1.sam")
      bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL1.fastq")
      bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
      
      bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
      cat(bowtie2_var,"\n")
      system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
      cat("\n")
      
      file.remove(paste0(sequence, "_KIR2DL1b.1.fastq"))
      file.remove(paste0(sequence, "_KIR2DL1b.2.fastq"))
    }
  }
  
  # Filter out contaminating reads from SAM file
  filter_contam_reads_from_sam_file(paste0(sequence,"_2DL1.sam"), sequence)
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DL1.sam")
  st_out <- paste0("-o ",sequence,"_2DL1.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DL1.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DL1.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DL1/2DL1AClong.fas"
  st_in <- paste0(sequence,"_2DL1.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DL1/2DL1AClong.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DL1nuc.vcf")
  
  cat("\n")
  #st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  #system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DL1cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DL1.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL1.1.fastq"))
  #file.copy(paste0(sequence, "_2DL1.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL1.2.fastq"))
  
  file.copy(paste0(sequence, "_2DL1nuc.vcf"), paste0(results.directory, "Vcf_Genomic/", sequence, "_2DL1nuc.vcf"))
  file.copy(paste0(sequence, "_2DL1.sam"), paste0(results.directory, "Sam_Genomic/", sequence, "_2DL1.sam"))
  
  file.remove(paste0(sequence, "_2DL1nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DL1S1in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DL1S1in.2.fastq"))
  file.remove(paste0(sequence,"_2DL1.1.fastq"))
  file.remove(paste0(sequence,"_2DL1.2.fastq"))
  file.remove(paste0(sequence,"_2DL1.sam"))
  file.remove(paste0(sequence,"_2DL1.bam"))
  file.remove(paste0(sequence,"_2DL1.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DL1.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DL1.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}


KIR_2DL23 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  
  ## Positive filter
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL23/All2DL23"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL23in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  ## Negative filter
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 10"
  bt2_score_min <- "--score-min L,0,-0.135"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL23/Not2DL23b"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DL23in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DL23in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL23.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  ## Removes everything that could be 2DS2
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 35"
  bt2_score_min <- "--score-min L,0,-0.09"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL23/All2DS2"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL23.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL23.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL23_notS2.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  ## Remove 2DL3 tail
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL23/2DL3long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL23_notS2.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL23_notS2.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_2DL2_notL3tail.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  ## 2DL2 alignment
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL23/2DL2long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_2DL2_notL3tail.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_2DL2_notL3tail.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DL2.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL2.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  ## 2DL2 VCF generation
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DL2.sam")
  st_out <- paste0("-o ",sequence,"_2DL2.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DL2.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DL2.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DL23/2DL2long.fas"
  st_in <- paste0(sequence,"_2DL2.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DL23/2DL2longtail.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DL2nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DL2cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DL2.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL2.1.fastq"))
  #file.copy(paste0(sequence, "_2DL2.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL2.2.fastq"))
  
  file.copy(paste0(sequence, "_2DL2nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DL2nuc.vcf"))
  
  
  ## 2DL3 alignment
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL23/2DL3long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL23_notS2.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL23_notS2.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DL23.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL23.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  ## 2DL23 VCF generation
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DL23.sam")
  st_out <- paste0("-o ",sequence,"_2DL23.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DL23.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DL23.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DL23/2DL3long.fas"
  st_in <- paste0(sequence,"_2DL23.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DL23/2DL3longshort.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DL23nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DL23cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DL23.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL23.1.fastq"))
  #file.copy(paste0(sequence, "_2DL23.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL23.2.fastq"))
  
  file.copy(paste0(sequence, "_2DL23nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DL23nuc.vcf"))
  
  
  ## Remove 2DL2 tail
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL23/2DL2long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL23_notS2.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL23_notS2.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_2DL3_notL2tail.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  ## 2DL3 alignment
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL23/2DL3long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_2DL3_notL2tail.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_2DL3_notL2tail.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DL3tail.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL3.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  ## 2DL3 VCF generation
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DL3tail.sam")
  st_out <- paste0("-o ",sequence,"_2DL3tail.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DL3tail.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DL3tail.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DL23/2DL3long.fas"
  st_in <- paste0(sequence,"_2DL3tail.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DL23/2DL3longtail.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ",sequence, "_2DL3nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DL3cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DL3.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL3.1.fastq"))
  #file.copy(paste0(sequence, "_2DL3.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL3.2.fastq"))
  
  file.copy(paste0(sequence, "_2DL3nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DL3nuc.vcf"))
  
  
  ## Cleanup
  file.remove(paste0(sequence, "_2DL2nuc.vcf"))
  file.remove(paste0(sequence, "_2DL23nuc.vcf"))
  file.remove(paste0(sequence, "_2DL3nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DL23in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DL23in.2.fastq"))
  file.remove(paste0(sequence,"_2DL2.1.fastq"))
  file.remove(paste0(sequence,"_2DL2.2.fastq"))
  file.remove(paste0(sequence,"_2DL23.1.fastq"))
  file.remove(paste0(sequence,"_2DL23.2.fastq"))
  file.remove(paste0(sequence,"_2DL3.1.fastq"))
  file.remove(paste0(sequence,"_2DL3.2.fastq"))
  file.remove(paste0(sequence,"_2DL2.sam"))
  file.remove(paste0(sequence,"_2DL23.sam"))
  file.remove(paste0(sequence,"_2DL3tail.sam"))
  file.remove(paste0(sequence,"_2DL2.bam"))
  file.remove(paste0(sequence,"_2DL23.bam"))
  file.remove(paste0(sequence,"_2DL3tail.bam"))
  file.remove(paste0(sequence,"_2DL2.sorted.bam"))
  file.remove(paste0(sequence,"_2DL23.sorted.bam"))
  file.remove(paste0(sequence,"_2DL3tail.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DL23.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DL23.2.fastq"))
  file.remove(paste0(sequence,"_KIR2DL23_notS2.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DL23_notS2.2.fastq"))
  file.remove(paste0(sequence,"_2DL2_notL3tail.1.fastq"))
  file.remove(paste0(sequence,"_2DL2_notL3tail.2.fastq"))
  file.remove(paste0(sequence,"_2DL3_notL2tail.1.fastq"))
  file.remove(paste0(sequence,"_2DL3_notL2tail.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
  
}

KIR_2DL4 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL4/All2DL4gen"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL4in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.2"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL4/not2DL4"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DL4in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DL4in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL4.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL4/2DL4FH5"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL4.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL4.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DL4.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL4.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DL4.sam")
  st_out <- paste0("-o ",sequence,"_2DL4.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DL4.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DL4.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DL4/2DL4FH5.fas"
  st_in <- paste0(sequence,"_2DL4.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DL4/2DL4FH5.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DL4nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DL4cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DL4.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL4.1.fastq"))
  #file.copy(paste0(sequence, "_2DL4.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL4.2.fastq"))
  
  file.copy(paste0(sequence, "_2DL4nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DL4nuc.vcf"))
  
  file.remove(paste0(sequence, "_2DL4nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DL4in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DL4in.2.fastq"))
  file.remove(paste0(sequence,"_2DL4.1.fastq"))
  file.remove(paste0(sequence,"_2DL4.2.fastq"))
  file.remove(paste0(sequence,"_2DL4.sam"))
  file.remove(paste0(sequence,"_2DL4.bam"))
  file.remove(paste0(sequence,"_2DL4.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DL4.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DL4.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}

KIR_2DL5 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 280000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL5/All2DL5"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DL5in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.1"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL5/not2DL5"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DL5in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DL5in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DL5.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DL5/2DL5B"
  #    bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DL5.1.fastq")
  #    bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DL5.2.fastq")
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DL5in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DL5in.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DL5.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DL5.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DL5.sam")
  st_out <- paste0("-o ",sequence,"_2DL5.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DL5.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DL5.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DL5/2DL5B.fas"
  st_in <- paste0(sequence,"_2DL5.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DL5/2DL5B.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DL5nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DL5cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DL5.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL5.1.fastq"))
  #file.copy(paste0(sequence, "_2DL5.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DL5.2.fastq"))
  
  file.copy(paste0(sequence, "_2DL5nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DL5nuc.vcf"))
  
  file.remove(paste0(sequence, "_2DL5nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DL5in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DL5in.2.fastq"))
  file.remove(paste0(sequence,"_2DL5.1.fastq"))
  file.remove(paste0(sequence,"_2DL5.2.fastq"))
  file.remove(paste0(sequence,"_2DL5.sam"))
  file.remove(paste0(sequence,"_2DL5.bam"))
  file.remove(paste0(sequence,"_2DL5.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DL5.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DL5.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}

KIR_2DP1 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/2DP1/All2DP1gen"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DP1in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.2"
  bt2_index <- "-x Resources/caller_resources/Filters/2DP1/not2DP1"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DP1in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DP1in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DP1.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DP1/2DP1"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DP1.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DP1.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DP1.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DP1.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DP1.sam")
  st_out <- paste0("-o ",sequence,"_2DP1.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DP1.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DP1.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DP1/2DP1.fas"
  st_in <- paste0(sequence,"_2DP1.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DP1/2DP1.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DP1nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DP1cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DP1.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DP1.1.fastq"))
  #file.copy(paste0(sequence, "_2DP1.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DP1.2.fastq"))
  
  file.copy(paste0(sequence, "_2DP1nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DP1nuc.vcf"))
  
  file.remove(paste0(sequence, "_2DP1nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DP1in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DP1in.2.fastq"))
  file.remove(paste0(sequence,"_2DP1.1.fastq"))
  file.remove(paste0(sequence,"_2DP1.2.fastq"))
  file.remove(paste0(sequence,"_2DP1.sam"))
  file.remove(paste0(sequence,"_2DP1.bam"))
  file.remove(paste0(sequence,"_2DP1.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DP1.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DP1.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}

KIR_2DS3 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS35/All2DS35_Gen"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DS35in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.153"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS35/not2DS35"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DS35in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DS35in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DS35.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS35/2DS3"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DS35.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DS35.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DS3.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DS3.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DS3.sam")
  st_out <- paste0("-o ",sequence,"_2DS3.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DS3.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DS3.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DS35/2DS3.fas"
  st_in <- paste0(sequence,"_2DS3.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DS35/2DS3.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DS3nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DS3cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DS3.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DS3.1.fastq"))
  #file.copy(paste0(sequence, "_2DS3.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DS3.2.fastq"))
  
  file.copy(paste0(sequence, "_2DS3nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DS3nuc.vcf"))
  
  file.remove(paste0(sequence, "_2DS3nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DS35in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DS35in.2.fastq"))
  file.remove(paste0(sequence,"_2DS3.1.fastq"))
  file.remove(paste0(sequence,"_2DS3.2.fastq"))
  file.remove(paste0(sequence,"_2DS3.sam"))
  file.remove(paste0(sequence,"_2DS3.bam"))
  file.remove(paste0(sequence,"_2DS3.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DS35.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DS35.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}

KIR_2DS4 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS4/All2DS4"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DS4in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.1"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS4/not2DS4"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DS4in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DS4in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DS4.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  bt2_score_min <- "--score-min L,0,-0.17"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS4/2DS4long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DS4.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DS4.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DS4.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DS4.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DS4.sam")
  st_out <- paste0("-o ",sequence,"_2DS4.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DS4.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DS4.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DS4/2DS4long.fas"
  st_in <- paste0(sequence,"_2DS4.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DS4/2DS4noE5.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DS4nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DS4cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DS4.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DS4.1.fastq"))
  #file.copy(paste0(sequence, "_2DS4.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DS4.2.fastq"))
  
  file.copy(paste0(sequence, "_2DS4nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DS4nuc.vcf"))
  
  file.remove(paste0(sequence, "_2DS4nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DS4in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DS4in.2.fastq"))
  file.remove(paste0(sequence,"_2DS4.1.fastq"))
  file.remove(paste0(sequence,"_2DS4.2.fastq"))
  file.remove(paste0(sequence,"_2DS4.sam"))
  file.remove(paste0(sequence,"_2DS4.bam"))
  file.remove(paste0(sequence,"_2DS4.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DS4.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DS4.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}

KIR_2DS35 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 240000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS35/All2DS35_Gen"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_2DS35in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.18"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS35/not2DS35"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_2DS35in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_2DS35in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR2DS35.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  bt2_score_min <- "--score-min L,0,-0.5"
  bt2_index <- "-x Resources/caller_resources/Filters/2DS35/2DS5long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR2DS35.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR2DS35.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_2DS35.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_2DS35.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_2DS35.sam")
  st_out <- paste0("-o ",sequence,"_2DS35.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_2DS35.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_2DS35.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/2DS35/2DS5long.fas"
  st_in <- paste0(sequence,"_2DS35.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/2DS35/2DS5long.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_2DS35nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call, bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_2DS35cons.fastq")))
  
  #file.copy(paste0(sequence, "_2DS35.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DS35.1.fastq"))
  #file.copy(paste0(sequence, "_2DS35.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_2DS35.2.fastq"))
  
  file.copy(paste0(sequence, "_2DS35nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_2DS35nuc.vcf"))
  
  file.remove(paste0(sequence, "_2DS35nuc.vcf"))
  file.remove(paste0("r",sequence,"_2DS35in.1.fastq"))
  file.remove(paste0("r",sequence,"_2DS35in.2.fastq"))
  file.remove(paste0(sequence,"_2DS35.1.fastq"))
  file.remove(paste0(sequence,"_2DS35.2.fastq"))
  file.remove(paste0(sequence,"_2DS35.sam"))
  file.remove(paste0(sequence,"_2DS35.bam"))
  file.remove(paste0(sequence,"_2DS35.sorted.bam"))
  file.remove(paste0(sequence,"_KIR2DS35.1.fastq"))
  file.remove(paste0(sequence,"_KIR2DS35.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}

KIR_3DL1 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/All3DL1andS1"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_3DL12in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.2"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/not3DL1S1"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_3DL12in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_3DL12in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR3DL1S1.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  
  ## 3DL1
  bt2_local <- "--local"
  bt2_N <- "-N 1"
  bt2_score_min <- "--score-min L,0,0.6"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/3DL1longii"
  bt2_t <- "-t"
  bt2_X <- "-X 750"
  bt2_no_u <- "--no-unal"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR3DL1S1.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR3DL1S1.2.fastq")
  bt2_sam <- paste0("-S ", sequence, "_3DL1.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_3DL1.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_local, bt2_N, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_no_u, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_local, bt2_N, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_no_u, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  
  ## 3DL1
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_3DL1.sam")
  st_out <- paste0("-o ",sequence,"_3DL1.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ", sequence, "_3DL1.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_3DL1.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/3DL1/3DL1longii.fas"
  st_in <- paste0(sequence,"_3DL1.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/3DL1/3DL1longii.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_3DL1nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call,bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_3DL1cons.fastq")))
  
  #file.copy(paste0(sequence, "_3DL1.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DL1.1.fastq"))
  #file.copy(paste0(sequence, "_3DL1.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DL1.2.fastq"))
  
  file.copy(paste0(sequence, "_3DL1nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_3DL1nuc.vcf"))
  
  
  ## Removing files
  
  file.remove(paste0("r",sequence,"_3DL12in.1.fastq"))
  file.remove(paste0("r",sequence,"_3DL12in.2.fastq"))
  
  file.remove(paste0(sequence,"_3DL1.sam"))
  file.remove(paste0(sequence,"_3DL1.bam"))
  file.remove(paste0(sequence,"_3DL1.sorted.bam"))
  
  file.remove(paste0(sequence, "_3DL1nuc.vcf"))
  
  file.remove(paste0(sequence, "_3DL1.1.fastq"))
  file.remove(paste0(sequence, "_3DL1.2.fastq"))
  
  file.remove(paste0(sequence,"_KIR3DL1S1.1.fastq"))
  file.remove(paste0(sequence,"_KIR3DL1S1.2.fastq"))
  
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}

KIR_3DL1S1 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/All3DL1andS1"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_3DL12in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.2"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/not3DL1S1"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_3DL12in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_3DL12in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR3DL1S1.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  
  ## 3DL1het
  bt2_local <- "--local"
  bt2_N <- "-N 1"
  bt2_score_min <- "--score-min L,0,0.6"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/3DL1longCAT"
  bt2_t <- "-t"
  bt2_X <- "-X 750"
  bt2_no_u <- "--no-unal"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR3DL1S1.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR3DL1S1.2.fastq")
  bt2_sam <- paste0("-S ", sequence, "_3DL1het.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_3DL1het.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_local, bt2_N, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_no_u, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_local, bt2_N, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_no_u, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  
  ## 3DS1het
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.17"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/not3DL1S1"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_3DL12in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_3DL12in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_KIR3DS1.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  
  bt2_local <- "--local"
  bt2_N <- "-N 1"
  bt2_score_min <- "--score-min L,0,0.6"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/3DS1longCAT"
  bt2_t <- "-t"
  bt2_X <- "-X 750"
  bt2_no_u <- "--no-unal"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR3DS1.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR3DS1.2.fastq")
  bt2_sam <- paste0("-S ", sequence, "_3DS1.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_3DS1.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_local, bt2_N, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_no_u, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_local, bt2_N, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_no_u, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  
  ## 3DL1het
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_3DL1het.sam")
  st_out <- paste0("-o ",sequence,"_3DL1het.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_3DL1het.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_3DL1het.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/3DL1/3DL1longCAT.fas"
  st_in <- paste0(sequence,"_3DL1het.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/3DL1/3DL1longCAT.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_3DL1hetnuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call,bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_3DL1hetcons.fastq")))
  
  #file.copy(paste0(sequence, "_3DL1het.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DL1het.1.fastq"))
  #file.copy(paste0(sequence, "_3DL1het.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DL1het.2.fastq"))
  
  file.copy(paste0(sequence, "_3DL1hetnuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_3DL1nuc.vcf"))
  
  
  ## 3DS1
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_3DS1.sam")
  st_out <- paste0("-o ",sequence,"_3DS1.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_3DS1.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_3DS1.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/3DL1/3DS1longCAT.fas"
  st_in <- paste0(sequence,"_3DS1.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/3DL1/3DS1longCAT.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_3DS1nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call,bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_3DS1cons.fastq")))
  
  #file.copy(paste0(sequence, "_3DS1.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DS1.1.fastq"))
  #file.copy(paste0(sequence, "_3DS1.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DS1.2.fastq"))
  
  file.copy(paste0(sequence, "_3DS1nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_3DS1nuc.vcf"))
  
  ## Consensus
  create_consensus <- function(sample, vcf.file, filter.directory, results.directory){
    
    create_consesus(paste0(results.directory, "Vcf/", sequence, "_3DS1nuc.vcf"))
  }
  
  ## Removing files
  file.remove(paste0("r",sequence,"_3DL12in.1.fastq"))
  file.remove(paste0("r",sequence,"_3DL12in.2.fastq"))
  
  file.remove(paste0(sequence,"_3DL1het.sam"))
  file.remove(paste0(sequence,"_3DL1het.bam"))
  file.remove(paste0(sequence,"_3DL1het.sorted.bam"))
  file.remove(paste0(sequence,"_3DS1.sam"))
  file.remove(paste0(sequence,"_3DS1.bam"))
  file.remove(paste0(sequence,"_3DS1.sorted.bam"))
  
  file.remove(paste0(sequence, "_3DL1hetnuc.vcf"))
  file.remove(paste0(sequence, "_3DS1nuc.vcf"))
  
  file.remove(paste0(sequence, "_3DL1het.1.fastq"))
  file.remove(paste0(sequence, "_3DL1het.2.fastq"))
  file.remove(paste0(sequence, "_3DS1.1.fastq"))
  file.remove(paste0(sequence, "_3DS1.2.fastq"))
  
  file.remove(paste0(sequence,"_KIR3DL1S1.1.fastq"))
  file.remove(paste0(sequence,"_KIR3DL1S1.2.fastq"))
  file.remove(paste0(sequence,"_KIR3DS1.1.fastq"))
  file.remove(paste0(sequence,"_KIR3DS1.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
  
  
}

KIR_3DS1 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/All3DL1andS1"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_3DL12in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.2"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/not3DL1S1"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_3DL12in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_3DL12in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR3DL1S1.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  
  ## 3DS1
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.17"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/not3DL1S1"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_3DL12in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_3DL12in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_KIR3DS1.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  
  bt2_local <- "--local"
  bt2_N <- "-N 1"
  bt2_score_min <- "--score-min L,0,0.6"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL1/3DS1longCAT"
  bt2_t <- "-t"
  bt2_X <- "-X 750"
  bt2_no_u <- "--no-unal"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR3DS1.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR3DS1.2.fastq")
  bt2_sam <- paste0("-S ", sequence, "_3DS1.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_3DS1.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_local, bt2_N, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_no_u, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_local, bt2_N, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_no_u, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  
  ## 3DS1
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_3DS1.sam")
  st_out <- paste0("-o ",sequence,"_3DS1.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_3DS1.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_3DS1.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/3DL1/3DS1longCAT.fas"
  st_in <- paste0(sequence,"_3DS1.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/3DL1/3DS1longCAT.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_3DS1nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call,bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_3DS1cons.fastq")))
  
  #file.copy(paste0(sequence, "_3DS1.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DS1.1.fastq"))
  #file.copy(paste0(sequence, "_3DS1.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DS1.2.fastq"))
  
  file.copy(paste0(sequence, "_3DS1nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_3DS1nuc.vcf"))
  
  
  ## Removing files
  file.remove(paste0("r",sequence,"_3DL12in.1.fastq"))
  file.remove(paste0("r",sequence,"_3DL12in.2.fastq"))
  
  file.remove(paste0(sequence,"_3DS1.sam"))
  file.remove(paste0(sequence,"_3DS1.bam"))
  file.remove(paste0(sequence,"_3DS1.sorted.bam"))
  
  file.remove(paste0(sequence, "_3DS1nuc.vcf"))
  
  file.remove(paste0(sequence, "_3DS1.1.fastq"))
  file.remove(paste0(sequence, "_3DS1.2.fastq"))
  
  file.remove(paste0(sequence,"_KIR3DL1S1.1.fastq"))
  file.remove(paste0(sequence,"_KIR3DL1S1.2.fastq"))
  file.remove(paste0(sequence,"_KIR3DS1.1.fastq"))
  file.remove(paste0(sequence,"_KIR3DS1.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
  
  
}

KIR_3DL2 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL2/All3DL2"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_3DL2in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.2"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL2/not3DL2"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_3DL2in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_3DL2in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR3DL2.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL2/3DL2long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR3DL2.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR3DL2.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_3DL2.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_3DL2.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_3DL2.sam")
  st_out <- paste0("-o ",sequence,"_3DL2.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_3DL2.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_3DL2.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/3DL2/3DL2long.fas"
  st_in <- paste0(sequence,"_3DL2.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/3DL2/3DL2long.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_3DL2nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call,bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_3DL2cons.fastq")))
  
  #file.copy(paste0(sequence, "_3DL2.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DL2.1.fastq"))
  #file.copy(paste0(sequence, "_3DL2.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DL2.2.fastq"))
  
  file.copy(paste0(sequence, "_3DL2nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_3DL2nuc.vcf"))
  
  file.remove(paste0(sequence, "_3DL2nuc.vcf"))
  file.remove(paste0("r",sequence,"_3DL2in.1.fastq"))
  file.remove(paste0("r",sequence,"_3DL2in.2.fastq"))
  file.remove(paste0(sequence,"_3DL2.1.fastq"))
  file.remove(paste0(sequence,"_3DL2.2.fastq"))
  file.remove(paste0(sequence,"_3DL2.sam"))
  file.remove(paste0(sequence,"_3DL2.bam"))
  file.remove(paste0(sequence,"_3DL2.sorted.bam"))
  file.remove(paste0(sequence,"_KIR3DL2.1.fastq"))
  file.remove(paste0(sequence,"_KIR3DL2.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}

KIR_3DL3 <- function(sequence, is_gz) {
  
  if(is_gz){
    fastq.pattern.1 <- unlist(strsplit(fastq.pattern.1, ".gz", fixed = TRUE))
    fastq.pattern.2 <- unlist(strsplit(fastq.pattern.2, ".gz", fixed = TRUE))
    
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2, ".gz"), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }else{
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.1), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1), is_gz)
    cut_fastq(paste0(sample.location, sequence, fastq.pattern.2), read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2), is_gz)
  }
  
  
  #sapply(paste0(sample.location, sequence, fastq.pattern.1), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.1))
  #sapply(paste0(sample.location, sequence, fastq.pattern.2), cut_fastq, read.cap = 120000, post.file.name = paste0(sequence, fastq.pattern.2))
  
  bt2_fastq <- "-q"
  bt2_thread_parameter <- paste0("-p", bowtie.threads)
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_L <- "-L 20"
  bt2_i <- "-i S,1,0.5"
  bt2_score_min <- "--score-min L,0,-0.187"
  bt2_I <- "-I 75"
  bt2_X <- "-X 1000"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL3/All3DL3"
  bt2_sequence_1 <- paste0("-1 ", sequence, fastq.pattern.1)
  bt2_sequence_2 <- paste0("-2 ", sequence, fastq.pattern.2)
  bt2_al_conc <- paste0("--al-conc ","r",sequence, "_3DL3in.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2",bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam)
  cat("\n", bowtie2_var, "\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_al_conc, bt2_sam))
  cat("\n")
  
  bt2_5 <- "-5 3"
  bt2_3 <- "-3 7"
  bt2_score_min <- "--score-min L,0,-0.2"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL3/not3DL3"
  bt2_sequence_1 <- paste0("-1 ","r",sequence,"_3DL3in.1.fastq")
  bt2_sequence_2 <- paste0("-2 ","r",sequence,"_3DL3in.2.fastq")
  bt2_un_conc <- paste0("--un-conc ",sequence,"_KIR3DL3.fastq")
  bt2_sam <- paste0("-S ", sequence, ".temp")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_un_conc, bt2_sam))
  cat("\n")
  
  bt2_score_min <- "--score-min L,0,-0.35"
  bt2_index <- "-x Resources/caller_resources/Filters/3DL3/3DL3long"
  bt2_sequence_1 <- paste0("-1 ", sequence,"_KIR3DL3.1.fastq")
  bt2_sequence_2 <- paste0("-2 ", sequence,"_KIR3DL3.2.fastq")
  bt2_sam <- paste0("-S ",sequence,"_3DL3.sam")
  bt2_al_conc <- paste0("--al-conc ", sequence, "_3DL3.fastq")
  bt2_un_conc <- paste0("--un-conc ", sequence, "_notmapped.fastq")
  
  bowtie2_var <- paste("bowtie2", bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc)
  cat(bowtie2_var,"\n")
  system2("bowtie2", c(bt2_fastq, bt2_thread_parameter, bt2_5, bt2_3, bt2_L, bt2_i, bt2_score_min, bt2_I, bt2_X, bt2_index, bt2_sequence_1, bt2_sequence_2, bt2_sam, bt2_al_conc, bt2_un_conc))
  cat("\n")
  
  st_b <- "-b"
  st_q <- "-q10"
  st_in <- paste0(sequence,"_3DL3.sam")
  st_out <- paste0("-o ",sequence,"_3DL3.bam")
  
  st_view <- paste("samtools view",st_b,st_q,st_in,st_out)
  cat(st_view,"\n")
  system2("samtools", c("view", st_b, st_q, st_in, st_out))
  cat("\n")
  
  st_out <- paste0("-o ",sequence,"_3DL3.sorted.bam")
  st_T <- paste0("-T temp")
  st_in <- paste0(sequence,"_3DL3.bam")
  
  st_sort <- paste("samtools sort",st_out,st_T,st_in)
  cat(st_sort,"\n")
  system2("samtools", c("sort", st_out, st_T, st_in))
  cat("\n")
  
  st_m <- "-m 3"
  st_F <- "-F 0.0002"
  st_u <- "-u"
  st_f <- "-f Resources/caller_resources/Filters/3DL3/3DL3long.fas"
  st_in <- paste0(sequence,"_3DL3.sorted.bam")
  st_l <- "-l Resources/caller_resources/Filters/3DL3/3DL3long.bed"
  st_break <- "|"
  bcf_call <- "bcftools call"
  bcf_multi_al <- "--multiallelic-caller"
  bcf_O <- "-O v"
  bcf_out <- paste0("-o ", sequence, "_3DL3nuc.vcf")
  
  st_mpileup_bcf_call <- paste("samtools mpileup",st_m,st_F,st_u,st_f,st_in,st_l,st_break,bcf_call,bcf_multi_al,bcf_O,bcf_out)
  cat(st_mpileup_bcf_call,"\n")
  system2("samtools", c("mpileup", st_m, st_F, st_u, st_f, st_in, st_l, st_break, bcf_call, bcf_multi_al, bcf_O, bcf_out))
  cat("\n")
  
  system2("samtools", c("mpileup", st_m, st_f, st_u, st_f, st_in, st_break, bcf_call, "-c", st_break, "vcfutils.pl", "vcf2fq", ">", paste0(results.directory, "Fastq/", sequence, "_3DL3cons.fastq")))
  
  #file.copy(paste0(sequence, "_3DL3.1.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DL3.1.fastq"))
  #file.copy(paste0(sequence, "_3DL3.2.fastq"), paste0(results.directory, "Fastq/", sequence, "_3DL3.2.fastq"))
  
  file.copy(paste0(sequence, "_3DL3nuc.vcf"), paste0(results.directory, "Vcf/", sequence, "_3DL3nuc.vcf"))
  
  file.remove(paste0(sequence, "_3DL3nuc.vcf"))
  file.remove(paste0("r",sequence,"_3DL3in.1.fastq"))
  file.remove(paste0("r",sequence,"_3DL3in.2.fastq"))
  file.remove(paste0(sequence,"_3DL3.1.fastq"))
  file.remove(paste0(sequence,"_3DL3.2.fastq"))
  file.remove(paste0(sequence,"_3DL3.sam"))
  file.remove(paste0(sequence,"_3DL3.bam"))
  file.remove(paste0(sequence,"_3DL3.sorted.bam"))
  file.remove(paste0(sequence,"_KIR3DL3.1.fastq"))
  file.remove(paste0(sequence,"_KIR3DL3.2.fastq"))
  file.remove(paste0(sequence,".temp"))
  file.remove(paste0(sequence,"_notmapped.1.fastq"))
  file.remove(paste0(sequence,"_notmapped.2.fastq"))
}
