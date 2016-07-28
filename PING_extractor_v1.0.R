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

ping_extractor <- function(
  sample.location = "../Sequences/", 
  fastq.pattern.1 = "_1.fastq.gz", 
  fastq.pattern.2 = "_2.fastq.gz", 
  bowtie.threads = 18
  ) {
  
  ping.ready <- function() {
    cat("----- Getting PING_grabber ready -----\n\n")
    
    dir.create("PING_sequences", showWarnings = F)
    
    cat("PING_sequences/ directory created.\n\n")
  }
  
  get_sequence_list <- function(folder.name = sample.location, file.pattern = fastq.pattern.1) {
    
    sequence_list = list.files(file.path(folder.name), pattern = file.pattern)
    
    if (is.na(sequence_list[1])) {
      stop("No sequences found, please place fastq files in the Sequences folder.")
    } else {
      sequence_list <- gsub(file.pattern, "", sequence_list)
      cat(paste("Found sequences: ", paste(sequence_list, collapse = "\n"), sep = "\n"))
      cat("\n")
      return(sequence_list)
    }
  }
  
  ping.mrG <- function(sequence.list) {
    
    # Pull out reads that match any KIR ---------------------------------------
    
    grabber <- function(sequence) {
      
      bt2_p <- paste0("-p", bowtie.threads)
      bt2_5 <- "--trim5 3"
      bt2_3 <- "--trim3 7"
      bt2_L <- "-L 20"
      bt2_i <- "-i S,1,0.5"
      bt2_min_score <- "--score-min L,0,-0.187"
      bt2_I <- "-I 75"
      bt2_X <- "-X 1000"
      bt2_x <- "-x Resources/grabber_resources/Filters/mrG/output"
      
      bt2_1 <- paste0("-1 ", sample.location, sequence, fastq.pattern.1)
      bt2_2 <- paste0("-2 ", sample.location, sequence, fastq.pattern.2)
      bt2_stream <- paste0("-S ", sequence, ".temp")
      
      if(is_gz){
        bt2_al_conc <- paste0("--al-conc-gz ", "PING_sequences/", sequence, "_KIR_%.fastq.gz")
      }else{
        bt2_al_conc <- paste0("--al-conc ", "PING_sequences/", sequence, "_KIR_%.fastq")
      }
      
      bt2_un <- "--un dump.me"
      
      cat(sequence,"\n\n")
      cat("bowtie2", bt2_p, bt2_5, bt2_3, bt2_L, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_stream, bt2_al_conc, bt2_un)
      system2("bowtie2", c(bt2_p, bt2_5, bt2_3, bt2_L, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_stream, bt2_al_conc, bt2_un))
      cat("\n\n")
      
      file.remove(paste0(sequence, ".temp"))
      file.remove(paste0("dump.me"))
    }
      
    for(sample in sequence.list) {
      cat("\n\nRunning MrGrabWaller on: ")
      grabber(sample)
    }
    
    cat("MrGrabwaller is complete. Extracted reads are deposited in the Sequences folder.\n")
    cat("fastq.patterns have been adjusted to _KIR_1.fastq(.gz) and _KIR_2.fastq(.gz).\n")
  }
  
  ping.ready()
  
  sequence_list <- get_sequence_list()
  
  is_gz <- last(unlist(strsplit(fastq.pattern.1, ".", fixed = T))) == "gz"
  
  ping.mrG(sequence_list)
}