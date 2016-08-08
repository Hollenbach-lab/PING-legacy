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

ping_run_all <- function(
  sample.location = "Sequences/",
  fastq.pattern.1 = "_1.fastq",
  fastq.pattern.2 = "_2.fastq",
  bowtie.threads = 4,
  results.directory = ""
  ){
  
  source("PING_allele_caller_v1.0.R")
  source("PING_extractor_v1.0.R")
  source("PING_gc_caller_v1.1.R")
  
  library(data.table)
  
  # Creates results directory
  results_directory <- function() {
    cat("----- Getting PING ready -----\n\n")
    
    if(results.directory != ""){
      save_to <- results.directory
    }else{
      save_to <- paste0("All_results_", format(Sys.time(), "%Y_%m_%d_%H:%M"), "/")
      
      count <- 1
      while(file.exists(save_to)) {
        save_to <- paste0("All_results_", format(Sys.time(), "%Y_%m_%d_%H:%M"), "_", count, "/")
        count <- count + 1
      }
    }
    dir.create(save_to)
    
    cat(paste("Results being saved to", save_to, "\n\n"))
    
    return(save_to)
  }
  
  ping_extractor(sample.location = sample.location, 
                 fastq.pattern.1 = fastq.pattern.1, 
                 fastq.pattern.2 = fastq.pattern.2, 
                 bowtie.threads = bowtie.threads)
  
  sample.location <- "PING_sequences/"
  
  is_gz <- last(unlist(strsplit(fastq.pattern.1, ".", fixed = T))) == "gz"
  
  if(is_gz){
    fastq.pattern.1 <- "_1.fastq.gz"
    fastq.pattern.2 <- "_2.fastq.gz"
  }else{
    fastq.pattern.1 <- "_1.fastq"
    fastq.pattern.2 <- "_2.fastq"
  }

  results <- results_directory()
  
  ping_gc_caller(sample.location = sample.location, results.directory = results, read.cap = 40000)
  
  ping_gc_output <- paste0(results, "Combined_results.csv")
  
  ping_allele_caller(sample.location = sample.location,
                     fastq.pattern.1 = fastq.pattern.1,
                     fastq.pattern.2 = fastq.pattern.2,
                     bowtie.threads = bowtie.threads,
                     ping.gc.output = ping_gc_output,
                     results.directory = results)
  
}
