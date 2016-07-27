#### Upcoming stuff
## INDEL functionality
## Code parallel option
#### This part is going to be much harder than expected, maybe for last
## Add 3DL3
## Test MIRA replace by implementing 2DP1, 2DS4, 3DL1S1

ping_extractor <- function(sample.location = "../Sequences/", fastq.pattern.1 = "_R1_001.fastq", fastq.pattern.2 = "_R2_001.fastq", bowtie.threads = 18, threshold.kff = 0.2) {
  
  ping.ready <- function() {
    cat("----- Getting PING_grabber ready -----\n\n")
    
    dir.create("PING_sequences", showWarnings = F)
    
    cat("PING_sequences/ directory created.\n\n")
  }
  
  ping.mrG <- function() {
    
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
      bt2_al_conc <- paste0("--al-conc ", "PING_sequences/", sequence, "_KIR_%.fastq")
      bt2_un <- "--un dump.me"
      
      cat(sequence,"\n\n")
      cat("bowtie2", bt2_p, bt2_5, bt2_3, bt2_L, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_stream, bt2_al_conc, bt2_un)
      system2("bowtie2", c(bt2_p, bt2_5, bt2_3, bt2_L, bt2_i, bt2_min_score, bt2_I, bt2_X, bt2_x, bt2_1, bt2_2, bt2_stream, bt2_al_conc, bt2_un))
      cat("\n\n")
      
      file.remove(paste0(sequence, ".temp"))
      file.remove(paste0("dump.me"))
    }
    
   
    sequence_list <- get_sequence_list()
      
    for(i in 1:length(sequence_list)) {
      cat("\n\nRunning MrGrabWaller on: ")
      grabber(sequence_list[i])
    }
    
    cat("MrGrabwaller is complete. Extracted reads are deposited in the Sequences folder.\n")
    cat("fastq.patterns have been adjusted to _KIR_1.fastq and _KIR_2.fastq.\n")
  }
  
  ping.ready()
  
  ping.mrG()
}