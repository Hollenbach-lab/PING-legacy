IPD KIR Resources Used
======================

## Directory: IPD_KIR_alignments
* Contains multiple alignments for each KIR locus
* Each multiple alignment output was collected directly from IPD-KIR 

## Directory: IPD_KIR_MSF_Files
* Contains MSF files for each KIR locus
* Each file was downloaded from the IPD-KIR FTP site

## Formatting MSF files: format_ipdkir_msf_file.pl
* DISCLAIMER: This formatting step is a temporary step until the FTP resources are updated
* This perl script will format all the unknown positions in each MSF file
* Currently, IPD-KIR MSF files denote unknown and deletion positions with '*', this perl script changes that based on IPD-KIR multiple alignements 
