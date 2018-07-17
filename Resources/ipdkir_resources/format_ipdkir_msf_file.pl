#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use DataBrowser;

#=======================================================================#
# PURPOSE:                                                              #
#  - Convert syntax for non nucleotide characters in MSF files          #
#  - Currently the MSF files on IPD-KIR are not formatted correctly     #
#    so we are unable to tell what sequence is unknown                  #
#  - As a result we have to compare the MSF output with the multiple    #
#    alignments generated through the IPD-KIR website to determine this # 
#  - This Script is TEMPORARY, once IPD-KIR has been updated this code  #
#    can be deleted                                                     #
#=======================================================================#

die "
Usage: $0 <KIR locus> <Copied alignment output from IPD-KIR>
 Ex:
   $0 2DL1 KIR2DL1_mult_alignment_CDS.6-6-18.ipdkir
" unless @ARGV ==2;
my ($LOCUS, $IPDKIR_ALN) = @ARGV;

# Download most updated alleles (CDS only) for this particular KIR locus from IPD-KIR (MSF)
my $msf_file = "./KIR$LOCUS\_nuc.msf";
if (! -e $msf_file) {
	system("wget ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/msf/KIR$LOCUS\_nuc.msf") == 0 or die;
}

# Store postion info from IPD-KIR multiple alignment
my %allele_hash;

my $allele = "";
my $seq    = "";
open(my $IN1, "<", $IPDKIR_ALN) or die "Error: cannot read file\n";
while (my $line = <$IN1>) {
	chomp($line);
	next if $line !~ /\w+/;
	next if $line =~ /cDNA/;
	
	my @info = split(/\s+/, $line);
	if ($line =~ /$LOCUS/) {
		$allele = $info[1];
		$allele = "KIR" . $allele;
		
		$allele_hash{$allele} = "" if !defined $allele_hash{$allele};
		
		for (my $i = 2; $i < @info; $i++) {
			$allele_hash{$allele} .= $info[$i];
		}
	}
	else {
		next;
	}
	
}
close $IN1;

#browse(\%allele_hash);

# Format sequence strings and Extract Reference Allele
my $ref_allele;
my $seq_ref_allele;
my $max_allele_name_len = 0;
foreach my $allele (keys %allele_hash) {
	my $seq = $allele_hash{$allele};
	$seq =~ s/\|//g;
	$allele_hash{$allele} = $seq;
	
	# Get length of allele with the longest name (used for spacing in the end)
	my $name_length = length($allele);
	$max_allele_name_len = $name_length if $name_length > $max_allele_name_len;
	
	# Establish which allele is the reference
	if ($seq !~ /\-/) {
		$ref_allele     = $allele;
		$seq_ref_allele = $seq;
	}
}


# Correct the MSF file output it 
my %current_positions;

open(my $IN2, "<", $msf_file) or die "Error: cannot read MSF file for $LOCUS\n";
while (my $line = <$IN2>) {
	chomp($line);
	if ($line =~ /$LOCUS/ and $line !~ /Name:/) {
		my @info = split(/\s+/, $line);
		my $start = 1;
		my $allele = $info[0];
		if (length($allele) == 0) {
			$start = 2;
			$allele = $info[1];
		}
		for (my $i = $start; $i < @info; $i++) {
			my $bp_block      = $info[$i];
			my $bp_block_size = length($bp_block);
			
			$current_positions{$allele} = 0 if !defined $current_positions{$allele};
			my $pos = $current_positions{$allele};
			
			my $known_seq = $allele_hash{$allele};
			my $known_seq_block = substr($known_seq,$pos,$bp_block_size);
			
			# If there is unknown sequence in the MSF file, make sure they are denoted with '*'s
			my $formatted_bp_block = format_bp_block($bp_block,$known_seq_block);
# 			print "$allele\t$pos\t$bp_block\t$known_seq_block\t$formatted_bp_block\n";
			
			$info[$i] = $formatted_bp_block if $formatted_bp_block =~ /\*/;
			
			$current_positions{$allele} += $bp_block_size
		}
# 		print "\n";
		
		### Output formatted MSF line
		my $formatted_line = join(' ', @info);
		# Space formatting to make the output look pretty
		my $num_spaces = $max_allele_name_len - length($allele);
		my $space = '';
		if ($num_spaces > 0) {
			$num_spaces -= 1 if $start == 2;  # Already 1 space at the beginning
			for (my $s=0; $s < $num_spaces;$s++) {$space .= ' ';}
		}
		print "$space$formatted_line\n";
	}
	else {
		print "$line\n";
	}
}
close $IN2;

#===================#
# SUBROUTINES       #
#===================#

# Convert BP block in MSF file if it is supposed to have astrix chars
sub format_bp_block{
	my ($current_block, $known_block) = @_;
	
	for (my $i = 0; $i < length($known_block); $i++) {
		my $bp = substr($known_block,$i,1);
		if ($bp eq '*') {substr($current_block,$i,1) = $bp;}
		else            {next;}
	}
	
	return $current_block;
}




