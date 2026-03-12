#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;



my ($suffix);
my ($maxRank);
my ($minIden, $minPos);
my ($minSize, $maxSize);
my ($UnShift, $UnStop);
my (@RefPep);
my ($CDSRatioR2Q, $CDSRatioQ2R);
my ($Verbose);
my ($help);


GetOptions(
	"suffix|sfx|s:i"	=>	\$suffix,
	"maxRank|maxR|xR:i"	=>	\$maxRank,
	"minIden|minI|mI:f"	=>	\$minIden,
	"minPos|minP|mP:f"	=>  \$minPos,
	"minSize|minS|mS:i"	=>	\$minSize,
	"maxSize|maxS|xS:i"	=> 	\$maxSize,
	"UnShift|usft"	=>	\$UnShift,
	"UnStop|ustp"	=>	\$UnStop,
	"RefPep|pep|rp:s"	=>	\@RefPep,
	"CDSRatioR2Q|R2Q:i"	=>	\$CDSRatioR2Q,
	"CDSRatioQ2R|Q2R:i"	=>	\$CDSRatioQ2R,
	"Verbose|V"	=>	\$Verbose,
	"help|h"	=>	\$help
);






## Default values

$suffix ||= "M";
$maxRank ||= 0 ;
$minIden ||= 0;
$minPos ||= 0;
$minSize ||= 0;
$maxSize ||= 1000000;
$UnShift ||= 0;
$UnStop ||= "0";
$CDSRatioR2Q ||= 0;
$CDSRatioQ2R ||= 0;
$Verbose ||= 0;

## Help message

($help)  && (help() && exit 0);
(!@ARGV && -t STDIN) && (help() && exit 0);

sub help {
	print <<HELP;
Usage: perl $0 <miniprot ouput gff file> [options]


Options:
	-h, --help		print help message
	-s, --suffix		suffix of ID
	-xR, --maxRank		maximum Rank
	-mI, --minIden		minimum Identity value
	-mP, --minPos		minimum Positive value
	-mS, --minSize		minimum size(length) of alignemnt
	-xS, --maxSize		maximum size(length) of alignemnt
	-usft, --UnShift	No frameshift events in alignment
	-ustp, --UnStop		No inframe stop codon
	-pep, --RefPep		Reference protein (FASTA)
	-R2Q, --CDSRatioR2Q	CDS length ratio (Reference/Query)
	-Q2R, --CDSRatioQ2R	CDS length ratio (Query/Reference)
	-v, --Verbose

Examples:
	# input from file
	perl $0 miniprot_output.gff [options]

	perl $0 -v miniprot_output.gff.gz [options]

	# input from STDIN
	miniprot genome.fa protein.fa | perl $0 [options]
	cat miniprot_output.gff | perl $0 [options]

	# Redirection
	perl $0 [options] < miniprot_output.gff

Author: shidishen
HELP
}

if (!@ARGV && -t STDIN) {
    print_help();
    exit 0;
}




## Read reference protein and get length

sub FastaLen {

	my %Len;
	my @peps = shift;
	foreach my $fa(@peps){
		if($fa =~ /\.gz$/){
			open IN,"gunzip -c $fa | ";
		}else{
			open IN,$fa;
		}

		$/=">";<IN>;$/="\n";
		while(<IN>){
			chomp(my $id = (split /\s+/,$_)[0]);
			$/=">";
			chomp(my $seq = <IN>);
			$/="\n";

			$seq =~ s/>//;
			$seq =~ s/\s+//g;
			$Len{$id} = length($seq);
		}


	}

	return %Len;

}

my (%RefLen, %RefPepInfo);

if(@RefPep){
	%RefLen = FastaLen(@RefPep);
	%RefPepInfo = (
		"RefLen"	=>	\%RefLen,
		"CDSRatioR2Q"	=>	$CDSRatioR2Q,
		"CDSRatioQ2R"	=>	$CDSRatioQ2R
	);
}



## Define dispatch table

my %Fvalue = (
	"maxRank"	=>	$maxRank,
	"minIden"	=>	$minIden,
	"minPos"	=>	$minPos,
	"minSize"	=>	$minSize,
	"maxSize"	=>	$maxSize,
	"UnShift"	=>	$UnShift,
	"UnStop"	=>	$UnStop
);

@RefPep && ($Fvalue{"RefPep"} = \%RefPepInfo);

my %Filter = (
	"maxRank"	=>	\&filter_by_Rank,
	"minIden"	=>	\&filter_by_Identity,
	"minPos"	=>	\&filter_by_Positive,
	"minSize"	=>	\&filter_by_minSize,
	"maxSize"	=>	\&filter_by_maxSize,
	"UnShift"	=>	\&filter_by_UnShift,
	"UnStop"	=>	\&filter_by_UnStop,
	"RefPep"	=>	\&filter_by_RefPep
);




## Define filters

sub filter_by_Rank{
	my $block = shift;
	my $minRank = shift;

	my @info = split(/\n/,$block);
	my $Rank = $1 if $info[0] =~ /Rank=(\d+)/;
	($Rank > $minRank) ? return 1 : return 0 ;
}

sub filter_by_Identity{
	my $block = shift;
	my $minIden = shift;

	my @info = split(/\n/,$block);
	my $Identity = $1 if $info[0] =~ /Identity=([^\;]*);/;
	($Identity > $minIden) ? return 1 : return 0 ;
}

sub filter_by_Positive{
	my $block = shift;
	my $minPos = shift;

	my @info = split(/\n/,$block);
	my $Positive = $1 if $info[0] =~ /Positive=([^\;]*);/;
	($Positive > $minPos) ? return 1 : return 0 ;
}


sub filter_by_minSize{
	my $block = shift;
	my $minSize = shift;

	my @info = split(/\n/,$block);

	my @line1 = split(/\t/,$info[0]);
	my $TranscriptLength = abs($line1[4] - $line1[3] + 1);

	($TranscriptLength > $minSize) ? return 1 : return 0 ;

}

sub filter_by_maxSize{
	my $block = shift;
	my $minSize = shift;

	my @info = split(/\n/,$block);

	my @line1 = split(/\t/,$info[0]);

	my $TranscriptLength = abs($line1[4] - $line1[3] + 1);

	($TranscriptLength < $maxSize) ? return 1 : return 0 ;

}

sub filter_by_UnStop{
	my $block = shift;
	my $UnStop = shift;
	my @info = split(/\n/,$block);

	my $Stop = $1 if $info[0] =~ /StopCodon=(\d)/;
	$Stop ||= 0;

	

	($Stop && $UnStop) ? return 0 : return 1 ;
}


sub filter_by_UnShift{
	my $block = shift;
	my $UnShift = shift;

	my @info = split(/\n/,$block);

	my $Shift = $1 if $info[0] =~ /Frameshift=(\d)/;
	$Shift ||= 0;

	($Shift && $UnShift) ? return 0 : return 1 ;
}


sub filter_by_RefPep{

	my $block = shift;
	my $RefPepInfo = shift;
	
	my $RefLen = $RefPepInfo{"RefLen"};
	my $CDSRatioR2Q = $RefPepInfo{"CDSRatioR2Q"};
	my $CDSRatioQ2R = $RefPepInfo{"CDSRatioQ2R"};


	my @info = split(/\n/,$block);
	my $RefID = $1 if $info[0] =~ /Target=(\S+)/;

	my $QryCDSLen = 0;
	foreach my $cdsline (@info[1..$#info]){
		my @line = split(/\t/, $cdsline);
		my $cdslen = abs($line[4] - $line[3] + 1);
		$QryCDSLen += $cdslen;
	}

	my $RefCDSLen = $RefLen{$RefID}*3;

	( ($RefCDSLen/$QryCDSLen > $CDSRatioR2Q) && ($QryCDSLen/$RefCDSLen > $CDSRatioQ2R) ) ? return 1 : return 0;


}






## Define miniprot output formatter

sub miniprotFormatter{

	my($block, $suffix, $Verbose) = @_;
	my $RefID = $1 if $block =~ /Target=(\S+)/;
	my $Rank = $1 if $block =~ /Rank=(\d+)/;
	
	my $ID = "$RefID.$suffix$Rank";

	($block =~ /(ID|Parent)=([^;]*)/) && $block =~ s/$2/$ID/g;

	if($Verbose == 0){
		my @lines = split(/\n/, $block);
		@lines = map {$_ = (split(/;/, $_))[0]} @lines;
		$block = join("\n", @lines)."\n";
	}

	return $block;

}




## While loopping gff file


if(@ARGV){
	($ARGV[0]=~/\.gz$/) ? (open IN,"<::gzip",$ARGV[0]) : (open IN,$ARGV[0]) ;
}else{
	open IN, "<&", \*STDIN;
}



$/="##";<IN>;$/="\n";
while(<IN>){


	my $paf = $_;
	$/="##";
	my $block = <IN>;$paf ||= "";
	$block =~ s/##//;
	$/="\n";
	$block || next;



	## For looping every block to filters
	my @keeps;
	foreach my $V(keys %Fvalue){
		my $keep = $Filter{$V}($block, $Fvalue{$V});
		push(@keeps, $keep);
	}

	## Filter out unfit block
	(grep { $_ == 0 } @keeps) && next;


	## Formatting

	$block = miniprotFormatter($block, $suffix, $Verbose);
	#print("##".$paf);
	print($block);

}
