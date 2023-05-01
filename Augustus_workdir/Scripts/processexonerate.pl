#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

#

my $config = configure(scalar @ARGV);

open (CMDS, ">$config->{'exonerate_runcmds'}") or die $!;
open (CIGAR, ">$config->{'exonerate_cigar'}") or die $!;
open (QGFF, ">$config->{'exonerate_querygff3'}") or die $!;
open (TGGFF, ">$config->{'exonerate_targetgenesgff3'}") or die $!;
open (TOGFF, ">$config->{'exonerate_targetothergff3'}") or die $!;
open (HINT, ">$config->{'exonerate_hints'}") or die $!;
 

print QGFF "##gff-version 3\n";
print TGGFF "##gff-version 3\n";
print TOGFF "##gff-version 3\n";

my %files2delete = ();

my $qalnid = 0;

for (my $chunkid = 1; $chunkid <= $config->{'querychunktotal'}; $chunkid++) {
	my $exonerate_out = $config->{'exonerate_outputbase'} . "_$chunkid". ".exonerate.out";
	
	##house keeping
	die "$exonerate_out.err is not empty. Please check the error before proceeding\n" unless (-z "$exonerate_out.err");
	$files2delete{"$exonerate_out.err"} = "";
	$files2delete{"$exonerate_out.seq"} = "";
	my ($qgff, $tgff);
	my @gffblock = ();
	open (EOUT, "<$exonerate_out") or die $!;
	while (my $eline = <EOUT>) {
		chomp $eline;
		next if ($eline =~ /^#$/ || $eline =~ /^##/ || $eline =~ /^# seqname/ || $eline =~ /^# --- START OF GFF DUMP ---$/ || $eline =~ /^-- completed exonerate analysis$/);
		if ($eline =~ /^Command line:/){
			print CMDS "$eline";
			$qgff = 1;
			$tgff = 0;
		}
		elsif ($eline =~ /^Hostname:/){
			print CMDS ", $eline\n";
		}			
		elsif ($eline =~ /^cigar/){
			print CIGAR "$eline\n";
		}
		elsif ($eline =~ /^# --- END OF GFF DUMP ---$/){
			## process collected gff block
			if ($qgff == 1 && $tgff == 0) {
				$qalnid++;
				my $qgfflines = process_qgffblock(\@gffblock, $qalnid);
				map { print QGFF "$_\n" } @$qgfflines;
			}
			elsif ($tgff == 1 && $qgff == 0) {
				my $genegfflines = process_tgffblock(\@gffblock, $qalnid);
				map { print TGGFF "$_\n" } @$genegfflines;
				map { print TOGFF "$_; alignment_id=$qalnid\n" } @gffblock;
				my $hintgfflines = process_hint_block(\@gffblock, $qalnid);
				map { print HINT "$_\n" } @$hintgfflines;
			}

			##assign gff block to query or target
			if ($qgff == 1 && $tgff == 0) {
				$tgff = 1;
				$qgff = 0;
			}
			elsif ($tgff == 1 && $qgff == 0){
				$qgff = 1;
				$tgff = 0;
			}
			## make current gff block empty
			@gffblock = ();
		}
		else {
			push (@gffblock, $eline);
		}
	}	
	close EOUT;
}

close CMDS;
close CIGAR;
close QGFF;
close TGGFF;
close TOGFF;
close HINT;

map { unlink($_) } keys %files2delete;

sub process_tgffblock {
	my $gfflines = shift;
	my $thisalnid = shift;
	my ($geneid, $transcriptid, $orientation);
	my @modifiedgfflines = ();
	my $exonrank=0;
	foreach my $gffline (@$gfflines){
		my @a = split ("\t",$gffline);
		if ($a[2] eq "gene") {
			if ($a[8] =~ /sequence ((\S+)_i\d+) ; gene_orientation (\S)/){
				$transcriptid = $1;
				$geneid = $2;
				$orientation = $3;
			}
			$a[8] = "ID=$geneid;gene_orientation=$orientation;alignment_id=$thisalnid;from=$transcriptid;";
			push(@modifiedgfflines, join("\t",@a));
			$a[2] = "mRNA";
			$a[8] = "ID=$transcriptid;Parent=$geneid;transcript_orientation=$orientation";
			push(@modifiedgfflines, join("\t",@a));
		}
		elsif ($a[2] eq "exon"){
			$exonrank++;
			$a[8] = "ID=$transcriptid.e$exonrank;Parent=$transcriptid";
			push(@modifiedgfflines, join("\t",@a));			
		}
	}
	return(\@modifiedgfflines);
}

sub process_hint_block {
	my $gfflines = shift;
	my $thisalnid = shift;
	my ($geneid, $transcriptid, $orientation);
	my @modifiedgfflines = ();
	my $exonrank=0;
	my $intronrank=0;
	foreach my $gffline (@$gfflines){
		my @a = split ("\t",$gffline);
		if ($a[2] eq "gene") {
			if ($a[8] =~ /sequence ((\S+)_i\d+) ; gene_orientation (\S)/) {
				$transcriptid = $1;
				$geneid = $2;
				$orientation = $3;
			}
			$a[8] = "ID=$geneid;gene_orientation=$orientation;alignment_id=$thisalnid";
			$a[2] = "mRNA";
			$a[8] = "ID=$transcriptid;Parent=$geneid;transcript_orientation=$orientation";
		}
		elsif ($a[2] eq "exon"){
			$exonrank++;
			$a[8] = "grp=$transcriptid;pri=$config->{'hintpriority'};src=$config->{'hintsource'}";
			push(@modifiedgfflines, join("\t",@a));
		}
		elsif ($a[2] eq "intron"){
			$intronrank++;
			$a[8] = "grp=$transcriptid;pri=$config->{'hintpriority'};src=$config->{'hintsource'}";
			push(@modifiedgfflines, join("\t",@a));
		}
	}
	return(\@modifiedgfflines);
}

sub process_qgffblock {
	my $gfflines = shift;
	my $thisalnid = shift;
	my @modifiedgfflines = ();
	foreach my $gffline (@$gfflines){
		my @a = split ("\t",$gffline);
		my @b = split (";", $a[8]);
		my @c = ();
		foreach my $attribute (@b){
			$attribute =~ s/^ +//;
			$attribute =~ s/ +$//;
			$attribute =~ s/Target /Target=/;
			$attribute =~ s/Align /Align=/;
			$attribute =~ s/alignment_id 0/alignment_id=$thisalnid/;
			push (@c, $attribute);
		}
		$a[8] = join(";", @c);
		push (@modifiedgfflines, join("\t",@a));
	}
	return(\@modifiedgfflines);
}


sub configure {
  my $args = shift;
  my $config = {};
  GetOptions(
  $config,
  'inputfasta|i=s',
  'outputdir|o=s',
  'targetgenome|r=s',
  'querychunktotal|c:i',
  'hintpriority|h:i',
  'hintsource|s=s',
  ) or usage(1);

  unless ($config->{'outputdir'} && -d $config->{'outputdir'}) {
    print "\nERROR: Provide exonerate output directory path.";
    usage(1);
  }
	else {
		$config->{'outputdir'} = abs_path($config->{'outputdir'});
	}

  if ($config->{'inputfasta'} && -e $config->{'inputfasta'}) {
    $config->{'inputfasta'} = abs_path($config->{'inputfasta'});
  }
  else{
    print "\nERROR: Provide path to transcriptome fasta used for exonerate input.";
    usage(1);
  }

  if ($config->{'targetgenome'} && -e $config->{'targetgenome'}) {
    $config->{'targetgenome'} = abs_path($config->{'targetgenome'});
  }
  else{
    print "\nERROR: Provide path to the genome used as target for exonerate.";
    usage(1);
  }

  unless($config->{'querychunktotal'}){
    print "\nERROR: Provide the total number of query chunks defined for exonerate (integer).";
    usage(1);
  }

  unless($config->{'hintpriority'}){
	print "\nERROR: Provide the hint priority level for your hints (integer).";
	usage(1);
  }

  unless($config->{'hintsource'}){
	print "\nERROR: Provide the hint source for your hints (string).";
	usage(1);
  }

	##${outputdir}/$(basename ${inputfasta} .fasta)_{}.exonerate.out
  $config->{'exonerate_outputbase'} = $config->{'outputdir'} . "/" . basename($config->{'inputfasta'}, ".fasta");
  $config->{'exonerate_runcmds'} = $config->{'exonerate_outputbase'}.".exonerate.cmds";
  $config->{'exonerate_cigar'} = $config->{'exonerate_outputbase'}.".exonerate.cigar";
  $config->{'exonerate_querygff3'} = $config->{'exonerate_outputbase'}.".exonerate.query.gff3";
  $config->{'exonerate_targetgenesgff3'} = $config->{'exonerate_outputbase'}.".exonerate.target.genes.gff3";
  $config->{'exonerate_targetothergff3'} = $config->{'exonerate_outputbase'}.".exonerate.target.other.gff3";
  $config->{'exonerate_hints'} = $config->{'exonerate_outputbase'}.".exonerate.target.hints.gff3";
  return $config;
}


sub usage {
  my $exit_code = shift;
	print "\n\n-----\n\n";
  print <<USAGEMSG;
USAGE:
  processexonerate.pl -inputfasta input.fasta -outputdir outputdir -targetgenome genome.fasta -querychunktotal 1
  This utility will take the output of embarassingly parallel tasks of exonerate and merge into meaningful files for downstream uses.
Options:
  -inputfasta/-i       Input transcriptome fasta used for exonerate.
  -outputdir/-o        Path to output directory of exonerate.
  -querychunktotal/-c  Number of chunks used for running embarrasingly parallel exonerate tasks.
  -targetgenome/-r     Reference genome used as target for exonerate.
  -hintpriority/-h     Hint priority level (integer)
  -hintsource/-s       Hint source (E for EST and P for Protein when using default Augustus extrinsic file, using other string is possible as long as they are defined in your custom Augustus extrinsic file)
USAGEMSG
  exit($exit_code);
}