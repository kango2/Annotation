#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

#
#perl script to append 9th column of stitched gene
#change it to ${workingdir}/Augustus/annotation_stitch/temp
my $input = $ARGV[0];
my $geneid = $ARGV[1];

if (!defined $input || !defined $geneid) {
    die "Usage: perl $0 \${workingdir}/Augustus/annotation_stitch/temp \${CurrentGeneId}\n";
}

my ($transcriptid, $endcoord);
my $CDSrank=0;
my $intronrank=0;
my $TSSrank=0;
my $TTSrank=0;
my $STARTrank=0;
my $STOPrank=0;
my $UTR5rank=0;
my $UTR3rank=0;
open INPUT, $input or die "$!: $input";
while(my $line = <INPUT>){
    chomp $line;
    my @modifiedgfflines = ();
    my @a = split("\t",$line);
    if ($a[2] eq "gene") {
        if ($a[8] =~ /ID=(\S+);/) {
            $transcriptid = "$geneid.t1";
        }
        $a[8] = "ID=$geneid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "mRNA") {
        $a[8] = "ID=$transcriptid;Parent=$geneid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "CDS") {
        $CDSrank++;
        $a[8] = "ID=$transcriptid.CDS$CDSrank;Parent=$transcriptid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "intron") {
        $intronrank++;
        $a[8] = "ID=$transcriptid.intron$intronrank;Parent=$transcriptid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "transcription_start_site") {
        $TSSrank++;
        $a[8] = "ID=$transcriptid.tss$TSSrank;Parent=$transcriptid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "transcription_end_site") {
        $TTSrank++;
        $a[8] = "ID=$transcriptid.tts$TTSrank;Parent=$transcriptid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "five_prime_utr") {
        $UTR5rank++;
        $a[8] = "ID=$transcriptid.5UTR$UTR5rank;Parent=$transcriptid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "three_prime_utr") {
        $UTR3rank++;
        $a[8] = "ID=$transcriptid.3UTR$UTR3rank;Parent=$transcriptid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "start_codon") {
        $STARTrank++;
        $a[8] = "ID=$transcriptid.start$STARTrank;Parent=$transcriptid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    elsif ($a[2] eq "stop_codon") {
        $STOPrank++;
        $a[8] = "ID=$transcriptid.stop$STOPrank;Parent=$transcriptid;";
        push(@modifiedgfflines, join("\t",@a));
    }
    print "@modifiedgfflines\n";
}
close INPUT;
