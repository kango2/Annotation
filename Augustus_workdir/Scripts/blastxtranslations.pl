#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my %genetic_code  = (
'TCA' => 'S','TCC' => 'S','AGC' => 'S','AGT' => 'S','TCG' => 'S','TCT' => 'S',
'TTC' => 'F','TTT' => 'F',
'TTA' => 'L','TTG' => 'L','CTC' => 'L','CTG' => 'L','CTT' => 'L','CTA' => 'L',
'TAC' => 'Y','TAT' => 'Y',
'TAA' => '*','TAG' => '*','TGA' => '*',
'TGC' => 'C','TGT' => 'C',
'TGG' => 'W',
'CCA' => 'P','CCC' => 'P','CCG' => 'P','CCT' => 'P',
'CAT' => 'H','CAC' => 'H',
'CAA' => 'Q','CAG' => 'Q',
'CGA' => 'R','CGC' => 'R','AGA' => 'R','AGG' => 'R','CGG' => 'R','CGT' => 'R',
'ATA' => 'I','ATC' => 'I','ATT' => 'I',
'ATG' => 'M',
'ACA' => 'T','ACC' => 'T','ACG' => 'T','ACT' => 'T',
'AAC' => 'N','AAT' => 'N',
'AAA' => 'K','AAG' => 'K',
'GTA' => 'V','GTC' => 'V','GTG' => 'V','GTT' => 'V',
'GCA' => 'A','GCC' => 'A','GCG' => 'A','GCT' => 'A',
'GAC' => 'D','GAT' => 'D',
'GAA' => 'E','GAG' => 'E',
'GGA' => 'G','GGC' => 'G','GGG' => 'G','GGT' => 'G',
);

my ($queryfasta, $reffasta, $blastx, $prefix) = @ARGV;

my $primary_id = "";
my $description = "";
my %cdnaseq = ();
open (F, "<$queryfasta") or die $!;
while (<F>){
  chomp $_;
  if ($_ =~ />(\S+)\s+?(.*)/){
    $primary_id = $1;
    $cdnaseq{$primary_id}{'description'} = $2;
  }
  else{
    $cdnaseq{$primary_id}{'seq'} .= $_;
  }
}
close F;

my %peptideseq = ();
open (F, "<$reffasta") or die $!;
while (<F>){
  chomp $_;
  if ($_ =~ />\S+?\|(\S+?)\|.*\s+?(.*)/){
    $primary_id = $1;
		
    $peptideseq{$primary_id}{'description'} = $2;
  }
  else{
    $peptideseq{$primary_id}{'seq'} .= $_;
  }
}
close F;

#my %peptideseq = ();
#open (F, "<$reffasta") or die $!;
#while (<F>){
#  chomp $_;
#  if ($_ =~ />(\S+)\s+?(.*)/){
#    $primary_id = $1;
#    $peptideseq{$primary_id}{'description'} = $2;
#  }
#  else{
#    $peptideseq{$primary_id}{'seq'} .= $_;
#  }
#}
#close F;


###/home/hardip/reference/uniprot_sprot.fasta

open (CDS,   ">$prefix.cds.all.fa") or die $!;
open (CDNA,  ">$prefix.cdna.all.fa") or die $!;
open (AASEQ, ">$prefix.pep.all.fa") or die $!;


my %translationseen = ();
open (F, "<$blastx") or die $!;
while (<F>){
  chomp $_;
  my @a = split ("\t", $_);
  my @b = split (/\|/, $a[1]);
  $a[11] =~ s/\s*//g;
  #next if ($a[3] <= 30 || $a[10] > 0.001);
  next if (exists $translationseen{$a[0]});
  my ($cds, $translation, $cds_start, $cds_end, $cdna, $strand, $startcodon, $stopcodon, $rellen) = get_translation(@a);
  next if (length($translation) < 50);
  print CDS   ">$a[0] $cdnaseq{$a[0]}{'description'} $cds_start:$cds_end:$strand:$startcodon:$stopcodon:$rellen\n$cds\n";
  print CDNA  ">$a[0] $cdnaseq{$a[0]}{'description'}\n$cdna\n";
  print AASEQ ">$a[0] $cdnaseq{$a[0]}{'description'}\n$translation\n";
  $translationseen{$a[0]}="";
}
close F;

foreach my $seq (keys %cdnaseq){
  unless (exists $translationseen{$seq}){
    my ($cds, $translation, $cds_start, $cds_end, $cdna, $strand, $startcodon, $stopcodon, $rellen) = longestORF($cdnaseq{$seq}{'seq'});
    if (length($translation) >= 50){
      print CDS   ">$seq $cdnaseq{$seq}{'description'} $cds_start:$cds_end:$strand:$startcodon:$stopcodon:$rellen\n$cds\n";
      print CDNA  ">$seq $cdnaseq{$seq}{'description'}\n$cdna\n";
      print AASEQ ">$seq $cdnaseq{$seq}{'description'}\n$translation\n";
    }
    else{
      print CDNA ">$seq $cdnaseq{$seq}{'description'}\n$cdnaseq{$seq}{'seq'}\n";
    }
  }
}

exit;


sub get_translation {
  my ($qid, $sid, $pid, $alen, $mm, $gapo, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = @_;
  my $frame = ($qstart < $qend) ? 1 : -1;
  
  ###check if you can look upstream/downstream or not
  my $cdna            = $cdnaseq{$qid}{'seq'};
  my $refseq          = $peptideseq{$sid}{'seq'};
  my $aln_cds_start   = $qstart;
  my $aln_cds_end     = $qend;
  my $has_stop_codon  = "no";
  my $has_start_codon = "no";
  ###reverse the sequence for -ve frame
  if ($frame < 0){
    $cdna = reverse($cdna);
    $cdna =~ tr/ACGT/TGCA/;
    $aln_cds_start = length($cdna) - $qstart + 1;
    $aln_cds_end   = length($cdna) - $qend   + 1;
  }
  ####get the translation as defined by the alignments
  my $aln_translation = "";
  my $aln_cds         = substr($cdna, $aln_cds_start - 1, $aln_cds_end - $aln_cds_start + 1);
  
  for (my $i=0;$i<=length($aln_cds)-3;$i+=3){
    $aln_translation .= (exists $genetic_code{substr($aln_cds, $i, 3)}) ? $genetic_code{substr($aln_cds, $i, 3)} : "X";
  }
  
  ###rules for extension are that it needs to go up to the stop codon and then search for start codon within it
  my $five_extra_translation = undef;
  my $newstart = undef;
  my %five_translations = ();
  for (my $i=$aln_cds_start - 4;$i>=0;$i-=3){
    my $this_aa = (exists $genetic_code{substr($cdna, $i, 3)}) ? $genetic_code{substr($cdna, $i, 3)} : "X";
    if ($this_aa eq "M"){
      $five_extra_translation = (defined $five_extra_translation) ? "$this_aa$five_extra_translation" : "$this_aa";
      $newstart = $i + 1;
      $five_translations{'M'}{$newstart} = $five_extra_translation;
    }
    elsif ($this_aa eq "*"){
      $newstart = $i + 4;
      $five_translations{'other'}{$newstart} = (defined $five_extra_translation) ? "$five_extra_translation" : "";
      goto STOPFOUND;
    }
    else {
      $five_extra_translation = (defined $five_extra_translation) ? "$this_aa$five_extra_translation" : "$this_aa";
      $newstart = $i+1;
    }
  }
  $five_translations{'other'}{$newstart} = $five_extra_translation if (defined $five_extra_translation && defined $newstart);
STOPFOUND:
  
  ###if the extra translation has found the start codon then well and good, if not then see if there is a start codon in the aln_translation
  if (exists $five_translations{'M'}){
    foreach my $mposition (sort {$a<=>$b} keys %{$five_translations{'M'}}){
      $aln_cds_start = $mposition;
      $aln_cds = substr($cdna, $aln_cds_start - 1, $aln_cds_end - $aln_cds_start + 1);
      $aln_translation = "$five_translations{'M'}{$mposition}$aln_translation";
      $has_start_codon = "yes";
      last;
    }
  }
  else{
    my $methionine_position = index($aln_translation, "M");
    ###start codon found within the alignment range
    if ($methionine_position >= 0 && $sstart <= length($peptideseq{$sid}) / 10 && $methionine_position <= length($aln_translation) / 10){
      $aln_cds_start = $aln_cds_start + ($methionine_position * 3);
      $aln_cds = substr($cdna, $aln_cds_start - 1, $aln_cds_end - $aln_cds_start + 1);
      $aln_translation = substr($aln_translation, $methionine_position);
      $has_start_codon = "yes";
    }
    ###no start codon found
    else{
      if (exists $five_translations{'other'}){
        foreach my $mposition (sort {$a<=>$b} keys %{$five_translations{'other'}}){
          $aln_cds_start = $mposition;
          $aln_cds = substr($cdna, $aln_cds_start - 1, $aln_cds_end - $aln_cds_start + 1);
          $aln_translation = "$five_translations{'other'}{$mposition}"."$aln_translation";
          last;
        }
      }
    }
  }
  #####check if the stop codon is already present in the translation, if so then adjust the coordinates and sequences
  my $stopcodon_position = index($aln_translation, "*");
  if ($stopcodon_position >= 0){
    $has_stop_codon = "yes";
    $aln_cds_end = $aln_cds_start + ($stopcodon_position * 3 - 1) + 3; # + 3; ###added three to include the stop codon
    $aln_cds = substr($cdna, $aln_cds_start - 1, $aln_cds_end - $aln_cds_start + 1);
    $aln_translation = substr($aln_translation, 0, $stopcodon_position+1); #+1); ##added 1 to length for including the stop codon
  }
  ###no stop codon found so search for it downstream of the existing sequence
  else{
    my $three_extra_translation = undef;
    my $newend = undef;
    for (my $i=$aln_cds_end; $i<=length($cdna) - 3; $i+=3){
      my $this_aa = (exists $genetic_code{substr($cdna, $i, 3)}) ? $genetic_code{substr($cdna, $i, 3)} : "X";
      if ($this_aa eq "*"){
        ###mask the following after testing
        $has_stop_codon = "yes";
        $three_extra_translation .= $this_aa;
        $newend = $i + 3; # + 3;
        last;
      }
      else{
        $three_extra_translation .= $this_aa;
        $newend = $i + 3;
      }
    }
    if (defined $three_extra_translation && defined $newend){
      $aln_cds_end = $newend;
      $aln_cds = substr($cdna, $aln_cds_start - 1, $aln_cds_end - $aln_cds_start + 1);
      $aln_translation .= $three_extra_translation;
    }
  }
  return ($aln_cds, $aln_translation, $aln_cds_start, $aln_cds_end, $cdna, $frame, $has_start_codon, $has_stop_codon, sprintf("%.2f", length($aln_translation)/length($refseq)));
}

sub longestORF {
  # longorf.pl v0208020920
  # (c) Dan Kortschak 2002
  my $best=0;
  my ($bests,$beste,$beststrand)=(-1,-1,0);
  my $bestorf="";

  my $dna= $_[0];
  my $revcom = reverse($dna);
  $revcom =~ tr/ACGT/TGCA/;

  my %strand=(
  '+'=>$dna,
  '-'=>$revcom
  );

  foreach my $direction (keys %strand) {
    my @starts=();
    my @ends=();
    for (my $frame=0;$frame<3;$frame++) {
      unless ($strand{$direction}=~m/^.{$frame}(taa|tga|tag)/i) {
        push @starts,$frame+1;
      }
    }
    while ($strand{$direction}=~m/(atg)/gi) {
      push @starts,pos($strand{$direction})-2;
    }

    while ($strand{$direction}=~m/(taa|tga|tag)/gi) {
      #push @ends,pos($strand{$direction})-2;
      push @ends,pos($strand{$direction})+1;
    }
    push @ends,(length($dna)-2,length($dna)-1,length($dna));

    for my $s (@starts) {
      for my $e (@ends) {
        if ($e%3==$s%3 and $e>$s) {
          if ($e-$s>$best) {
            $best=$e-$s;
            ($bests,$beste,$beststrand)=($s,$e,$direction);
            $bestorf=substr($strand{$direction}, $s-1, $e-$s);
          }
          last
        } else {
          next
        }
      }
    }
  }
  my $translation = "";
  for (my $i=0;$i<length($bestorf);$i+=3){
    $translation .= $genetic_code{substr($bestorf, $i, 3)};
  }
  my $cds_start         = $bests;
  my $cds_end           = $beste;
  my $cdna              = $dna;
  my $translation_frame = 1;
  if ($beststrand eq "-"){
    $cds_start = length($dna) - $beste + 1;
    $cds_end   = length($dna) - $bests + 1;
    $cdna      = $revcom;
    $translation_frame = -1;
  }
  return ($bestorf, $translation, $cds_start, $cds_end, $cdna, $translation_frame, (($translation =~ /^M/) ? "yes" : "no"), (($translation =~ /\*$/) ? "yes" : "no"), "NA");
}
