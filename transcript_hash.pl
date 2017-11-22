#!/usr/bin/perl

use strict;
#use warnings;
use File::Basename;

die "perl $0 <input_file> <input_file2> <output_file1>" unless @ARGV==3;

my $in_1=shift;
my $in_2=shift;
my $out=shift;
my %hash=();
my %hash_ID=();
my %Chr = (
    'NC_000001.10'      =>      'chr1',
    'NC_000002.11'      =>      'chr2',
    'NC_000003.11'      =>      'chr3',
    'NC_000004.11'      =>      'chr4',
    'NC_000005.9'       =>      'chr5',
    'NC_000006.11'      =>      'chr6',
    'NC_000007.13'      =>      'chr7',
    'NC_000008.10'      =>      'chr8',
    'NC_000009.11'      =>      'chr9',
    'NC_000010.10'      =>      'chr10',
    'NC_000011.9'       =>      'chr11',
    'NC_000012.11'      =>      'chr12',
    'NC_000013.10'      =>      'chr13',
    'NC_000014.8'       =>      'chr14',
    'NC_000015.9'       =>      'chr15',
    'NC_000016.9'       =>      'chr16',
    'NC_000017.10'      =>      'chr17',
    'NC_000018.9'       =>      'chr18',
    'NC_000019.9'       =>      'chr19',
    'NC_000020.10'      =>      'chr20',
    'NC_000021.8'       =>      'chr21',
    'NC_000022.10'      =>      'chr22',
    'NC_000023.10'      =>      'chrX',
    'NC_000024.9'       =>      'chrY',
    'NC_012920.1'       =>      'chrM',
   );

open(OUT,">$out");
################get transript ID
my @file=split/\./,$in_1;
if(($file[-1] eq "gtf") || ($file[-1] eq "xls")){
    open(IN,"$in_1") || die "Can't open file1 $in_1\n";
    while(my $line=<IN>){
        chomp($line);
        my @line1=split/\t/,$line;
        my @line2=split/NM_/,$line1[2];
        for(my $i=1;$i<@line2;$i++){
            my @line3=split/:/,$line2[$i];
            my $Genbank="Genbank:NM_".$line3[0];
            my $ID=$Genbank."\t".$line3[2];
            $hash{$ID}=0;
            $hash_ID{$Genbank}=1;
        }
    }
}elsif($file[-1] eq "txt"){
    open(IN,"$in_1") || die "Can't open file1 $in_1\n";
    while(my $line=<IN>){
        chomp($line);
        my @line1=split/\t/,$line;
        my $Genbank="Genbank:".$line1[0];
        my $ID=$Genbank."\t".$line1[1];
        $hash{$ID}=0;
        $hash_ID{$Genbank}=1;
    }
}elsif(($in_1=~/NM_/) && ($in_1=~/c\./)){
    my @line=split/:/,$in_1;
    my $Genbank="Genbank:".$line[0];
    my $ID=$Genbank."\t".$line[1];
    $hash{$ID}=0;
    $hash_ID{$Genbank}=1;
}
my $y;
my %hash_NM=();
my %hash_Pa=();
my $Parent;
open(IN2,"$in_2") || die "Can't open file2 $in_2\n";
while(my $Line=<IN2>){
    my @trans=split/,/,$Line;
    my @Trans=split/\./,$trans[1];
    my @Line=split/\t/,$Line;
    if(($hash_ID{$Trans[0]}) && ($Line[2] eq "exon")){
        $Line=~/\;Parent=rna(\d+)\;/;
        $Parent="Parent=rna".$1;
        $hash_NM{$Trans[0]}=$Parent;
    }
    if($Line=~/\;$Parent\;/){
        $hash_Pa{$Parent}.=$Line;
    }
}
my $j;
my @key=keys %hash;
foreach $j(@key){
    my $x=0;
    my @pos;
    my @J=split/\t/,$j;
    my %hash_length=();
    my %hash_chr=();
    $J[1]=~/(\d+)/;
    my $num=$1;
    if($hash_NM{$J[0]}){
        my @values=split/\n/,$hash_Pa{$hash_NM{$J[0]}};
        my $length;
        for(my $i=0;$i<@values;$i++){
            my @CDS=split/\t/,$values[$i];
            if(($CDS[2] eq "CDS") && ($CDS[8]=~/\;$hash_NM{$J[0]}\;/)){
                $length+=$CDS[4]-$CDS[3]+1;
                $pos[$x]=$length;
                $x++;
                my $length_strand=$length."\t".$CDS[6];
                if($CDS[6] eq "+"){
                    $hash_length{$length_strand}=$CDS[3];
                }else{
                    $hash_length{$length_strand}=$CDS[4];
                }
                $hash_chr{$J[0]}=$Chr{$CDS[0]};
            } 
        }
    }
    for(my $q=0;$q<@pos;$q++){
        if(($num>$pos[$q]) && ($num<$pos[$q+1])){
            my $CDS_length=$num-$pos[$q];
            my @ID_strand=keys %hash_length;
            my @strand=split/\t/,$ID_strand[1];
            if($strand[1] eq "+"){
                 my $chr_pos=$hash_length{$pos[$q+1]."\t".$strand[1]}+$CDS_length-1;
                 print OUT $hash_chr{$J[0]}."\t".$j."\t".$strand[1]."\t".$chr_pos."\n";
            }else{
                 my $Chr_Pos=$hash_length{$pos[$q+1]."\t".$strand[1]}-$CDS_length+1;
                 print OUT $hash_chr{$J[0]}."\t".$j."\t".$strand[1]."\t".$Chr_Pos."\n";
            }
        }
    }
}
close IN;
close IN2;
close OUT;
