#!/usr/bin/perl
use warnings;
use strict;
use utf8;
use File::Basename;
use List::Util qw/sum/;
use List::Util qw/max min/;
use Getopt::Long;

# ## getoptions
# my ($inputdir,$outdir,$percentage_threshold);
# GetOptions (
#     "h|?" =>\&help,
#     "inputdir:s" =>\$inputdir,
#     "outdir:s" => \$outdir,
#     "percentage_threshold:s" => \$percentage_threshold,
# ) || &help;
# &help unless ($inputdir && $outdir);
# `mkdir $outdir` unless (-e $outdir);
# $percentage_threshold||="15"; 

# sub help
# {
# 	print <<"	Usage End.";
# 	Usage:
# 		--inputdir	inputfile dir   path of *arm.txt file
# 		--outdir	output dir  must be given
#       --percentage_threshold  [percentage_threshold]    default=5
# 		-h		Help    document
# 	Usage End.
# 	exit;
# }

my $inputdir = shift @ARGV;
my $outdir   = shift @ARGV;
my $percentage_threshold = shift @ARGV;



my @files = glob "$inputdir/*_arm.txt";
my %hash1;
my %hash2;

my $fo = "$outdir/arm.txt";
open my $O , '>' , $fo or die "$0 : failed to open output file '$fo': $!\n";
select $O;
my @cols;
foreach my $fi (@files) {
    open my $I , '<' , $fi or die "$0 : failed to open input file '%fi': $!\n";
    my $colname = basename($fi);
    $colname =~ s/General_//;
    $colname =~ s/_arm.txt//;
    my (%dhash,%nhash,%amhash);
    while (<$I>) {
        chomp;
        my @f = split /\t/;
        if (/^arm/) {
            push (@cols,$colname);
        }
        else {
            my $arm = $f[0];
            my $gene = $f[1];
            my $del = $f[2];
            my $nor = $f[3];
            my $amp = $f[4];
            push @{$dhash{$arm}},$del;
            push @{$nhash{$arm}},$nor;
            push @{$amhash{$arm}},$amp;
        }
    }

    foreach my $arm (sort keys %dhash) {
        my @dels = @{$dhash{$arm}};
        my @nors = @{$nhash{$arm}};
        my @amps = @{$amhash{$arm}};
        my $dn = @dels;
        my $nn = @nors;
        my $an = @amps;
        my $dmean = sprintf ("%.2f",(sum @dels)/$dn);
        my $nmean = sprintf ("%.2f",(sum @nors)/$nn);
        my $amean = sprintf ("%.2f",(sum @amps)/$an);
        my $arm_amp = $arm."+";
        my $arm_del = $arm."-";
        if ($amean < $percentage_threshold) {    ## 占比<阈值认为变化不明显
            $amean = 0;
            $hash1{$arm_amp}{$colname} = $amean;
        }
        else {
            $hash1{$arm_amp}{$colname} = $amean;
        }        
        if ($dmean < $percentage_threshold) {
            $dmean = 0;
            $hash2{$arm_del}{$colname} = $dmean;
        }
        else {
            $hash2{$arm_del}{$colname} = $dmean;
        }
    }
}
my @sort_cols = sort @cols;
my @sc = join("\t",@sort_cols);
print "arm\t@sc\n";
foreach my $arm_amp (sort keys %hash1) {
    my @array;
    foreach my $col (sort keys %{ $hash1{$arm_amp} }) {
#    foreach my $col (sort keys $hash1{$arm_amp}) {
        my $value = $hash1{$arm_amp}{$col};
        push (@array,$value);
    }
    my @os = join ("\t",@array);
    print "$arm_amp\t@os\n";
}

foreach my $arm_del (sort keys %hash2) {
    my @array2;
    foreach my $col2 (sort keys %{ $hash2{$arm_del} }) {
#    foreach my $col2 (sort keys $hash2{$arm_del}) {
        my $value2 = $hash2{$arm_del}{$col2};
        push (@array2,$value2);
    }
    my @os2 = join ("\t",@array2);
    print "$arm_del\t@os2\n";
}


my $out = "$outdir/cnv_stats.txt";
system "sort -k 1n $fo > $out";
system "rm $fo";

# sub usage {
#     print "perl general_arm_stat.pl output_dir(must be same as step1 output dir) [percentage_threshold]\n";
#     exit;
# }
