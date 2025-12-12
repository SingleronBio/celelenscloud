#!/usr/bin/perl
use warnings;
use utf8;
use strict;
use File::Basename;
use Getopt::Long;
use File::Copy;

## getoptions
my ($inputdir,$outdir,$fgene,$del_threshold,$amp_threshold,$armfile,$percentage_threshold,$cnvfilter,$armfilter);
GetOptions (
    "h|?" =>\&help,
    "inputdir:s" =>\$inputdir,
    "outdir:s" => \$outdir,
    "fgene:s" => \$fgene,
    "del_threshold:s" => \$del_threshold,
    "amp_threshold:s" => \$amp_threshold,
    "armfile:s" => \$armfile,
    "percentage_threshold:s" => \$percentage_threshold,
    "cnvfilter:s" => \$cnvfilter,
    "armfilter:s" => \$armfilter,

) || &help;
&help unless ($inputdir && $outdir);
`mkdir $outdir` unless (-e $outdir);
$del_threshold||="0.99"; 
$amp_threshold||="1.01";
$percentage_threshold||="15";
$armfile||="/Public/Script/shouhou/inferCNV/ucsc_chrom_arm/ucsc_hg38_chrom_arm.txt";
$cnvfilter||="0.5";
$armfilter||="0.9";

sub help
{
	print <<"	Usage End.";
	Usage:
		--inputdir	inputfile dir   path of *_members.txt file
		--outdir	output dir  must be given
		--fgene	gene used in HMM (HMM*.genes_used.dat)  must be give
        --del_threshold  [delation_threshold]    default=0.99
        --amp_threshold [amplition_threshold]    default=1.01
        --armfile [chromtain_arm_file] default: ucsc_hg38_chrom_arm.txt
        --percentage_threshold  [percentage_threshold]    default=15
        --cnvfilter [filter low heatmap cnv state ]   deafult=0.5
        --armfilter [filter arm ]    default=0.9
		-h		Help    document
	Usage End.
	exit;
}

my $dir = $inputdir;
my $dio = $outdir;
my $fg = $fgene;
my $fa = $armfile;
print "$fa\n";
## ucsc hg38 chrom arm file
open my $A , '<' , $fa or die "$0 : failed to open input file '$fa' : $!\n";
my %ahash;
while (<$A>) {
    chomp;
    my @f = split /\t/;
    my $chr = $f[0];
    $chr =~ s/chr//g;
    my $start = $f[1];
    my $end = $f[2];
    my $arm = $f[3];
    my $v = join ("\t",$start,$end,$arm);
    push @{$ahash{$chr}},$v;
}
## HMM used gene information
#my $fg = "$dir/HMM_CNV_predictions.HMMi6.qnorm.hmm_mode-subclusters.Pnorm_0.5.genes_used.dat";
open my $G , '<', $fg or die "$0 : failed to open input file '$fg' : $!\n";
my %ghash;
while (<$G>) {
    chomp;
    unless (/chr/) {
        my @f = split /\t/;
 	my $gene = $f[0];
        my $chr = $f[1];
        my $start = $f[2];
        my $end = $f[3];
        my $v = join ("\t",$chr,$start,$end);
        $ghash{$gene}=$v;
    }
}

### cluster files
my @files = glob "$dir/*_members.txt";
foreach my $fi (@files) {
    open my $I , '<' , $fi or die "$0 : failed to open input file '$fi' : $!\n";
    my $bn = basename($fi);
    $bn =~ s/members.txt/arm.txt/;   
    my $fo = "$dio/$bn";
    open my $O , '>' , $fo or die "$0 : failed to open output file '$fo' : $!\n";
    select $O ;
    print "arm\tgene\tdeletion\tnormal\tamplification\n";
    my %hash;
    my @genes;
    while (<$I>) {
        chomp;
        if ($. == 1) {
            @genes = split (/\s+/,$_);  
        }
        else {
            my @f = split (/\s+/,$_);
            for (my $j = 0; $j < $#f; $j++) {
                my $col = $j+1;
                push @{$hash{$genes[$j]}}, $f[$col];
            }
        }
    }
    foreach my $gene (sort keys %hash) {
        my $del = 0;
        my $nor = 0;
        my $amp = 0;
        my @vs = @{$hash{$gene}};
        $gene =~ s/"//g;
        foreach my $v (@vs) {
            if ($v <= $del_threshold) {
                $del = $del + 1;
            }
            elsif ($v >= $amp_threshold) {
                $amp = $amp + 1;
            }
            else {
                $nor = $nor + 1;
            }
        }
        my $total = $del + $nor + $amp;
        my $d = sprintf ("%.4f",$del/$total) * 100;
        my $n = sprintf ("%.4f",$nor/$total) * 100;
        my $a = sprintf ("%.4f",$amp/$total) * 100;
        if (exists $ghash{$gene}) {
            my $v = $ghash{$gene};
            my @p = split (/\t/,$v);
            my $chr = $p[0];
            my $start = $p[1];
            my $end = $p[2];
            if (exists $ahash{$chr}) {
                my @arms = @{$ahash{$chr}};
                foreach my $arm (@arms) {
                    my @a = split (/\t/,$arm);
                    my $as = $a[0];
                    my $ae = $a[1];
                    my $ca = $a[2];
                    if ($start > $as && $end < $ae) {
                        my $chrarm = $chr.$ca;
                        print "$chrarm\t$gene\t$d\t$n\t$a\n";
                    }
                }
            }
        }        
    }
    close $I or warn "$0 : failed to close input file '$fi' :$!\n";
    close $O or warn "$0 : failed to close output file '$fo' :$!\n";
}

# sub usage {
#     print "perl gene_stats.pl file_input_dir output_dir HMM*.genes_used.dat [delation_threshold] [amplition_threshold] [chromtain_arm_file]\n";
#     exit;
# }


system "perl /opt/lims_module_script/Heterogeneity_inferCNV/general_arm_stat.pl $dio $dio $percentage_threshold";
my $outh = "$dio/heatmap";
mkdir $outh unless -d $outh;
system "Rscript /opt/lims_module_script/Heterogeneity_inferCNV/heatmap.R $dio $outh $cnvfilter $armfilter";


