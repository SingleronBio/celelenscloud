#!/usr/bin/perl
use warnings;
use strict;
use utf8;
use File::Basename;
use Getopt::Long;

## getoptions
my ($inputdir,$outdir,$type,$annofile);
GetOptions (
    "h|?" =>\&help,
    "inputdir:s" =>\$inputdir,
    "outdir:s" => \$outdir,
    "type:s" => \$type,
    "annofile:s" => \$annofile,
) || &help;
&help unless ($inputdir && $outdir);
`mkdir $outdir` unless (-e $outdir);
## $opt||="No"; (设置默认参数)

sub help
{
	print <<"	Usage End.";
	Usage:
		--inputdir	infercnv result dir   must be given
		--outdir	output dir  must be given
		--type	cancer celltype file or sample name file; one line split by,  must be give
        --annofile  infercnv annofile
		-h		Help    document
	Usage End.
	exit;
}


open my $C , '<' , $type or die "$0 : failed to open input file '$type' :$!\n";
my @cells=split (/\,/,<$C>);
close $C or warn "$0 : failed to close input file '$type' :$!\n";
print "@cells";

my %ahash;
open my $A , '<' , $annofile or die "$0 : failed to open input file '$annofile' : $!\n";
while (<$A>) {
    chomp;
    my @f= split /\t/;
    my $barcode = $f[0];
    my $ct = $f[1];
    if (grep /^$ct$/,@cells) {
        print "$ct\n";
        push @{$ahash{$ct}},$barcode;
    }
}
close $A or warn "$0 : failed to close input file '$annofile' :$!\n";
# my $hcl = "$inputdir/General_HCL_1_members.txt";
# if (-e $hcl) {
#     my @files = glob "$inputdir/General_HCL*";
#     foreach my $ct (sort keys %ahash) {
#         my @cellbarcodes = @{$ahash{$ct}};
#         my $fo = "$outdir/${ct}_members.txt";
#         open my $O , '>' , $fo or die "$0 : failed to open output file '$fo' :$!\n";
#         select $O;
#         my $ft = $files[0];
#         open my $T , '<' , $ft or die "$0 :failed to open input file '$ft' :$!\n";
#         my $firstLine = <$T>;
#         print "$firstLine";
#         close $T or warn "$0 : failed to close input file '$ft' :$!\n";
#         foreach my $fi (@files) {
#             open my $I , '<' , $fi or die "$0 : failed to open input file '$fi' :$!\n";
#             while (<$I>) {
#                 chomp;
#                 if ($. != 1) {
#                     my @f = split /\s+/;
#                     my $barcode = $f[0];
#                     $barcode =~ s/\"//g;
#                     if (grep /^$barcode$/, @cellbarcodes) {
#                         print "$_\n";
#                     }
#                 }
#             }
#             close $I or warn "$0 : failed to close input file '$fi' :$!\n";
#         }
#     }    
# }
# else {
    system ("Rscript /opt/lims_module_script/Heterogeneity_inferCNV/dataTransfer.R $inputdir");
    foreach my $ct (sort keys %ahash) {
        my @cellbarcodes = @{$ahash{$ct}};
        my $fo = "$outdir/${ct}_members.txt";
        open my $O , '>' , $fo or die "$0 : failed to open output file '$fo' :$!\n";
        select $O;
        my $fi = "$inputdir/infercnv.observations.references.covert.txt";
        open my $I , '<' , $fi or die "$0 : failed to open input file '$fi' :$!\n";
        while (<$I>) {
            chomp;
            if ($. == 1) {
                print "$_\n";
            }
            else {
                my @f = split /\s+/;
                my $barcode = $f[0];
                $barcode =~ s/\"//g;
                if (grep /^$barcode$/, @cellbarcodes) {
                    print "$_\n";
                }
            }
        }
        close $I or warn "$0 : failed to close input file '$fi' :$!\n";
    }
#}
