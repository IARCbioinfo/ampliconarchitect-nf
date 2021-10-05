


my $bed=$ARGV[0];
#my $bam=$ARGV[];
my $sample=$ARGV[1];

my $cram=`grep $sample MESOMICSID_genomic.rename.bam.tsv | awk '{print \$2}'  | xargs basename | sed 's/bam/cram/'`;
chomp($cram);
my $bam="/data/gcs/mesomics/files/WGS/cram/".$cram;
my $AA="/home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py";
my $ref="/home/digenovaa/scratch/mesomics/AMPARCHITECT/GRCh38";
my $file="coverage.stats";
my $root=`pwd`;
chomp $root;

print $cram."\n";
# we count the number of segments

my $n_s=`wc -l $bed | awk '{print \$1}'`;

if($n_s == 0 or $n_s > 70){
	print "Incorrect number of seeds incorrect $n_s\n";
	system("touch $sample.aa.done");
	exit(0);
}
my $dir=$sample.".aa";
my $repo=$root."/".$dir;

my $script = <<EOF;
mkdir -p $dir
cd $dir
ln -s $ref .
#we export the main variables
MOSEKLM_LICENSE_FILE=/data/scratch/digenovaa/mesomics/AMPARCHITECT/mosek/mosek.lic
export MOSEKLM_LICENSE_FILE
export AA_DATA_REPO=$repo
touch $file
python2 /home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py --bed $root/$bed --bam $bam --out $sample --ref GRCh38
cd ..
touch $sample.aa.done
EOF

open(FILE, ">".$sample.".run_aa.sh") or die "cannot create file";


print  FILE $script;
my $cmd="sh ".$sample.".run_aa.sh";
system($cmd);



