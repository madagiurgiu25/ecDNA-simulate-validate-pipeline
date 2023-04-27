# default environment
conda create -n snakemake python=3.7
conda activate snakemake
conda install snakemake mamba

# 1. pbsim2
git clone https://github.com/madagiurgiu25/pbsim2.git
cd  pbsim2
./configure
make
make install
export PATH=pbsim2/bin/:$PATH

# 2. cnv smurfseq
git clone https://github.com/smithlabcode/smurfseq_scripts.git
pd=`pwd`
export PATH=$pd/smurfseq_scripts/scripts:$PATH

# 3. append unicycler
git clone https://github.com/rrwick/Unicycler.git
pd=`pwd`
export PATH=$pd/Unicycler:$PATH

# binning the genome
cd scripts/cnv_data

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

python bingenome.py
# generates hg38_bins_5k_size.bed, hg38_bins_50k_size.bed
# uses as input hg38.chrom.sizes


# gcBias

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.gc5Base.bw
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz

conda install -c bioconda -c conda-forge  ucsc-bigwigtowig bedops

bigWigToWig hg38.gc5Base.bw hg38.gc5Base.wig
wig2bed --multisplit bar --keep-header  < hg38.gc5Base.wig > hg38.gc5Base.bed
rm hg38.gc5Base.bw hg38.gc5Base.wig


#deadzone genomic regions
#exclude regions with low mappability or repeats

python bingenome.py
# will generate map_hg38_500kb.wig.bed_relative.bed
# takes as  input map_hg38_500kb.wig.bed and hg38.chrom.sizes

bedtools intersect -wa -wb -a hg38_bins_5k_size.bed -b map_hg38_500kb.wig.bed_relative.bed | awk '{print $1"\t"$2"\t"$3"\t"$7}' > hg38_bins_5k_size_mappability.bed
bedtools intersect -wa -wb -a hg38_bins_50k_size.bed -b map_hg38_500kb.wig.bed_relative.bed | awk '{print $1"\t"$2"\t"$3"\t"$7}' > hg38_bins_50k_size_mappability.bed

# get indexes for all the regiosn with mappability <0.7
cat hg38_bins_5k_size_mappability.bed | awk '$4 <= 0.7 {print NR}' > hg38_bins_5k_size_exclude.txt
cat hg38_bins_50k_size_mappability.bed | awk '$4 <= 0.7 {print NR}' > hg38_bins_50k_size_exclude.txt

# nanogladiator
tar -xvf NanoGLADIATOR_1.0.tar.gz

