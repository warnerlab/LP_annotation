# Lp annotation collaboration with Lyons lab

#
# Trimming Illumina reads
#

#this is where we work
seq_dir=/Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/illumina_seqs
trimmo_dir=/Users/warnerj/Documents/Warner_lab/Tools/QAQC/trimmomatic/

#trimming 
java -jar $trimmo_dir/trimmomatic-0.38.jar PE -threads 4 -phred33 \
$seq_dir/Lytechinu_pictus_S19_L005_R1_001.fastq.gz \
$seq_dir/Lytechinu_pictus_S19_L005_R2_001.fastq.gz \
$seq_dir/Lp_S19_R1_trimmed_paired.fq.gz \
$seq_dir/Lp_S19_R1_trimmed_unpaired.fq.gz \
$seq_dir/Lp_S19_R2_trimmed_paired.fq.gz \
$seq_dir/Lp_S19_R2_trimmed_unpaired.fq.gz \
ILLUMINACLIP:$seq_dir/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmo_dir/trimmomatic-0.38.jar PE -threads 4 -phred33 \
$seq_dir/Lytechinu_pictus_S24_L008_R1_001.fastq.gz \
$seq_dir/Lytechinu_pictus_S24_L008_R2_001.fastq.gz \
$seq_dir/Lp_S24_R1_trimmed_paired.fq.gz \
$seq_dir/Lp_S24_R1_trimmed_unpaired.fq.gz \
$seq_dir/Lp_S24_R2_trimmed_paired.fq.gz \
$seq_dir/Lp_S24_R2_trimmed_unpaired.fq.gz \
ILLUMINACLIP:$seq_dir/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmo_dir/trimmomatic-0.38.jar PE -threads 4 -phred33 \
$seq_dir/Lytechinu_pictus_S28_L005_R1_001.fastq.gz \
$seq_dir/Lytechinu_pictus_S28_L005_R2_001.fastq.gz \
$seq_dir/Lp_S28_R1_trimmed_paired.fq.gz \
$seq_dir/Lp_S28_R1_trimmed_unpaired.fq.gz \
$seq_dir/Lp_S28_R2_trimmed_paired.fq.gz \
$seq_dir/Lp_S28_R2_trimmed_unpaired.fq.gz \
ILLUMINACLIP:$seq_dir/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#fastqc loop
for file in $seq_dir/*trimmed_paired*;
	do /Users/warnerj/Documents/Warner_lab/tools/QAQC/FastQC/fastqc "$file" --outdir=$seq_dir/fastqc_reports;
	done

###
### Stringite
###

# HISAT and STRINGTIE were on run on cyverse
# stringtie merge was run on cyvers




##
## Repeat Modeler
##

## Repeat Modeler install

# perl is already installed

## Repeat Masker install:

#Installing RMBlast
#blast from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz
#path from http://www.repeatmasker.org/isb-2.6.0+-changes-vers2.patch.gz
tar zxvf ncbi-blast-2.6.0+-src.tar.gz 
gunzip isb-2.6.0+-changes-vers2.patch.gz 
cd ncbi-blast-2.6.0+-src 
patch -p1 < ../isb-2.6.0+-changes-vers2.patch 
cd ncbi-blast-2.6.0+-src/c++ 
./configure --with-mt --prefix=/usr/local/bin/rmblast --without-debug
make
make install

#installing TRF
# install from http://tandem.bu.edu/trf/trf404.mac.download.html
cp trf404.mac-leopard /usr/local/bin/trf404.mac-leopard
chmod ug+xr trf404.mac-leopard
#Repeat databast


#intsll RECON
#from http://www.repeatmasker.org/RECON-1.08.tar.gz


## Going to try to move on with stringtie.

# need to get transcripts from the stringtie merged gtf
/Users/warnerj/Documents/Warner_lab/tools/seqtools/gffread \
-w /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fa \
-g /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/lytechinus_pictus_30Nov2018_OWxax.fasta \
/Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/merged.out.gtf


#proceeding with transdecoder workflow

# get GTF to fasta
/Users/warnerj/Documents/Warner_lab/tools/seqtools/TransDecoderv5.5.0/util/gtf_genome_to_cdna_fasta.pl \
/Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/merged.out.gtf \
/Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/lytechinus_pictus_30Nov2018_OWxax.fasta \
> /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta 

#convert to GFF3
/Users/warnerj/Documents/Warner_lab/tools/seqtools/TransDecoderv5.5.0/util/gtf_to_alignment_gff3.pl \
/Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/merged.out.gtf \
> /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.gff3

# Get best ORF
cd /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/
/Users/warnerj/Documents/Warner_lab/tools/seqtools/TransDecoderv5.5.0/TransDecoder.LongOrfs -t /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta

#count em up
grep -c ">" /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep
#   71409 transcripts
grep -o "^>MSTRG.[0-9]*." /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep | uniq | wc
#   33915   genes

/Users/warnerj/Documents/Warner_lab/tools/seqtools/TransDecoderv5.5.0/util/cdna_alignment_orf_to_genome_orf.pl \
     /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.gff3 \
     /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.gff3 \
     /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta  \
     > /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder.genome.gff3

##
## Lots of warnings:
##Warning [16125], shouldn't have a minus-strand ORF on a spliced transcript structure

#from Brian Haas https://groups.google.com/forum/#!topic/transdecoder-users/P4U4IG7Jr-s
#You can ignore those warnings.  Basically, it just means that the ORF that was predicted is found to be in the antisense orientation when taking into account the spliced orientation of the transcript sequence on the genome.

#If you run TransDecoder in strand-specific mode, then you should have ORFs that only match the transcribed orientation of the stringtie transcripts. 


#Move the files and rename them to something meaningful:
mv /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder.longest_orfs.pep
mv /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.cds /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder.longest_orfs.gff3
mv /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.gff3 /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder.longest_orfs.cds



cd /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder_dir/
makeblastdb -in /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep \
-dbtype 'prot' \
-out 'longest_orfs.pep'
-parse_seqids

/Users/warnerj/Documents/Warner_lab/tools/seqtools/gffread \
-w /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/LP_stringtie_transcripts.fa \
-g /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/lytechinus_pictus_30Nov2018_OWxax.fasta \
/Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/stringtie/merged.out.gtf

cd /Users/warnerj/Documents/Warner_lab/dbs/
makeblastdb -in /Users/warnerj/Documents/Warner_lab/dbs/uniprot_trembl.fasta \
-dbtype 'prot' \
-out 'trembl' \
-parse_seqids



## Trying RepeatModeler
/usr/local/RepeatModeler-open-1.0.11/BuildDatabase -name Lp.split.0 -engine ncbi /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/repeatmodeler/Lp.fa.split.0.fa

# Had to split up the input file into 30 chunks to avoid memory issues.
cat repeatmodeler.sh
#cd /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/split/;
#for file in *; do
#	export prefix=$(echo $file | grep -o 'Lp.fa.split.[0-9]\+' -); 
#	/usr/local/RepeatModeler-open-1.0.11/BuildDatabase -name $prefix $file ;
#	/usr/local/RepeatModeler-open-1.0.11/RepeatModeler -pa 4 -database $prefix >& $prefix.run.out;
#	echo "Done with $file";
#done

gtime --verbose repeatmodeler.sh


#Re-running some fails, and splitting the large last file
cd /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/split/;
/usr/local/RepeatModeler-open-1.0.11/RepeatModeler -pa 4 -recoverDir RM_87606.SatJul60907202019 -database Lp.fa.split.60000 >& Lp.fa.split.60000.run.out;
/usr/local/RepeatModeler-open-1.0.11/RepeatClassifier -consensi /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/split/RM_87606.SatJul60907202019/consensi.fa


/usr/local/RepeatModeler-open-1.0.11/RepeatModeler -pa 4 -database Lp.fa.split.50000 >& Lp.fa.split.50000.run.out;
/usr/local/RepeatModeler-open-1.0.11/RepeatModeler -pa 4 -database Lp.fa.split.46000 >& Lp.fa.split.46000.run.out;
/usr/local/RepeatModeler-open-1.0.11/BuildDatabase -name Lp.fa.split.8000 Lp.fa.split.8000.fa
/usr/local/RepeatModeler-open-1.0.11/RepeatModeler -pa 4 -database Lp.fa.split.8000 >& Lp.fa.split.8000.run.out;

cd /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/split/split_deeper
/Users/warnerj/Documents/Warner_lab/tools/seqtools/fasta-splitter.pl --n-parts 10 /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/split/Lp.fa.split.62000.fa

sh repeatmodeler_splitdeeper.sh

#build the library
touch lp_hic_repeats.fa
for file in /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/split/RM*/consensi.fa.classified; 
	do cat $file >> /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/lp_hic_repeats.fa
	done

for file in /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/split/split_deeper/RM*/consensi.fa.classified; 
	do cat $file >> /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/lp_hic_repeats.fa;
	done

grep -c ">" /Users/warnerj/Documents/Warner_lab/Collabs/Lp_annotation/assemblies/HiC/lp_hic_repeats.fa;


#Maker run on ATMOSPHERE


### Ok here's a real run:
ssh-keygen
cat ~/.ssh/id_rsa.pub

sudo chown -hR $USER /vol_b
sudo chgrp -hR $USER /vol_b
cd /vol_b
mkdir wq_maker_run
cd wq_maker_run

mkdir data
cd data
#get data
wget https://de.cyverse.org/dl/d/15182853-7732-479B-8832-43159FAF02B5/lytp.transcriptome.nucseq.fasta
wget https://de.cyverse.org/dl/d/419D177C-3778-4C54-915E-3EF32B20724E/lp_hic_repeats.fa
wget https://de.cyverse.org/dl/d/CFD6227C-7240-4E67-8EA3-E03DDA22E556/Lv_Sp_proteins.fa
wget https://de.cyverse.org/dl/d/29AD1DB2-77B1-4056-A837-2BAAB4B4008F/LP_stringtie_transcripts.fasta

cd ..

maker -CTL
rm maker_opts.ctl
curl https://raw.githubusercontent.com/warnerlab/LP_annotation/master/MAKER%20run/maker_opts.ctl > maker_opts.ctl

nohup wq_maker -contigs-per-split 1 -cores 1 -memory 2048 -disk 4096 -N wq_maker1_${USER} -d all -o master.dbg -debug_size_limit=0 -stats maker_out_stats.txt > log_file.txt 2>&1 &

#in each worker:
nohup work_queue_worker -N wq_maker1_${USER} --cores all --debug-rotate-max=0 -d all -o worker.dbg > log_file_2.txt 2>&1 &
work_queue_status -M wq_maker1_${USER}


#make output foloder
mkdir round1_results

#merge the outputs
nohup gff3_merge -n -s -d lytechinus_pictus_30Nov2018_OWxax.maker.output/lytechinus_pictus_30Nov2018_OWxax_master_datastore_index.log > round1_results/lytechinus_pictus.all.noseq.gff &
nohup fasta_merge -d lytechinus_pictus_30Nov2018_OWxax.maker.output/lytechinus_pictus_30Nov2018_OWxax_master_datastore_index.log &
nohup gff3_merge -s -d lytechinus_pictus_30Nov2018_OWxax.maker.output/lytechinus_pictus_30Nov2018_OWxax_master_datastore_index.log > round1_results/lytechinus_pictus.all.gff &

mv lytechinus_pictus_30Nov2018_OWxax.all.maker.proteins.fasta round1_results/lytechinus_pictus.all.maker.proteins.fasta
mv lytechinus_pictus_30Nov2018_OWxax.all.maker.transcripts.fasta round1_results/lytechinus_pictus.all.maker.transcripts.fasta
mv lytechinus_pictus_30Nov2018_OWxax.all.maker.trnascan.transcripts.fasta round1_results/lytechinus_pictus.all.maker.trnascan.transcripts.fasta

#count the gene models from round 1
cat round1_results/lytechinus_pictus.all.gff  | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
# 78275 2978.6
awk '{print $3}' round1_results/lytechinus_pictus.all.gff | sort | uniq -c

# 21189798 
#  201671 CDS
#   61299 contig
#  218324 exon
#  772492 expressed_sequence_match
#   20774 five_prime_UTR
#   78275 gene
# 1135763 match
# 7756174 match_part
#   53117 mRNA
# 1095335 protein_match
#   12988 three_prime_UTR
#   27713 tRNA

grep "^>" round1_results/lytechinus_pictus.all.maker.proteins.fasta | wc -l
# 53117
grep "^>" round1_results/lytechinus_pictus.all.maker.transcripts.fasta | wc -l
# 53117

#visualize AED
curl https://raw.githubusercontent.com/mscampbell/Genome_annotation/master/AED_cdf_generator.pl > AED_cdf_generator.pl 
perl AED_cdf_generator.pl -b 0.025 lytechinus_pictus_30Nov2018_OWxax.all.gff > round1_results/AED_round1_txt

### TRAINING SNAP
mkdir snap
cd snap
mkdir round1
cd round1
#export zff just those with AED > 0.25 and longer than 50 AA
maker2zff -x 0.25 -l 50 -d ../../lytechinus_pictus_30Nov2018_OWxax.maker.output/lytechinus_pictus_30Nov2018_OWxax_master_datastore_index.log &
rename 's/genome/Lp_rnd1.zff.length50_aed0.25/g' *

# gather some stats and validate
fathom Lp_rnd1.zff.length50_aed0.25.ann Lp_rnd1.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1 &
cat gene-stats.log

# MODEL31497 skipped due to errors
# 2834 sequences
# 0.363203 avg GC fraction (min=0.277724 max=0.498339)
# 15991 genes (plus=7917 minus=8074)
# 743 (0.046464) single-exon
# 15248 (0.953536) multi-exon
# 179.024094 mean exon (min=1 max=12686)
# 1283.358154 mean intron (min=20 max=23248)

fathom Lp_rnd1.zff.length50_aed0.25.ann Lp_rnd1.zff.length50_aed0.25.dna -validate > validate.log 2>&1 &
# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom Lp_rnd1.zff.length50_aed0.25.ann Lp_rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1 &
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1 &
# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1 &
cd ..
# assembly the HMM
hmm-assembler.pl Lp_rnd1.zff.length50_aed0.25 params > Lp_rnd1.zff.length50_aed0.25.hmm

cd ../..


## TRAINING AUGUSTUS
nohup awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' round1_results/lytechinus_pictus.all.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi data/lytechinus_pictus_30Nov2018_OWxax.fasta -bed - -fo round1_results/Lp.all.maker.transcripts1000.fasta &

#MOVED AUGUSTUS CONFIG
export AUGUSTUS_CONFIG_PATH="/home/scijake/busco/augustus_config/config/"

nohup python /home/scijake/busco/scripts/run_BUSCO.py \
-i ../round1_results/Lp.all.maker.transcripts1000.fasta -o Lp_rnd1_maker -l /home/scijake/busco/metazoa_odb9/ \
-m genome -c 9 --long -sp zebrafish -z --augustus_parameters='--progress=true' > busco.log 2>&1 &

cat ../augustus/run_Lp_rnd1_maker/short_summary_Lp_rnd1_maker.txt 
# BUSCO version is: 3.1.0 
# The lineage dataset is: metazoa_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 978)
# To reproduce this run: python /home/scijake/busco/scripts/run_BUSCO.py -i ../round1_results/Lp.all.maker.transcripts1000.fasta -o Lp_rnd1_maker -l /home/scijake/busco/metazoa_odb9/ -m genome -c 9 --long -z -sp zebrafish --augustus_parameters '--progress=true'
#
# Summarized benchmarking in BUSCO notation for file ../round1_results/Lp.all.maker.transcripts1000.fasta
# BUSCO was run in mode: genome

# C:71.4%[S:53.7%,D:17.7%],F:9.8%,M:18.8%,n:978

# 698Complete BUSCOs (C)
# 525Complete and single-copy BUSCOs (S)
# 173Complete and duplicated BUSCOs (D)
# 96Fragmented BUSCOs (F)
# 184Missing BUSCOs (M)
# 978Total BUSCO groups searched


#compare to previous:
nohup python /home/scijake/busco/scripts/run_BUSCO.py \
-i ../data/LP_stringtie_transcripts.fasta  -o Lp_stringtie -l /home/scijake/busco/metazoa_odb9/ \
-m transcriptome -c 9 -sp zebrafish -z --augustus_parameters='--progress=true' > busco.log 2>&1 &

cat run_Lp_stringtie/short_summary_Lp_stringtie.txt 
# BUSCO version is: 3.1.0 
# The lineage dataset is: metazoa_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 978)
# To reproduce this run: python /home/scijake/busco/scripts/run_BUSCO.py -i ../data/LP_stringtie_transcripts.fasta -o Lp_stringtie -l /home/scijake/busco/metazoa_odb9/ -m transcriptome -c 9 -z
#
# Summarized benchmarking in BUSCO notation for file ../data/LP_stringtie_transcripts.fasta
# BUSCO was run in mode: transcriptome

# C:89.5%[S:52.7%,D:36.8%],F:8.1%,M:2.4%,n:978

# 875Complete BUSCOs (C)
# 515Complete and single-copy BUSCOs (S)
# 360Complete and duplicated BUSCOs (D)
# 79Fragmented BUSCOs (F)
# 24Missing BUSCOs (M)
# 978Total BUSCO groups searched

nohup python /home/scijake/busco/scripts/run_BUSCO.py \
 -i ../round1_results/lytechinus_pictus.all.maker.transcripts.fasta -o Lp_rnd1_maker_transcripts -l /home/scijake/busco/metazoa_odb9/ \
 -m genome -c 9 -sp zebrafish -z --augustus_parameters='--progress=true' > busco.log 2>&1 &


cat run_Lp_rnd1_maker_transcripts/short_summary_Lp_rnd1_maker_transcripts.txt 
# # BUSCO version is: 3.1.0 
# # The lineage dataset is: metazoa_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 978)
# # To reproduce this run: python /home/scijake/busco/scripts/run_BUSCO.py -i ../round1_results/lytechinus_pictus.all.maker.transcripts.fasta -o Lp_rnd1_maker_transcripts -l /home/scijake/busco/metazoa_odb9/ -m genome -c 9 -z -sp zebrafish --augustus_parameters '--progress=true'
# #
# # Summarized benchmarking in BUSCO notation for file ../round1_results/lytechinus_pictus.all.maker.transcripts.fasta
# # BUSCO was run in mode: genome

# C:68.9%[S:48.6%,D:20.3%],F:9.7%,M:21.4%,n:978

# 674Complete BUSCOs (C)
# 475Complete and single-copy BUSCOs (S)
# 199Complete and duplicated BUSCOs (D)
# 95Fragmented BUSCOs (F)
# 209Missing BUSCOs (M)
# 978Total BUSCO groups searched

nohup python /home/scijake/busco/scripts/run_BUSCO.py \
-i ../data/lytechinus_pictus_30Nov2018_OWxax.fasta -o Lp_genome -l /home/scijake/busco/metazoa_odb9/ \
-m genome -c 9 -sp zebrafish -z --augustus_parameters='--progress=true' > busco.log 2>&1 &

cat buscos/run_Lp_genome/short_summary_Lp_genome.txt 
# BUSCO version is: 3.1.0 
# The lineage dataset is: metazoa_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 978)
# To reproduce this run: python /home/scijake/busco/scripts/run_BUSCO.py -i ../data/lytechinus_pictus_30Nov2018_OWxax.fasta -o Lp_genome -l /home/scijake/busco/metazoa_odb9/ -m genome -c 9 -z -sp zebrafish --augustus_parameters '--progress=true'
#
# Summarized benchmarking in BUSCO notation for file ../data/lytechinus_pictus_30Nov2018_OWxax.fasta
# BUSCO was run in mode: genome

# C:87.5%[S:70.4%,D:17.1%],F:4.3%,M:8.2%,n:978

# 856Complete BUSCOs (C)
# 689Complete and single-copy BUSCOs (S)
# 167Complete and duplicated BUSCOs (D)
# 42Fragmented BUSCOs (F)
# 80Missing BUSCOs (M)
# 978Total BUSCO groups searched

#moving on...
#rename the training files to useful
cd /vol_b/wqmaker/augustus/run_Lp_rnd1_maker/augustus_output/retraining_parameters
rename 's/BUSCO_Lp_rnd1_maker_3330639781/Lytechinus_pictus/g' *
sed -i 's/BUSCO_Lp_rnd1_maker_3330639781/Lytechinus_pictus/g' Lytechinus_pictus_parameters.cfg
sed -i 's/BUSCO_Lp_rnd1_maker_3330639781/Lytechinus_pictus/g' Lytechinus_pictus_parameters.cfg.orig1

cp Lytechinus_pictus* /vol_b/wqmaker/augustus/config/species/Lytechinus_pictus

#moving this back home
cp -r /home/scijake/busco/augustus_config/config /vol_b/wqmaker/augustus/
export AUGUSTUS_CONFIG_PATH="/vol_b/wqmaker/augustus/config/"

#add the trained paramters
mkdir /vol_b/wqmaker/augustus/config/species/Lytechinus_pictus

cp Lytechinus_pictus* /vol_b/wqmaker/augustus/config/species/Lytechinus_pictus

cd /vol_b/wqmaker/round1_results
# transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' lytechinus_pictus.all.noseq.gff > lytechinus_pictus_rnd1.all.maker.est2genome.gff
# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' lytechinus_pictus.all.noseq.gff > lytechinus_pictus_rnd1.all.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' lytechinus_pictus.all.noseq.gff > lytechinus_pictus_rnd1.all.maker.repeats.gff


export AUGUSTUS_CONFIG_PATH="/vol_b/wqmaker/augustus/config/"

####
#
# MAKER round 2
#
#####


#MAKER WANTS THE AUGUSTUS FILE HERE:
sudo cp -r augustus/config /vol_b/wqmaker/

curl https://raw.githubusercontent.com/warnerlab/LP_annotation/master/MAKER%20run/maker_opts_rnd2.ctl > maker_opts.ctl

nohup wq_maker -contigs-per-split 1 -cores 1 -memory 2048 -disk 4096 -N wq_maker2_${USER} -d all -o master.dbg -debug_size_limit=0 -stats maker_out_stats.txt > log_file.txt 2>&1 &

nohup work_queue_worker -N wq_maker2_${USER} --cores all --debug-rotate-max=0 -d all -o worker.dbg > log_file_2.txt 2>&1 &
work_queue_status -M wq_maker2_${USER}

mkdir round2_results

#merge the outputs
nohup gff3_merge -n -s -d lytechinus_pictus_30Nov2018_OWxax.maker.output/lytechinus_pictus_30Nov2018_OWxax_master_datastore_index.log > round2_results/lytechinus_pictus.all.noseq.gff &
nohup fasta_merge -d lytechinus_pictus_30Nov2018_OWxax.maker.output/lytechinus_pictus_30Nov2018_OWxax_master_datastore_index.log &
nohup gff3_merge -s -d lytechinus_pictus_30Nov2018_OWxax.maker.output/lytechinus_pictus_30Nov2018_OWxax_master_datastore_index.log > round2_results/lytechinus_pictus.all.gff &

mv lytechinus_pictus_30Nov2018_OWxax.all.maker.proteins.fasta round2_results/lytechinus_pictus.all.maker.proteins.fasta
mv lytechinus_pictus_30Nov2018_OWxax.all.maker.transcripts.fasta round2_results/lytechinus_pictus.all.maker.transcripts.fasta
mv lytechinus_pictus_30Nov2018_OWxax.all.maker.trnascan.transcripts.fasta round2_results/lytechinus_pictus.all.maker.trnascan.transcripts.fasta

#count the gene models from round 2
cat round2_results/lytechinus_pictus.all.gff  | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
# 61709 5807.37
awk '{print $3}' round2_results/lytechinus_pictus.all.gff | sort | uniq -c

# 22058893 
#  187238 CDS
#   61477 contig
#  224184 exon
#  432075 expressed_sequence_match
#   16290 five_prime_UTR
#   61709 gene
# 1207759 match
# 3448006 match_part
#   33069 mRNA
#  600296 protein_match
#   15016 three_prime_UTR
#   28640 tRNA

grep "^>" round2_results/lytechinus_pictus.all.maker.proteins.fasta | wc -l
# 33069
grep "^>" round2_results/lytechinus_pictus.all.maker.transcripts.fasta | wc -l
# 33069

#visualize AED
perl AED_cdf_generator.pl -b 0.025 round2_results/lytechinus_pictus.all.gff > round2_results/AED_round2.txt

### TRAINING SNAP
cd snap
mkdir round2
cd round2
#export zff just those with AED < 0.25 and longer than 50 AA
maker2zff -x 0.25 -l 50 -d ../../lytechinus_pictus_30Nov2018_OWxax.maker.output/lytechinus_pictus_30Nov2018_OWxax_master_datastore_index.log &
rename 's/genome/Lp_rnd2.zff.length50_aed0.25/g' *

# gather some stats and validate
fathom Lp_rnd2.zff.length50_aed0.25.ann Lp_rnd2.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1 &
cat gene-stats.log

# MODEL19914 skipped due to errors
# 2099 sequences
# 0.361559 avg GC fraction (min=0.297566 max=0.476190)
# 13576 genes (plus=6768 minus=6808)
# 218 (0.016058) single-exon
# 13358 (0.983942) multi-exon
# 166.917572 mean exon (min=1 max=10900)
# 1421.189575 mean intron (min=4 max=143458)

fathom Lp_rnd2.zff.length50_aed0.25.ann Lp_rnd2.zff.length50_aed0.25.dna -validate > validate.log 2>&1 &
# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom Lp_rnd2.zff.length50_aed0.25.ann Lp_rnd2.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1 &
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1 &
# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1 &
cd ..
# assembly the HMM
hmm-assembler.pl Lp_rnd2.zff.length50_aed0.25 params > Lp_rnd2.zff.length50_aed0.25.hmm

cd ../..


## TRAINING AUGUSTUS
nohup awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' round2_results/lytechinus_pictus.all.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi data/lytechinus_pictus_30Nov2018_OWxax.fasta -bed - -fo round2_results/Lp.all.maker.transcripts1000.fasta &

#Using the temporary path, this is going to get confusing
export AUGUSTUS_CONFIG_PATH="/vol_b/wqmaker/augustus/config/"

nohup python /opt/Busco/scripts/run_BUSCO.py \
-i ../round2_results/Lp.all.maker.transcripts1000.fasta -o Lp_rnd2_maker -l /opt/Busco/datasets/metazoa_obd9/ \
-m genome -c 9 --long -sp Lytechinus_pictus -z --augustus_parameters='--progress=true' > busco.log 2>&1 &

cat run_Lp_rnd2_maker/short_summary_Lp_rnd2_maker.txt

# # BUSCO version is: 3.1.0 
# # The lineage dataset is: metazoa_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 978)
# # To reproduce this run: python /opt/Busco/scripts/run_BUSCO.py -i ../round2_results/Lp.all.maker.transcripts1000.fasta -o Lp_rnd2_maker -l /opt/Busco/datasets/metazoa_obd9/ -m genome -c 9 --long -z -sp Lytechinus_pictus --augustus_parameters '--progress=true'
# #
# # Summarized benchmarking in BUSCO notation for file ../round2_results/Lp.all.maker.transcripts1000.fasta
# # BUSCO was run in mode: genome

# 	C:77.1%[S:57.8%,D:19.3%],F:7.2%,M:15.7%,n:978

# 	754	Complete BUSCOs (C)
# 	565	Complete and single-copy BUSCOs (S)
# 	189	Complete and duplicated BUSCOs (D)
# 	70	Fragmented BUSCOs (F)
# 	154	Missing BUSCOs (M)
#	978	Total BUSCO groups searched

#compare to transcriptome
nohup python /opt/Busco/scripts/run_BUSCO.py \
 -i ../round2_results/lytechinus_pictus.all.maker.transcripts.fasta -o Lp_rnd2_maker_transcripts -l /opt/Busco/datasets/metazoa_obd9/ \
 -m transcriptome -c 9 -sp Lytechinus_pictus -z --augustus_parameters='--progress=true' > busco.log 2>&1 &

# # BUSCO version is: 3.1.0 
# # The lineage dataset is: metazoa_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 978)
# # To reproduce this run: python /opt/Busco/scripts/run_BUSCO.py -i ../round2_results/lytechinus_pictus.all.maker.transcripts.fasta -o Lp_rnd2_maker_transcripts -l /opt/Busco/datasets/metazoa_obd9/ -m transcriptome -c 9 -z
# #
# # Summarized benchmarking in BUSCO notation for file ../round2_results/lytechinus_pictus.all.maker.transcripts.fasta
# # BUSCO was run in mode: transcriptome

# 	C:75.2%[S:60.3%,D:14.9%],F:11.9%,M:12.9%,n:978

# 	736	Complete BUSCOs (C)
# 	590	Complete and single-copy BUSCOs (S)
# 	146	Complete and duplicated BUSCOs (D)
# 	116	Fragmented BUSCOs (F)
# 	126	Missing BUSCOs (M)
# 	978	Total BUSCO groups searched

cd run_Lp_rnd2_maker/augustus_output/retraining_parameters/

rename 's/BUSCO_Lp_rnd2_maker_3305258666/Lytechinus_pictus/g' *
sed -i 's/BUSCO_Lp_rnd2_maker_3305258666/Lytechinus_pictus/g' Lytechinus_pictus_parameters.cfg
sed -i 's/BUSCO_Lp_rnd2_maker_3305258666/Lytechinus_pictus/g' Lytechinus_pictus_parameters.cfg.orig1

#move round 2 file to stash
cp -r /vol_b/wqmaker/config/species/Lytechinus_pictus /vol_b/wqmaker/augustus/round2_config/species/

cp Lytechinus_pictus* /vol_b/wqmaker/config/species/Lytechinus_pictus/



###
#
#
# MAKER round 3
#
#
##


#trying trinity this time:
wget https://de.cyverse.org/dl/d/15182853-7732-479B-8832-43159FAF02B5/lytp.transcriptome.nucseq.fasta 

curl https://raw.githubusercontent.com/warnerlab/LP_annotation/master/MAKER%20run/maker_opts_rnd3.ctl > maker_opts.ctl

nohup wq_maker -contigs-per-split 1 -cores 1 -memory 2048 -disk 4096 -N wq_maker3_${USER} -d all -o master.dbg -debug_size_limit=0 -stats maker_out_stats.txt > log_file.txt 2>&1 &
nohup work_queue_worker -N wq_maker3_${USER} --cores all --debug-rotate-max=0 -d all -o worker.dbg > log_file_2.txt 2>&1 &
work_queue_status -M wq_maker3_${USER}



