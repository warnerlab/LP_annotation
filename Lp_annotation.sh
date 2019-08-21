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

#In MASTER instance m1.medium
sudo chown -hR $USER /vol_b
sudo chgrp -hR $USER /vol_b
cd /vol_b
mkdir wq_maker_run
cd wq_maker_run

#get data
wget https://de.cyverse.org/dl/d/9B78EAE0-FD47-4CD8-B5E1-30284591498B/lytechinus_pictus_30Nov2018_OWxax.fasta
wget https://de.cyverse.org/dl/d/419D177C-3778-4C54-915E-3EF32B20724E/lp_hic_repeats.fa
wget https://de.cyverse.org/dl/d/CFD6227C-7240-4E67-8EA3-E03DDA22E556/Lv_Sp_proteins.fa
wget https://de.cyverse.org/dl/d/29AD1DB2-77B1-4056-A837-2BAAB4B4008F/LP_stringtie_transcripts.fasta



