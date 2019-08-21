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

#round one:
cat maker_opts.ctl
# #-----Genome (these are always required)
# genome=./data/lytechinus_pictus_30Nov2018_OWxax.fasta  #genome sequence (fasta file or fasta embeded in GFF3 file)
# organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

# #-----Re-annotation Using MAKER Derived GFF3
# maker_gff= #MAKER derived GFF3 file
# est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
# altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
# protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
# rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
# model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
# pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
# other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

# #-----EST Evidence (for best results provide a file for at least one)
# est=./data/LP_stringtie_transcripts.fasta #set of ESTs or assembled mRNA-seq in fasta format
# altest= #EST/cDNA sequence file in fasta format from an alternate organism
# est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
# altest_gff= #aligned ESTs from a closly relate species in GFF3 format

# #-----Protein Homology Evidence (for best results provide a file for at least one)
# protein=./data/Lv_Sp_proteins.fa  #protein sequence file in fasta format (i.e. from mutiple oransisms)
# protein_gff=  #aligned protein homology evidence from an external GFF3 file

# #-----Repeat Masking (leave values blank to skip repeat masking)
# model_org= #select a model organism for RepBase masking in RepeatMasker
# rmlib=./data/lp_hic_repeats.fa #provide an organism specific repeat library in fasta format for RepeatMasker
# repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
# rm_gff= #pre-identified repeat elements from an external GFF3 file
# prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
# softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

# #-----Gene Prediction
# snaphmm= #SNAP HMM file
# gmhmm= #GeneMark HMM file
# augustus_species= #Augustus gene prediction species model
# fgenesh_par_file= #FGENESH parameter file
# pred_gff= #ab-initio predictions from an external GFF3 file
# model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
# est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
# protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
# trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
# snoscan_rrna=./test_data/Os-rRNA.fa #rRNA file to have Snoscan find snoRNAs
# unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

# #-----Other Annotation Feature Types (features MAKER doesn't recognize)
# other_gff= #extra features to pass-through to final MAKER generated GFF3 file

# #-----External Application Behavior Options
# alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
# cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

# #-----MAKER Behavior Options
# max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
# min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

# pred_flank=200 #flank for extending evidence clusters sent to gene predictors
# pred_stats=0 #report AED and QI statistics for all predictions as well as models
# AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
# min_protein=0 #require at least this many amino acids in predicted proteins
# alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
# always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
# map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
# keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

# split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
# single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
# single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
# correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

# tries=2 #number of times to try a contig if there is a failure for some reason
# clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
# clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
# TMP= #specify a directory other than the system default temporary directory for temporary files

nohup wq_maker -contigs-per-split 5234 -cores 1 -memory 2048 -disk 4096 -N wq_maker1_${USER} -d all -o master.dbg -debug_size_limit=0 -stats test_out_stats.txt > log_file.txt 2>&1 &

#make this:
cat ~/.ansible.cfg
[defaults]
host_key_checking = False

cp /opt/WQ-MAKER_example_data/maker-hosts .
echo "149.165.168.125" >> maker-hosts
echo "149.165.168.217" >> maker-hosts
echo "149.165.169.93" >> maker-hosts
echo "149.165.168.213" >> maker-hosts
echo "149.165.169.63" >> maker-hosts
echo "149.165.168.185" >> maker-hosts
echo "149.165.168.134" >> maker-hosts
echo "149.165.168.238" >> maker-hosts
echo "149.165.169.41" >> maker-hosts
echo "149.165.168.138" >> maker-hosts
echo "149.165.168.223" >> maker-hosts
echo "149.165.168.221" >> maker-hosts

cp /opt/WQ-MAKER_example_data/worker-launch.yml .
---
- hosts : workers
  environment:
    PATH: "{{ ansible_env.PATH }}:/home/upendra/bin:/home/upendra/.local/bin:/opt/icommands:/opt/icommands:/opt/exonerate-2.2.0-x86_64/bin/:/opt/cctools/bin:/opt/ncbi-blast-2.6.0+/bin/:/opt/snoscan-0.9.1/:/opt/tRNAscan-SE-1.3.1/:/opt/snap/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/opt/augustus-3.2.2/bin:/opt/maker/bin:/opt/RepeatMasker:/opt/snap"    
    PERL5LIB: "/opt/tRNAscan-SE-1.3.1::/opt/cctools/lib/perl5/site_perl"
  tasks :
  - name : Execute the script
    shell : /opt/cctools/bin/work_queue_worker -N wq_maker1_${USER} -s /home/${USER} --cores all --debug-rotate-max=0 -d all -o /home/${USER}/worker.dbg

