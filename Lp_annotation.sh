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

nohup ansible-playbook -u ${USER} -i maker-hosts worker-launch.yml > log_file_2.txt 2>&1 &

#annotation failed because I didn't have snoscan file.
#trying again with just the last contig:
#In MASTER instance m1.medium
sudo chown -hR $USER /vol_b
sudo chgrp -hR $USER /vol_b
cd /vol_b
mkdir wq_maker_run
cd wq_maker_run

mkdir data
cd data
#get data
wget https://de.cyverse.org/dl/d/9B78EAE0-FD47-4CD8-B5E1-30284591498B/lytechinus_pictus_30Nov2018_OWxax.fasta
wget https://de.cyverse.org/dl/d/419D177C-3778-4C54-915E-3EF32B20724E/lp_hic_repeats.fa
wget https://de.cyverse.org/dl/d/CFD6227C-7240-4E67-8EA3-E03DDA22E556/Lv_Sp_proteins.fa
wget https://de.cyverse.org/dl/d/29AD1DB2-77B1-4056-A837-2BAAB4B4008F/LP_stringtie_transcripts.fasta

awk '/Scaffold_62804;HRSCAF=67051/{flag=1;print $0;next}/^>/{flag=0}flag' ~/lytechinus_pictus_30Nov2018_OWxax.fasta >>test.fa
cd ..

maker -CTL
rm maker_opts.ctl
curl https://raw.githubusercontent.com/warnerlab/LP_annotation/master/MAKER%20run/maker_opts.ctl > maker_opts.ctl
#changed genome to test.fa

nohup wq_maker -contigs-per-split 1 -cores 1 -memory 2048 -disk 4096 -N wq_maker1_${USER} -d all -o master.dbg -debug_size_limit=0 -stats maker_out_stats.txt > log_file.txt 2>&1 &
touch ~/.ansible.cfg
echo "[defaults]" >> ~/.ansible.cfg
echo "host_key_checking = False" >> ~/.ansible.cfg

cp /opt/WQ-MAKER_example_data/maker-hosts .
# changed hosts file!
cp /opt/WQ-MAKER_example_data/worker-launch.yml .
#changed work-launch.yml!

nohup ansible-playbook -u ${USER} -i maker-hosts worker-launch.yml > log_file_2.txt 2>&1 &


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
# -i ../round1_results/lytechinus_pictus.all.maker.transcripts.fasta -o Lp_rnd1_maker_transcripts -l /home/scijake/busco/metazoa_odb9/ \
# -m genome -c 9 -sp zebrafish -z --augustus_parameters='--progress=true' > busco.log 2>&1 &

# cp Lytechinus_pictus* /vol_b/wqmaker/augustus/config/species/Lytechinus_pictus
# cat run_Lp_rnd1_maker_transcripts/short_summary_Lp_rnd1_maker_transcripts.txt 
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


#moving on...
cd /vol_b/wqmaker/augustus/run_Lp_rnd1_maker/augustus_output/retraining_parameters
rename 's/BUSCO_Lp_rnd1_maker_3330639781/Lytechinus_pictus/g' *
sed -i 's/BUSCO_Lp_rnd1_maker_3330639781/Lytechinus_pictus/g' Lytechinus_pictus_parameters.cfg
sed -i 's/BUSCO_Lp_rnd1_maker_3330639781/Lytechinus_pictus/g' Lytechinus_pictus_parameters.cfg.orig1

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

