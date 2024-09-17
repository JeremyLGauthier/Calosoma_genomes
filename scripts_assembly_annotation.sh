#! /bin/bash

############################################################
# Full script to perform the assembly and annotation
#Used for the process of C. wilkesii and C. tepidum genomes
#
# Jeremy Gauthier Museum of Geneva 2024
############################################################

#PacBio data assembly

#extract data
pbindex $sample.bam
bam2fastq -o $sample $sample.bam

#Hifiasm first tool tested
hifiasm -o $sample.asm -t 4 $sample.fastq.gz

#fasta of the final assembly
awk '/^S/{print ">"$2;print $3}' $sample.asm.bp.hap1.p_ctg.gfa > $sample.asm.bp.hap1.p_ctg.fasta

#purge dup
minimap2 -ax map-hifi $sample.asm.bp.hap1.p_ctg.fasta "$i" | gzip -c - > temp_file.paf.gz
pbcstat temp_file.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log
split_fa $sample.asm.bp.hap1.p_ctg.fasta > $sample.asm.bp.hap1.p_ctg.fasta_split
minimap2 -xasm5 -DP $sample.asm.bp.hap1.p_ctg.fasta_split $sample.asm.bp.hap1.p_ctg.fasta_split | gzip -c - > pri_asm.split.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed $sample.asm.bp.hap1.p_ctg.fasta
mv purged.fa purged_$sample.asm.bp.hap1.p_ctg.fasta

#Kmer evaluation
mkdir tmp
kmc -k21 -t2 -m8 -ci1 -cs10000 $sample.fastq.gz reads tmp/
kmc_tools transform reads histogram reads.histo -cx10000
genomescope.R -i reads.histo -o output_dir -k 21

#HiC assembly
bwa mem purged_$sample.asm.bp.hap1.p_ctg.fasta $sample.HiC-reads_R1.fastq.gz | samtools view -bh - | filter-chimeras.py - > r1.bam
bwa mem purged_$sample.asm.bp.hap1.p_ctg.fasta $sample.HiC-reads_R2.fastq.gz | samtools view -bh - | filter-chimeras.py - > r2.bam
combine_ends.py r1.bam r2.bam | samtools fixmate -m - - | samtools sort - | samtools markdup -r - combined.bam

yahs purged_$sample.asm.bp.hap1.p_ctg.fasta combined.bam

#Busco
busco -i $sample_yahs.out_scaffolds_final.fa -o busco_out -m genome -l insecta_odb10 

#Cov
minimap2 -ax $sample_yahs.out_scaffolds_final.fa $sample.fastq.gz > pb_on_ref.sam
samtools sort pb_on_ref.sam -o pb_on_ref_sort.bam
samtools coverage pb_on_ref_sort.bam

#gc content
source activate envemboss
infoseq $sample_yahs.out_scaffolds_final.fa

#blast
pyfasta split -n 1 -k 100000 $sample_yahs.out_scaffolds_final.fa
grep ">" $sample_yahs.out_scaffolds_final.split.100Kmer.fa | sed -e 's/>//g' > list_seq
while read a
	do
	echo $a > temp_seq
	fastaselect.pl $sample_yahs.out_scaffolds_final.split.100Kmer.fa temp_seq > temp_seq.fasta
	diamond blastx --query temp_seq.fasta --threads 2 --max-target-seqs 10 --db uniprot_sprot.fasta.dmnd --evalue 1e-25 --outfmt 6 --out temp_blast
	cat final_out temp_blast >> final_out
	done < list_seq

#Annotation

#masking
RepeatMasker $sample_yahs.out_scaffolds_final.fa -xsmall

#proteins mapping

prothint.py $sample_yahs.out_scaffolds_final.fa.masked reference_proteines.fasta
#output: prothint_augustus.gff

#braker on prot
singularity exec braker3.sif braker.pl --AUGUSTUS_CONFIG_PATH=/local/Augustus/config --AUGUSTUS_BIN_PATH=/local/Augustus/bin --AUGUSTUS_SCRIPTS_PATH=/local/Augustus/scripts --GENEMARK_PATH=/local/gmes_linux_64 --genome=$sample_yahs.out_scaffolds_final.fa.masked --threads=2 --species=$sample --hints=prothint_augustus.gff

#outputs: hintsfile.gff and augustus.hints.gtf

#long reads mapping

source activate envisoseq3
isoseq cluster flnc.fofn clustered.bam --verbose --use-qvs

minimap2 -ax splice:hq -uf $sample_yahs.out_scaffolds_final.fa.masked clustered.hq.fasta.gz > rnaseq_on_chrom.sam

sort -k 3,3 -k 4,4n rnaseq_on_chrom.sam > rnaseq_on_chrom.s.sam

source activate envcdna_cupcake2
collapse_isoforms_by_sam.py --input clustered.hq.fasta -s rnaseq_on_chrom.s.sam --dun-merge-5-shorter -o cupcake

source activate envaugustus3
stringtie2fa.py -g $sample_yahs.out_scaffolds_final.fa.masked -f cupcake.collapsed.gff -o cupcake.fa
gmst.pl --strand direct cupcake.fa.mrna --output gmst.out --format GFF
gmst2globalCoords.py -t cupcake.collapsed.gff -p gmst.out -o gmst.global.gtf -g $sample_yahs.out_scaffolds_final.fa.masked

#output: gmst.global.gtf

#COMBINING ALL with TSEBRA

tsebra.py -g augustus.hints.gtf -e hintsfile.gff -l gmst.global.gtf -c /home/jeremy/local/TSEBRA/config/long_reads.cfg -o tsebra.gtf
#inputs : augustus.hints.gtf and hintsfile.gff from first braker ; gmst.global.gtf from long read mapping
#output: tsebra.gtf

#extract prot and cds
source activate envgffread
gffread -g chrom_purged_231220_64366e_A01_AUGM_10.asm.bp.hap1.p_ctg.gfa.fasta.masked -x cds.fa --gtf tsebra.gtf -y prot.fa

#Orthofinder
orthofinder -b prot_folder


