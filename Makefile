DATDIR= /data/DSPR

PICARD= /home/REPLAY/software/picard/build/libs/picard.jar
GATK= /home/REPLAY/software/GenomeAnalysisTK.jar

A6SOFT= A6Q30.._A6_A7_joint_calls.vcf
A7SOFT= A7Q30.._A6_A7_joint_calls.vcf
A6MED= A6Q30.._A6_A7_joint_calls_medium_filtered_variants.vcf
A7MED= A7Q30.._A6_A7_joint_calls_medium_filtered_variants.vcf
A6HARD= A6Q30.._A6_A7_joint_calls_hard_filtered_variants.vcf
A7HARD= A7Q30.._A6_A7_joint_calls_hard_filtered_variants.vcf

TMPINDEXFILES= \
	dmel.fasta.amb \
	dmel.fasta.ann \
	dmel.fasta.bwt \
	dmel.fasta.pac \
	dmel.fasta.sa

all: A6-soft.fasta A7-soft.fasta A6-med.fasta A7-med.fasta \
	A6-hard.fasta A7-hard.fasta

# We do not use the genome from FlyBase. There is a difference in the
# genome of the mitochondrion compared to ENSEMBL.
#	wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.18_FB2017_05/fasta/dmel-all-chromosome-r6.18.fasta.gz
	
# Get Drosophila genome.
.INTERMEDIATE: dmel.fasta
dmel.fasta:
	zcat Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz > dmel.fasta	

# Pattern rule to produce BWA indexes.
%.amb %.ann %.bwt %.pac %.sa: %
	bwa index $^

# Index the cleaned genome.
.INTERMEDIATE: dmel.fa.fai
dmel.fasta.fai: $(TMPINDEXFILES)
	samtools faidx dmel.fasta

# Use PICARD to create dictionary.
.INTERMEDIATE: dmel.dict
dmel.dict: dmel.fasta dmel.fasta.fai
	java -jar $(PICARD) CreateSequenceDictionary \
		R=dmel.fasta \
		O=dmel.dict

# Create all the fasta files.
A6-soft.fasta: dmel.fasta dmel.fasta.fai dmel.dict
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel.fasta \
		-V $(DATDIR)/$(A6SOFT) \
		-o pre-A6-soft.fasta
	python normalize_headers.py pre-A6-soft.fasta A6s > A6-soft.fasta
	#rm pre-A6-soft.fasta

A7-soft.fasta: dmel.fasta dmel.fasta.fai dmel.dict
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel.fasta \
		-V $(DATDIR)/$(A7SOFT) \
		-o pre-A7-soft.fasta
	python normalize_headers.py pre-A7-soft.fasta A7s > A7-soft.fasta
	rm pre-A7-soft.fasta

A6-med.fasta: dmel.fasta dmel.fasta.fai dmel.dict
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel.fasta \
		-V $(DATDIR)/$(A6MED) \
		-o pre-A6-med.fasta
	python normalize_headers.py pre-A6-med.fasta A6m > A6-med.fasta
	rm pre-A6-med.fasta

A7-med.fasta: dmel.fasta dmel.fasta.fai dmel.dict
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel.fasta \
		-V $(DATDIR)/$(A7MED) \
		-o pre-A7-med.fasta
	python normalize_headers.py pre-A7-med.fasta A7m > A7-med.fasta
	rm pre-A7-med.fasta

A6-hard.fasta: dmel.fasta dmel.fasta.fai dmel.dict
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel.fasta \
		-V $(DATDIR)/$(A6HARD) \
		-o pre-A6-hard.fasta
	python normalize_headers.py pre-A6-hard.fasta A6h > A6-hard.fasta
	rm pre-A6-hard.fasta

A7-hard.fasta: dmel.fasta dmel.fasta.fai dmel.dict
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel.fasta \
		-V $(DATDIR)/$(A7HARD) \
		-o pre-A7-hard.fasta
	python normalize_headers.py pre-A7-hard.fasta A7h > A7-hard.fasta
	rm pre-A7-hard.fasta

clean:
	rm -rf *.vcf *.idx *.fai *.fasta *.dict $(TMPINDEXFILES)
