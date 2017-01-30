DATDIR= /data/DSPR

PICARD= /home/REPLAY/software/picard/build/libs/picard.jar
GATK= /home/REPLAY/software/GenomeAnalysisTK.jar

TMPINDEXFILES= \
	dmel-r5-clean.fasta.amb \
	dmel-r5-clean.fasta.ann \
	dmel-r5-clean.fasta.bwt \
	dmel-r5-clean.fasta.pac \
	dmel-r5-clean.fasta.sa

INDEXFILES= \
	A6.fasta.amb A7.fasta.amb \
	A6.fasta.ann A7.fasta.ann \
	A6.fasta.bwt A7.fasta.bwt \
	A6.fasta.pac A7.fasta.pac \
	A6.fasta.sa  A7.fasta.sa

all: $(INDEXFILES)

.INTERMEDIATE: dmel-all-chromosome-r5.57.fasta.gz
dmel-all-chromosome-r5.57.fasta.gz:
	wget ftp://flybase.org/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz

dmel-r5-clean.fasta: dmel-all-chromosome-r5.57.fasta.gz
	python clean_dmel_genome.py dmel-all-chromosome-r5.57.fasta.gz > dmel-r5-clean.fasta

%.amb %.ann %.bwt %.pac %.sa: %
	bwa index $^

dmel-r5-clean.fasta.fai: $(TMPINDEXFILES)
	samtools faidx dmel-r5-clean.fasta

dmel-r5-clean.dict: dmel-r5-clean.fasta dmel-r5-clean.fasta.fai
	java -jar $(PICARD) CreateSequenceDictionary \
		R=dmel-r5-clean.fasta \
		O=dmel-r5-clean.dict

.INTERMEDIATE: t7_all_sites_w_indels.vcf
t7_all_sites_w_indels.vcf:
	vcf-concat $(DATDIR)/t7_sites_full.vcf $(DATDIR)/t7_INDELS.vcf | \
		python clean_vcf.py > t7_all_sites_w_indels.vcf

.INTERMEDIATE: b3886_all_sites_w_indels.vcf
b3886_all_sites_w_indels.vcf:
	vcf-concat $(DATDIR)/b3886_sites_full.vcf $(DATDIR)/b3886_INDELS.vcf | \
		python clean_vcf.py > b3886_all_sites_w_indels.vcf

.INTERMEDIATE: t7_all_sites_w_indels_final.vcf
t7_all_sites_w_indels_final.vcf: t7_all_sites_w_indels.vcf dmel-r5-clean.dict
	java -jar $(PICARD) SortVcf \
		I=t7_all_sites_w_indels.vcf \
		O=t7_all_sites_w_indels_final.vcf \
		VERBOSITY=WARNING \
		SEQUENCE_DICTIONARY=dmel-r5-clean.dict

.INTERMEDIATE: b3886_all_sites_w_indels_final.vcf
b3886_all_sites_w_indels_final.vcf: b3886_all_sites_w_indels.vcf dmel-r5-clean.dict
	java -jar $(PICARD) SortVcf \
		I=b3886_all_sites_w_indels.vcf \
		O=b3886_all_sites_w_indels_final.vcf \
		VERBOSITY=WARNING \
		SEQUENCE_DICTIONARY=dmel-r5-clean.dict

A7.fasta: dmel-r5-clean.fasta t7_all_sites_w_indels_final.vcf
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel-r5-clean.fasta \
		-V t7_all_sites_w_indels_final.vcf \
		-o A7.fasta && \
	rm t7_all_sites_w_indels_final.vcf.idx


A6.fasta: dmel-r5-clean.fasta b3886_all_sites_w_indels_final.vcf
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel-r5-clean.fasta \
		-V b3886_all_sites_w_indels_final.vcf \
		-o A6.fasta && \
	rm b3886_all_sites_w_indels_final.vcf.idx

clean:
	rm -rf *.vcf *.idx *.fai *.fasta *.dict $(INDEXFILES) $(TMPINDEXFILES)
