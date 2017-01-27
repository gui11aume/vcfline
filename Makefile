all: dmel-r5-clean.fasta.fai dmel-r5-clean.dict X \
	t7_all_sites_w_indels.vcf b3886_all_sites_w_indels.vcf

DATDIR= /data/DSPR

PICARD= /home/REPLAY/software/picard/build/libs/picard.jar
GATK= /home/REPLAY/software/GenomeAnalysisTK.jar

INDEXFILES= \
	dmel-r5-clean.fasta.amb \
	dmel-r5-clean.fasta.ann \
	dmel-r5-clean.fasta.bwt \
	dmel-r5-clean.fasta.pac \
	dmel-r5-clean.fasta.sa

.INTERMEDIATE: dmel-all-chromosome-r5.57.fasta.gz
dmel-all-chromosome-r5.57.fasta.gz:
	wget ftp://flybase.org/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz

dmel-r5-clean.fasta: dmel-all-chromosome-r5.57.fasta.gz
	python clean_dmel_genome.py dmel-all-chromosome-r5.57.fasta.gz > dmel-r5-clean.fasta

%.amb %.ann %.bwt %.pac %.sa: %
	bwa index $^

dmel-r5-clean.fasta.fai: $(INDEXFILES)
	samtools faidx dmel-r5-clean.fasta

dmel-r5-clean.dict: dmel-r5-clean.fasta
	java -jar $(PICARD) CreateSequenceDictionary \
		R=dmel-r5-clean.fasta O=dmel-r5-clean.dict

t7_all_sites_w_indels.vcf:
	vcf-concat $(DATDIR)/t7_sites_full.vcf $(DATDIR)/t7_INDELS.vcf | \
		python clean_vcf.py | vcf-sort > t7_all_sites_w_indels.vcf

b3886_all_sites_w_indels.vcf:
	vcf-concat $(DATDIR)/b3886_sites_full.vcf $(DATDIR)/b3886_INDELS.vcf | \
		python clean_vcf.py | vcf-sort > b3886_all_sites_w_indels.vcf

X: dmel-r6-clean.fasta t7_all_sites_w_indels.vcf
	java -jar $(GATK) \
		-T FastaAlternateReferenceMaker \
		-R dmel-r5-clean.fasta \
		-V t7_all_sites_w_indels.vcf \
		-o t7_reference.fasta

clean:
	rm -rf *vcf $(INDEXFILES) dmel-r5-clean.fasta.fai dmel-r5-clean.fasta
