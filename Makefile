all: dmel-r6-clean.fasta.fai

INDEXFILES= \
	dmel-r6-clean.fasta.amb \
	dmel-r6-clean.fasta.ann \
	dmel-r6-clean.fasta.bwt \
	dmel-r6-clean.fasta.pac \
	dmel-r6-clean.fasta.sa


dmel-r6-clean.fasta: dmel-all-chromosome-r6.13.fasta
	python clean_dmel_genome.py dmel-all-chromosome-r6.13.fasta > dmel-r6-clean.fasta

%.amb %.ann %.bwt %.pac %.sa: %
	bwa index $^

dmel-r6-clean.fasta.fai: $(INDEXFILES)
	samtools faidx dmel-r6-clean.fasta

clean:
	rm -rf $(INDEXFILES) dmel-r6-clean.fasta.fai dmel-r6-clean.fasta
