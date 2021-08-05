#bwa mem -P -t12 -k17 -U0 -L0,0 -T25 A6_calibration.fasta /data/Undetermined_S0_L001_R1_001.fastq.gz /data/Undetermined_S0_L001_R2_001.fastq.gz > L1_A6_clibration.sam
#bwa mem -P -t12 -k17 -U0 -L0,0 -T25 db/A7_calibration.fasta /data/Undetermined_S0_L001_R1_001.fastq.gz /data/Undetermined_S0_L001_R2_001.fastq.gz > L1_A7_clibration.sam
#bwa mem -P -t12 -k17 -U0 -L0,0 -T25 A6.fasta /data/Undetermined_S0_L001_R1_001.fastq.gz /data/Undetermined_S0_L001_R2_001.fastq.gz > L1_A6.sam
#bwa mem -P -t12 -k17 -U0 -L0,0 -T25 db/A7.fasta /data/Undetermined_S0_L001_R1_001.fastq.gz /data/Undetermined_S0_L001_R2_001.fastq.gz > L1_A7.sam
bwa mem -P -t12 -k17 -U0 -L0,0 -T25 A6+A6_calibration.fasta /data/Undetermined_S0_L001_R1_001.fastq.gz /data/Undetermined_S0_L001_R2_001.fastq.gz | head -1000000 > L1_A6+A6_clibration.sam
