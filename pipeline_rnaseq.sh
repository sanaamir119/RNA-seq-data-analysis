#####check th equality of fastaqc
fastqc $1
fastqc $2
###trimming
cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out_1.fastq -p out_2.fastq $1 $2
####
mkdir STARgenome
cp refrence.fasta STARgenome
cp refernce.gtf STARgenome
###index the genome file
STAR --runMode genomeGenerate --runThreadN 2 --genomeDir STARgenome --genomeFastaFiles refrence.fasta --sjdbGTFfile refernce.gtf --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 49
####aligning
STAR --genomeDir STARgenome --runThreadN 2 --readFilesIn out_1.fastq out_2.fastq --outFileNamePrefix new_ --outFilterMultimapNmax 1 --outReadsUnmapped unmapped_new --outSAMtype BAM SortedByCoordinate
####index file from sam tools
samtools index new_Aligned.sortedByCoord.out.bam
#####count the reads
htseq-count -f bam new_Aligned.sortedByCoord.out.bam refernce.gtf > count_$3.txt
#######differential gene expression Dseq2
Rscript dseq2.R $4 $5

#####END#############
