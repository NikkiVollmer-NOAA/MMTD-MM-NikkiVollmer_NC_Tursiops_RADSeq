#because SNPsaurus ran every sample 3x, needed to merge all together into one bam per sample. SNPsaurus suggested to do the alignment 
#with all the fastq files first and then "sort and merge the resulting bam files into a single alignment file".
#Used the code below to do this. I didn't have a script so did each merge manually

samtools merge newname.bam bam1.bam bam2.bam bam3.bam

#after this was done did samtools depth and flagstat on the merged bams
