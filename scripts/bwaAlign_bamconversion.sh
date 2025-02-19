#WHAT SNPSAURUS DID: All reads were mapped to the reference with an alignment identity threshold of .95 using bbmap (BBMap tools). They used bbmap for alignment, not bwa.

#TO INDEX navigate to the folder where the reference genome is and type "bwa index pathtoref.gz"
bwa index /home/nicole.vollmer@nmfs.local/Tt_Ref/GCA_011762595.1_mTurTru1.mat.Y_genomic.fna.gz

#below is based on Ana Costa's code
#!/bin/bash
#to run script navigate to the folder with the script (am putting where my clean_fastq.gz files are) and type "nohup ./bwaAlign.sh &"
#adding & at the end sends the process to the background and nohup allows you to log out of the server and keep the process running
#all normal output to the screen will be saved in nohup.out file which you can go back and check to see if process has terminated/errored/completed


#the next line I am saying that I want to align all files that end with "_clean_fastq.gz" in the NC_Tursiops_trimmed folder and I am calling this request f1 for future use
for f1 in *_clean_fastq.gz

do

#the next line is going to rename all of the files called with f1 by removing "_clean_fastq.gz" from them all - so everything after the %
#and I am calling this request f2 for future use
f2=${f1%_clean_fastq.gz}

#below is the code for the mapping to the reference (GCA_011762595.1_mTurTru1.mat.Y_genomic.fna.gz)
#bwa index command followed by && which is used to chain commands together, such that the next command is run if and only if the preceding command exited without errors (or, more accurately, exits with a return code of 0).
#bwa mem command that calls the # of threads to use (-t), the RGline (-R) which is the readgroup line header (see more info below), the pathway to the genome file, $f1 is where to find the fastq files (see f1 line above)
#next it is calling samtools to work. When using only BWA MEM, it will generate a SAM file. So use SAMTools to convert SAM to BAM file. view: '-b' outputs in BAM format, the '-' is just to keep the same name as before,
#sort: sorts alignments by coordinates/chromosome order and creates output file name .sorted.bam (using code from f2 line),
#index: makes a compressed coordinate-sorted BAM file (.bai) for fast random-access for future programs (such as when doing variant calling) and named same as for 'sort'
#-t 20
#-R
#samtools stuff
#bwa index /home/nicole.vollmer@nmfs.local/Tt_Ref/GCA_011762595.1_mTurTru1.mat.Y_genomic.fna.gz &&      AM COMMENTING out this line because I already indexed my genome in a separate step
bwa mem -t 20 -R "@RG\tID:$f2\tPL:Illumina\tSM:$f2" /home/nicole.vollmer@nmfs.local/Tt_Ref/GCA_011762595.1_mTurTru1.mat.Y_genomic.fna.gz  $f1 | samtools view -b - | samtools sort - -o $f2.sorted.bam &&
samtools index $f2.sorted.bam

done


#Ana Costa's readgroup header line = @RG\tID:$f2\tLB:Library1\tPL:Illumina\tPU:lane1\tSM:$f2
#readgroup (RG) header line: all the \t you see will be converted to a TAB in output SAM files. The RG tags the sam/bam files. 
#It will give the sample ID to each BAM file which is useful when calling variants. 
#ID is the identifier and SM is the sample. For this code we are using the same name for both and using what we defined for f2 (see f2 line). 
#LB is DNA prep library identifier, can put something generic like 'library1' if all are from same library or can LEAVE OUT all together
#PL is platform which is ILLUMINA for SNPsaurus stuff
#PU is platform unit, can LEAVE OUT
