#this code uses samtools and the depth command.  If you use "samtools depth" on its own (for each sorted.bam) you will get a gigantic file
#with 3 columns: name of reference sequence, base index within the reference, depth of coverage at that base. With the code below, specifically with the 
#awk command, you are calculating the average of the 3rd column and thus the average read depth. You will get a .txt file
#(name-depth.txt) for every sorted.bam file. I am storing these txt files in a folder I created called Depth that is in the same directory as my .bam files.
#In each txt file you will see 2 numbers, the first is the average depth and the second is how many bases are in the file (e.g., 5.25475 294904660 
#where across the 294,904,660 bases that had reads mapped to them the average covereage was 5.25x; in other words, for bases that have reads mapped 
#(e.g. excluding the sites where no reads map), the average number of reads that map to these sites is 5.25. 

for sample in *.sorted.bam

do 

f1=${sample%.sorted.bam}

samtools depth $sample | awk '{sum+=$3;cnt++}END{print sum/cnt " "sum}'>Depth/$f1-depth.txt

done


