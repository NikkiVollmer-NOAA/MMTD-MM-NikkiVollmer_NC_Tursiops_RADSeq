for sample in *.sorted.bam

do 

f1=${sample%.sorted.bam}

samtools flagstat $sample>Flagstat/$f1-stats.txt

done


