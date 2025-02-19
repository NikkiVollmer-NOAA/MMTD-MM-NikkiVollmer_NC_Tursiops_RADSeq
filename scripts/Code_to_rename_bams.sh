#Use this code to rename files (bams) from SNPSaurus ID to Lab ID
#Need to make a txt (tab delimited file, here called "listfiles.txt") where first column is the SNPSaurus file name (e.g. 1842_ATGCCGCT-ACTGCATA_S119_L008_R1_001.sorted.bam) 
#and second column is the name to want for each file (e.g. 8Tt402b.bam)
#If this file is created in Excel (on PC) make sure it is opened in Notepad++ and all line breaks are removed before running code
#Save this txt file in the same place as the bams you are renaming and in terminal type in each line

while read from to
do
mv ${from}* $to
done < listfiles.txt

#If you want to first verify if the command above will do what you want before renaming the files, you can use:

while read from to
do
echo "mv ${from}* $to”
done < listfiles.txt
