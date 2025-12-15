for (( i =2;  i <=7; i++ ))
  do 
samtools view --threads 24 -S -b SRR2629708$i.sam  >  SRR2629708$i.bam
done
