for (( i=2; i <=7; i++ ))
  do 
fasterq-dump --threads 24 --verbose -t /archive/alotaibih/sehribanb/temp --outdir /archive/alotaibih/sehribanb/mouse_normal/d1 SRR2629708$i
done
