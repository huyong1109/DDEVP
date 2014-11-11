#! /bin/bash -f 
# This script is used to submit ensemble cases
# contains Temperature and Salinity perturbation of e-14, 4-12 and so on.
set ng0 = 1
for ngx in 1 2 4 5 8 10 16 20 32
do
for ngy in 1 2 4 8 16 32 
do
  ng0=$(( $ngx * $ngy ))
  echo $ng0
  sed -i "s/NGX0=[0-9]*/NGX0=${ngx}/g" ./resglo.h  
  sed -i "s/NGY0=[0-9]*/NGY0=${ngy}/g" ./resglo.h  
  sed -i "s/#BSUB -n [0-9]*/#BSUB -n $ng0/g" ./global.sh
  head -n 2 ./resglo.h
  head -n 2 ./global.sh
  bsub -K < ./global.sh
  wait
  
done
done
    
