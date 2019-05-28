#!/bin/bash

cd cache

for gzfl in bevpts/bev-patients*.csv.gz; do
  csvfl=$(echo $gzfl | sed 's/\.gz$//' | sed 's/bevpts\///')
  echo "Processed filename " $gzfl "to produce a csv filename called" $csvfl  
  flPrefix=$(echo $csvfl | sed 's/\.csv$//')
  echo "Processed filename " $csvfl "to produce a filename prefix called " $flPrefix 

  if [ -f "$flPrefix??.csv.gz" ] 
    then 
      echo $flPrefix "files exist, exiting this iteration"
      continue
  fi

  zcat $gzfl > $csvfl
  #mv $gzfl bk.$gzfl
  
  split -d -C 500M $csvfl $flPrefix.
  
  if [ -f "$csvfl" ] 
  then
    firstline=$(head -n 1 $csvfl)
    echo "Extracted header: " $firstline
  else 
    break 
  fi
 
  gzip $csvfl

  ls $flPrefix* | grep [0-9] | xargs -I {FL} mv {FL} {FL}.csv

  for newfile in $flPrefix*; do
    sed -i "1s/^/$firstline/" $newfile
    gzip $newfile
    echo "Finished processing" $newfile
  done
done
  

echo "Split complete" 

touch split


