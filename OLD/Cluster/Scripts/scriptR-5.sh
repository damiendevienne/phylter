#!/bin/sh 

vari=$1 
/panhome/comte/ecritRfile-5.sh $vari
/usr/remote/R-3.3.2/bin/R CMD BATCH --no-save /panhome/comte/Script_5_$vari.R
\rm /panhome/comte/Script_5_$vari.R
