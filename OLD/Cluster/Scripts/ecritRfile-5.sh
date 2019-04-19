#!/bin/sh 

vari=$1 
echo "i <- $vari" > /panhome/comte/Script_5_$vari.R
cat /panhome/comte/Simul2.R >> /panhome/comte/Script_5_$vari.R 

