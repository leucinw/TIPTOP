#!/bin/bash
#export PATH=/home/liuchw/Softwares/tinkers/tinker8.7-gnu:$PATH
export PATH=$TINKER89:$PATH
RES=$1
cp template.in protein.in
rm -rf $RES*
sed "s/PPP/$RES/g" -i protein.in
protein.x < protein.in 
rm -rf *.seq *.int protein.in