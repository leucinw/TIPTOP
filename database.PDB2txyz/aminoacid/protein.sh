#!/bin/bash
export PATH=$TINKER8C8:$PATH
RES=$1
cp template.in protein.in
rm -rf $RES*
sed "s/PPP/$RES/g" -i protein.in
protein.x < protein.in 
rm -rf *.seq *.int protein.in