#!/bin/bash
export PATH=$TINKER8C8:$PATH
RES=$1
cp template.in nucleic.in
rm -rf $RES*
sed "s/PPP/$RES/g" -i nucleic.in
nucleic < nucleic.in
rm -rf *.seq *.int nucleic.in