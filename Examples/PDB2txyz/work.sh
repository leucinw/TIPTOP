pwd=`pwd`
for pdb in 1BNA 1EHZ 1L2Y
do
	cd $pwd
	rm -rf $pdb
	mkdir $pdb
	cp $pdb.pdb $pdb
done
