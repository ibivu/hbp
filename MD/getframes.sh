#!/bin/sh
#cd /home/juami/LGA_FOLDPLOT_PACKAGE
mkdir frames
i=0
while [ $i -le 8000 ]
do
	#getting pdb files for frames 1-1000. For groups to be used, select "1) Protein" each time
	echo 1 | trjconv_d -f ../md_pull.trr -s ../md_pull.tpr -o frames/frame$i.pdb -b $i -e $i
	#set intervals of 20
	i=`expr $i + 20`
done
