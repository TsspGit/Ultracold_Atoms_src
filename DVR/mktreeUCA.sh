#!/bin/bash
new_folder=$1
iz=$2
folders=`find $dir -maxdepth 2 -type 'd'`
folders2=`echo $folders | sed 's/q1d/"$1"/g' | sed 's/iz50/iz"$2"/g'`
for i in $folders2;
do 
	echo $i;
	#mkdir -p $i;
done