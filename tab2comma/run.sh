#!/bin/bash

for file in *.dat
do
	outputfile=$(echo $file | sed 's/.dat/.csv/g')
	python tab2comma.py $file ../$outputfile
done 