#!/usr/bin/python
#tab2comma.py inputfile outputfile
import csv
import sys

reader = csv.reader(open(sys.argv[1], "rb"), delimiter='\t')
writer = csv.writer(open(sys.argv[2], "wb"), delimiter=',')
#uncomment next line if you want everything quoted
#writer = csv.writer(open(sys.argv[2], "wb"), delimiter='\t', quoting=csv.QUOTE_ALL)
writer.writerows(reader)