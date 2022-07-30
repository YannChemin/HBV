import csv
import numpy as np
qo = np.genfromtxt('dischargeobs.dat')
print qo
print qo.shape
#reader = csv.reader(open('dischargeobs.dat', 'rb'), delimiter='\t', quotechar='|')
#for line in reader:
#	print ', '.join(line)
#print reader.line_num
btot=8
print btot
#qo=np.zeros((reader.line_num,btot), dtype=float)
#print qo
#for line in reader:
#	if(col_num >= btot):
#		col_num=0
#	qo[reader.line_num][1]=line
#	col_num = col_num + 1

#print qo
