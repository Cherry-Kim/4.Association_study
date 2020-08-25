import string,sys,os

print '### STEP1. Construct merged.SetID ###'
#os.system('wget ?')
fp=open('hg19.ncbiRefSeq.sort.de.txt','r')
info={}
for line in fp:
	line_temp=line[:-1].split('\t') #['1', '11874', '14409', 'DDX11L1']
	chr=line_temp[0]
	if not chr in info.keys():
		info[chr]=[]
	info[chr].append([line_temp[1],line_temp[2],line_temp[3]])
fp.close()

fp=open('merged.bim','r')
fpout=open('merged.SetID','w')
for line in fp:
	line_temp=line.split('\t') #['1', '1:10250:A:C', '0', '10250', 'C', 'A\n']
	id=line_temp[0]
	if not id in info.keys():continue
	pos=int(line_temp[3])
	val=info[id]
	for a in val:
		s=int(a[0])
		e=int(a[1])
		g=a[2]
		if s<=pos and pos<=e:
			fpout.write(g+'\t'+line_temp[1]+'\n')
fp.close()
fpout.close()

print '### STEP2. Association study using SKAT ###'
from rpy2 import robjects as ro
r = ro.r
r.source("SKAT.R")

