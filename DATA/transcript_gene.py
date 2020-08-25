import string,sys,os

os.system('wget ?')

fp=open('hg19_refGene.txt','r')
fpout=open('hg19.ncbiRefSeq.txt','w')
for line in fp:
	line_temp=line.split('\t')
	if not ('_' in line_temp[0]):
		if 'transcript' == line_temp[2]:
			a=line_temp[8].split(';')[2]
			gene=a.split(' ')[3][1:-1]
			fpout.write(line_temp[0][3:]+'\t'+line_temp[3]+'\t'+line_temp[4]+'\t'+gene+'\n')
fp.close()
fpout.close()

os.system('sort -k1,1V -k2n hg19.ncbiRefSeq.txt > hg19.ncbiRefSeq.sort.txt')
