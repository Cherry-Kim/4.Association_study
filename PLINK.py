import string,sys,os,glob

POP='TEST'
POP1='YUHL'
POP2='YUHL_WGS'
file_list=os.listdir('/home/hykim/EWAS/YUHL')

'''
print ' ### STEP0. Delete chrX,Y,M, indel ###'
fp=glob.glob('*_sort.g.vcf')
for fname in fp:
	print '### Sample',fname
	fp1=open(fname,'r')
	sample=fname.split('_sort.g.vcf')[0]
	fpout=open(sample+'.filtred.g.vcf','w')
	for line in fp1:
		if line.startswith('#'):
			fpout.write(line)
		else:
			line_temp=line.split('\t')
			if ('chrM' != line_temp[0]) and ('chrX' != line_temp[0]) and ('chrY' != line_temp[0]):
				if (len(line_temp[3]) == 1) and (len(line_temp[4]) == 1):
					ID = line_temp[0][3:]+':'+line_temp[1]+':'+line_temp[3]+':'+line_temp[4]
					line_temp[2]=ID
					fpout.write('\t'.join(line_temp))
			
fp1.close()
fpout.close()
'''
print '### STEP1-1. Construct .bim .bed .fam files ###'
vcf_list=[file for file in file_list if file.endswith('.filtred.g.vcf')]
for i in vcf_list:
	print i
	sample=i.split('.filtred.g.vcf')[0]
	print sample
	os.system('plink2 --const-fid --vcf '+i+' --make-bed --out '+sample)

print ' ### STEP1-2. Revise .fam ###'
fp=open('merged.cov','r')
hd=fp.readline()
ID={}
for line in fp:
	line_temp=line[:-1].split('\t')
	if not line_temp[1] in ID.keys():
		ID[line_temp[1]]=line_temp[3],line_temp[4] #3:Sex 4:Pop
fp.close()

fp=glob.glob('*.fam')
for fname in fp:
	fp1=open(fname,'r')
	sample=fname.split('.fam')[0]
	fpout=open(sample+'2.fam','w')
	fpout1=open(sample+'.nosex','w')
	for line in fp1:
		line_temp=line.split(' ')
		sample=line_temp[1]
		if sample in ID.keys():
			fpout.write(line_temp[1]+' '+line_temp[1]+' '+line_temp[2]+' '+line_temp[3]+' '+ID[sample][0]+' '+ID[sample][1]+'\n')
			fpout1.write(ID[sample][0]+'\t'+line_temp[1]+'\n')
os.system('mv '+POP1+'2.fam '+POP1+'.fam')
os.system('mv '+POP2+'2.fam '+POP2+'.fam')

print ' ### STEP2. GWAS using PLINK2 ###'
os.system('plink2 --bfile '+POP1+' --bmerge '+POP2+'.bed '+POP2+'.bim '+POP2+'.fam --make-bed --out '+POP)
os.system('plink2 --bfile merged --linear --covar merged.cov --covar-name Age,Sex --pheno Protein.txt --pheno-name protein --out result --adjust')

print ' ### Haploview input files (.ped .info) ###'
os.system('plink2 --bfile '+POP+' --extract Gene_SNP_list.txt --recodeHV --out '+POP)
