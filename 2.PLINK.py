import string,sys,os,glob

POP='TEST'
POP1='TEST1'
POP2='TEST2'

print ' ### STEP1. Annotatino using SnpEff - rs ###'
os.system('java -jar /home/hykim/program/snpEff/SnpSift.jar annotate /home/hykim/REF/dbsnp_138.hg19.vcf TEST1_genotype.vcf > TEST1.dbSNP.vcf')

print ' ### STEP2-1. Delete chrX,Y,M, Indel ###'
fp=glob.glob('*dbSNP.vcf')
for fname in fp:
	print '### Sample',fname
	fp1=open(fname,'r')
	sample=fname.split('.vcf')[0]
	fpout=open(sample+'.filt.g.vcf','w')
	for line in fp1:
		if line.startswith('#'):
			fpout.write(line)
		else:
			line_temp=line.split('\t')
			if ('chrM' != line_temp[0]) and ('chrX' != line_temp[0]) and ('chrY' != line_temp[0]):
				if (len(line_temp[3]) == 1) and (len(line_temp[4]) == 1):
					fpout.write(line)

fp1.close()
fpout.close()

print '### STEP2-2. Construct .bim .bed .fam files ###'
file_list=os.listdir('/home/hykim/Rare_Var')
vcf_list=[file for file in file_list if file.endswith('.filt.g.vcf')]
for i in vcf_list:
	sample=i.split('.filt.g.vcf')[0]
	os.system('plink2 --const-fid --vcf '+i+' --make-bed --out '+sample)
os.system('plink2 --bfile '+POP1+' --bmerge '+POP2+'.bed '+POP2+'.bim '+POP2+'.fam --make-bed --out '+POP)
os.system('plink2 --bfile '+POP1+' --exclude '+POP+'-merge.missnp --make-bed --out '+POP1+'_tmp')
os.system('plink2 --bfile '+POP2+' --exclude '+POP+'-merge.missnp --make-bed --out '+POP2+'_tmp')
os.system('plink2 --bfile '+POP1+'_tmp --bmerge '+POP2+'_tmp --make-bed --out merged')
os.system('plink2 --bfile merged --freq --out merged_allele.txt')

print ' ### STEP3-1. GWAS using PLINK2 ###'
fp=open('TEST_cov.txt','r')
hd=fp.readline()
ID={}
for line in fp:
	line_temp=line[:-1].split('\t')
	if not line_temp[1] in ID.keys():
		ID[line_temp[1]]=line_temp[3],line_temp[4]
fp.close()

fp=open('merged.fam','r')
fpout=open('merged2.fam','w')
fpout1=open('merged.nosex','w')
for line in fp:
	line_temp=line.split(' ')
	sample=line_temp[1]
	if sample in ID.keys():
		fpout.write(line_temp[1]+' '+line_temp[1]+' '+line_temp[2]+' '+line_temp[3]+' '+ID[sample][0]+' '+ID[sample][1]+'\n')
		fpout1.write(ID[sample][0]+'\t'+line_temp[1]+'\n')
fp.close()
fpout.close()
os.system('mv merged2.fam merged.fam')

print ' ### STEP3-2. GWAS using PLINK2 ###'
os.system('plink2 --bfile merged --linear --covar TEST_cov.txt --covar-name sex,pop --pheno Protein.txt --pheno-name protein --out result --adjust')
