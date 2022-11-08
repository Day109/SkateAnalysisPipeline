import os,glob
import pandas as pd

def Trimmomatic():

	if not '1.Trimmomatic' in os.listdir('./'):os.system('mkdir 1.Trimmomatic')

	for R1 in RawDataList:

		SampleId = '_'.join(R1.split('/')[-1].split('_')[:3])

		FwRead = R1
		RvRead = R1.replace(Extension[0],Extension[1])

		PairedFw = FwRead.split('/')[-1].replace(Extension[0],'_P1.fastq.gz')
		UnpairedFw = FwRead.split('/')[-1].replace(Extension[0],'_U1.fastq.gz')
		PairedRv = RvRead.split('/')[-1].replace(Extension[1],'_P2.fastq.gz')
		UnpairedRv = RvRead.split('/')[-1].replace(Extension[1],'_U2.fastq.gz')

		Code = 'java -jar '+ProgramDir+'Trimmomatic-0.39/trimmomatic-0.39.jar '+ \
			   'PE '+ \
			   '-threads '+str(Thread)+' '+ \
			   '-phred33 '+ \
			   '-summary 1.Trimmomatic/'+SampleId+'Trimmomatic.summary '+ \
			   FwRead+' '+RvRead+' '+ \
			   '1.Trimmomatic/'+PairedFw+' 1.Trimmomatic/'+UnpairedFw+' '+ \
			   '1.Trimmomatic/'+PairedRv+' 1.Trimmomatic/'+UnpairedRv+' '+ \
			   'ILLUMINACLIP:'+ProgramDir+'Trimmomatic-0.39/adapters/'+ \
			   'TruSeq3-PE-2.fa:2:30:10 '+ \
			   'LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50'

		print(Code);os.system(Code)


def StarIndex():

	if not 'StarIndex' in os.listdir('./'):os.system('mkdir StarIndex')

	SysStr = 'STAR '+ \
			 '--runThreadN '+str(Thread)+' '+ \
			 '--runMode genomeGenerate '+ \
			 '--genomeDir StarIndex/ '+ \
			 '--genomeFastaFiles '+ReferenceFile+' '+ \
			 '--sjdbGTFfile '+GffFile
	print(SysStr);os.system(SysStr)


def Star():

	if not '2.2.Star' in os.listdir('./'):os.system('mkdir 2.2.Star')

	for R1 in sorted(glob.glob('1.Trimmomatic/*_P1.fastq.gz')):

		SampleId = '_'.join(R1.split('/')[-1].split('_')[:3])

		FwRead = R1
		RvRead = R1.replace('_P1.fastq.gz','_P2.fastq.gz')
		UnReads = [R1.replace('_P1.fastq.gz','_U1.fastq.gz'),
				   R1.replace('_P1.fastq.gz','_U2.fastq.gz')]

		SysStr = 'STAR '+ \
				 '--runThreadN '+str(Thread)+' '+ \
				 '--runMode alignReads '+ \
				 '--genomeDir StarIndex/ '+ \
				 '--readFilesIn '+FwRead+' '+RvRead+' '+ \
				 '--readFilesCommand zcat '+ \
				 '--outFileNamePrefix 2.2.Star/'+SampleId

		print(SysStr);os.system(SysStr)


def SummarizeAlignment():

	SummaryDf = {};Header = ['Sample','Total Reads','Uniquely mapped','Multple mapped','Too many mapped','unmapped-tooshort']
	for ColName in Header:SummaryDf.setdefault(ColName,[])

	LogFiles = sorted(glob.glob('2.2.Star/*Log.final.out'))
	for LogFile in LogFiles:

		SampleId = LogFile.split('/')[-1].replace('Log.final.out','')

		fin = open(LogFile,'r')
		for line in fin:
			if 'Number of input reads |' in line:
				TotalReads = line.strip().split('|')[1].split()[0]
			if 'Uniquely mapped reads % |' in line:
				Unique = line.strip().split('|')[1].split()[0]
			if '% of reads mapped to multiple loci |' in line:
				Multiple = line.strip().split('|')[1].split()[0]
			if '% of reads mapped to too many loci |' in line:
				TooMany = line.strip().split('|')[1].split()[0]
			if '% of reads unmapped: too short |' in line:
				TooShort = line.strip().split('|')[1].split()[0]
		fin.close()

		Elements = [SampleId,TotalReads,Unique,Multiple,TooMany,TooShort]
		for i in range(len(Header)):
			SummaryDf[Header[i]].append(Elements[i])
	
	SummaryDf = pd.DataFrame(SummaryDf,columns=Header)
	SummaryDf.to_csv('SummaryAlignment.txt',sep='\t')


def Samtools(Flist):
	
	for Sam in Flist:

		SortedSam = Sam.replace('Aligned.out.sam','_Sorted.sam')
		os.system('samtools sort '+Sam+' -O SAM -o '+SortedSam)
	
		## Variation 1 ##
		'''
		MapQFiltered = SortedSam.replace('_Sorted.sam','_MapQFilter.bam')
		os.system('samtools view -b -q 10 %s > %s' % (SortedSam,MapQFiltered))

		MtRemoved = SortedSam.replace('_Sorted.sam','_MapQFilter_MtRemoved_Ind.bam')
		Mt = SortedSam.replace('_Sorted.sam','_MapQFilter_Mt_Ind.bam')
		SysStr = 'samtools view -b '+ \
				 '-L MtLoci.txt '+ \
				 '-U %s %s > %s'
		os.system(SysStr % (MtRemoved,MapQFiltered,Mt))
		'''
		## Variation 2 ##

		MtRemoved = SortedSam.replace('_Sorted.sam','_Var2_MtRemoved_Ind.bam')
		Mt = SortedSam.replace('_Sorted.sam','_Var2_Mt_Ind.bam')
		SysStr = 'samtools view -b '+ \
				 '-L MtLoci.txt '+ \
				 '-U %s %s > %s'
		os.system(SysStr % (MtRemoved,SortedSam,Mt))


def Run_FeatureCount():

	Cmd = 'featureCounts -T '+str(Thread)+' '+ \
		  '-a '+GtfFile+' '+ \
		  '-o Mouse_CountDf.txt '+ \
		  '-g gene_name '+ \
		  ' 2.2.Star/*_MtRemoved_Ind.bam'

	print(Cmd)
	os.system(Cmd)


if __name__ == "__main__":

	Thread = 10
	Extension = ['_R1_001.fastq.gz','_R2_001.fastq.gz']
	ProgramDir = '/path/to/STAR_Samtools_and_etcs/'
	ReferenceFile = '/path/to/reference/Mus_musculus.GRCm38.dna.primary_assembly.fa'
	GtfFile = '/path/to/reference/Mus_musculus.GRCm38.102.gtf'
	RawDataList = glob.glob('0.RawData/*'+Extension[0])

	#Trimmomatic()

	#StarIndex()
	#Star()
	SummarizeAlignment()

	Samtools(sorted(glob.glob('2.2.Star/*Aligned.out.sam')))

	Run_FeatureCount()

	
