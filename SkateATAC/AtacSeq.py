import os, glob
import pandas as pd

def Trimmomatic():

	if not '1.Trimmomatic' in os.listdir('./'):os.system('mkdir 1.Trimmomatic')

	for R1 in RawDataList:

		SampleId = '_'.join(R1.split('/')[-1].split('_')[2:5])

		FwRead = R1
		RvRead = R1.replace('R1_001.fastq.gz','R2_001.fastq.gz')

		PairedFw = FwRead.split('/')[-1].replace('_R1_001.fastq.gz','_P1.fastq.gz')
		UnpairedFw = FwRead.split('/')[-1].replace('_R1_001.fastq.gz','_U1.fastq.gz')
		PairedRv = RvRead.split('/')[-1].replace('_R2_001.fastq.gz','_P2.fastq.gz')
		UnpairedRv = RvRead.split('/')[-1].replace('_R2_001.fastq.gz','_U2.fastq.gz')

		Code = 'java -jar '+ProgramDir+'Trimmomatic-0.39/trimmomatic-0.39.jar '+ \
			   'PE '+ \
			   '-threads '+str(Thread)+' '+ \
			   '-phred33 '+ \
			   '-summary 1.Trimmomatic/'+SampleId+'Trimmomatic.summary '+ \
			   FwRead+' '+RvRead+' '+ \
			   '1.Trimmomatic/'+PairedFw+' 1.Trimmomatic/'+UnpairedFw+' '+ \
			   '1.Trimmomatic/'+PairedRv+' 1.Trimmomatic/'+UnpairedRv+' '+ \
			   'ILLUMINACLIP:'+ProgramDir+'Trimmomatic-0.39/adapters/'+ \
			   'NexteraPE-PE.fa:2:30:10 '+ \
			   'LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30 '+ \
			   '2> 1.Trimmomatic/'+SampleId+'.log'

		print(Code);os.system(Code)


def BwaIndex():

	Code = ProgramDir+'bwa/bwa index '+ReferenceFile
	print(Code);os.system(Code)


def Bwa():

	if not '2.Bwa' in os.listdir('./'):os.system('mkdir 2.Bwa')

	for R1 in sorted(glob.glob('1.Trimmomatic/*_P1.fastq.gz')):

		SampleId = '_'.join(R1.split('/')[-1].split('_')[2:5])
		RunId = '_'.join(R1.split('/')[-1].split('_')[2:6])

		FwRead = R1
		RvRead = R1.replace('_P1.fastq.gz','_P2.fastq.gz')

		Code = ProgramDir+'bwa/bwa mem '+ \
			   '-M -t '+str(Thread)+' '+ \
			   '-R "@RG\\tID:'+RunId+'\\tPL:ILLUMINA\\tLB:TRUSEQ\\tSM:'+SampleId+'" '+ \
			   ReferenceFile+' '+ \
			   FwRead+' '+RvRead+ \
			   ' > 2.Bwa/'+RunId+'.sam'+ \
			   ' 2> 2.Bwa/'+RunId+'.log'

		print(Code);os.system(Code)


		BamFile = '2.Bwa/'+RunId+'.bam'
		Code = 'samtools view -S -b 2.Bwa/'+RunId+'.sam > '+BamFile
		print(Code);os.system(Code)
		os.system('rm 2.Bwa/'+RunId+'.sam')

		Code = 'samtools sort '+BamFile+' -o '+BamFile.replace('.bam','Sorted.bam')
		print(Code);os.system(Code)

		Code = 'samtools index '+BamFile.replace('.bam','Sorted.bam')
		print(Code);os.system(Code)


def MACS():

	if not '3.Macs' in os.listdir('./'):os.system('mkdir 3.Macs')

	#BamList = sorted(glob.glob('2.Bwa/*Sorted.bam'))
	#BamList = sorted(glob.glob('2.Bwa/*_DupRmv.bam'))
	BamList = sorted(glob.glob('2.Bwa/*_MapQFilter_MtRemoved.bam'))
	#BamList = sorted(glob.glob('2.Bwa/*_MapQFilter_MtRemoved_Ind.bam'))
	for BamFile in BamList:
		
		SampleId = BamFile.split('/')[-1].replace('_MapQFilter_MtRemoved.bam','')

		Code = 'macs2 callpeak '+ \
			   '-t '+BamFile+' '+ \
			   '-n '+SampleId+' '+ \
			   '--outdir 3.Macs/ '+ \
			   '-q 0.05 '+ \
			   '--nomodel '+ \
			   '-f BAMPE '+ \
			   '-g 1974810099 '+ \
			   '2> 3.Macs/'+SampleId+'_Cmd.log'

		print(Code);os.system(Code)


def RemoveHeader():

	FileList = sorted(glob.glob('3.Macs/*_peaks.xls'))
	for FileName in FileList:

		Sample = '_'.join(FileName.split('/')[-1].split('_')[:3])
		OutFile = '3.Macs/'+Sample+'_peaks.bed'

		fout = open(OutFile,'w')
		fin = open(FileName,'r')
		for line in fin:

			if '#' in line[0]:
				pass
			elif len(line.strip()) == 0:
				pass
			elif 'chr' == line.split('\t')[0]:
				pass
			else:
				fout.write(line)
		fin.close();fout.close()


def MergeBed():

	RemoveHeader()

	CatFile = '3.Macs/ConcatPeaks.bed'
	SortedCatFile = '3.Macs/ConcatPeaks_Sorted.bed'
	MergedFile = '3.Macs/MergedPeaks.bed'

	BedList = sorted(glob.glob('3.Macs/*_peaks.bed'))
	SampleList = [i.split('/')[-1].replace('_peaks.bed','') for i in BedList]

	## Sort chromosome, position ##

	SortedList = []
	for BedFile in BedList:
		
		SortedFile = BedFile.replace('_peaks.bed','_sorted.bed')
		SortedList.append(SortedFile)
		SysStr = 'sort -k1,1 -k2,2n %s > %s' % (BedFile,SortedFile)
		os.system(SysStr)

	## Concat BedFiles ##

	SysStr = 'cat '+' '.join(SortedList)+' > '+CatFile
	os.system(SysStr)

	## Sort Concatenated BedFile ##

	SysStr = 'sort -k1,1 -k2,2n %s > %s' % (CatFile,SortedCatFile)
	os.system(SysStr)

	## Merge Sorted Concatenated BedFile ##

	SysStr = 'bedtools merge -i %s > %s' % (SortedCatFile,MergedFile)
	os.system(SysStr)


	## Parse 3'Utr and Genic regions ##

	UtrDfFile = '../Reference/UtrDf.txt'
	GenicDfFile = '../Reference/GenicDf.txt'
	GeneDfFile = '../Reference/GeneDf.txt'

	MakeUtrCoDf(GffFile,UtrDfFile)
	MakeGenicDf(GffFile,GenicDfFile)
	MakeGeneCoDf(GffFile,GeneDfFile)

	SortedUtrDf = UtrDfFile.replace('.txt','_Sorted.txt')
	SortedGenicDf = GenicDfFile.replace('.txt','_Sorted.txt')
	SortedGeneDf = GeneDfFile.replace('.txt','_Sorted.txt')

	SysStr = 'sort -k1,1 -k2,2n %s > %s' % (UtrDfFile,SortedUtrDf)
	os.system(SysStr)

	SysStr = 'sort -k1,1 -k2,2n %s > %s' % (GenicDfFile,SortedGenicDf)
	os.system(SysStr)

	SysStr = 'sort -k1,1 -k2,2n %s > %s' % (GeneDfFile,SortedGeneDf)
	os.system(SysStr)

	## Intersect 3'Utr with Peaks ##

	UtrMergedBed = '3.Macs/MergedPeaks_utr.bed'
	SysStr = 'bedtools intersect '+ \
			 '-wa -wb '+ \
			 '-a '+MergedFile+' '+ \
			 '-b '+SortedUtrDf+' '+ \
			 '> '+UtrMergedBed
	os.system(SysStr)

	MergedUtrDf = pd.read_csv(UtrMergedBed,header=None,sep='\t')
	MergedUtrDf = MergedUtrDf[[0,1,2,6,7,8]]

	## Intersect Genic regions with Peaks ##

	GenicMergedBed = '3.Macs/MergedPeaks_Genic.bed'
	SysStr = 'bedtools intersect '+ \
			 '-wa -wb '+ \
			 '-a '+MergedFile+' '+ \
			 '-b '+SortedGenicDf+' '+ \
			 '> '+GenicMergedBed
	os.system(SysStr)

	MergedGenicDf = pd.read_csv(GenicMergedBed,header=None,sep='\t')
	MergedGenicDf = MergedGenicDf[[0,1,2,6,7,8]]

	## Intersect Gene coordinates with Peaks ##

	GeneMergedBed = '3.Macs/MergedPeaks_Gene.bed'
	SysStr = 'bedtools intersect '+ \
			 '-wa -wb '+ \
			 '-a '+MergedFile+' '+ \
			 '-b '+SortedGeneDf+' '+ \
			 '> '+GeneMergedBed
	os.system(SysStr)

	MergedGeneDf = pd.read_csv(GeneMergedBed,header=None,sep='\t')
	MergedGeneDf = MergedGeneDf[[0,1,2,6,7,8]]

	## Merge 3'Utr & Genic region & gene coordinates ##
	
	AnnotatedBed = '3.Macs/MergedPeaks_Annotated.bed'
	AnnotatedDf = MergedUtrDf.append(MergedGenicDf).sort_values(by=[0,1,2])
	AnnotatedDf = AnnotatedDf.append(MergedGeneDf).sort_values(by=[0,1,2])
	AnnotatedDf.to_csv(AnnotatedBed,sep='\t',index=False,header=False)

	AnnotatedMerged = '3.Macs/MergedPeaks_Annotated_Merged.bed'
	SysStr = 'bedtools merge '+ \
			 '-c 4,5,6 '+ \
			 '-o distinct,distinct,distinct '+ \
			 '-i %s > %s'
	os.system(SysStr % (AnnotatedBed,AnnotatedMerged))

	## Intersect each sample ##

	IntersectFile = '3.Macs/IntersectPeaks.bed'
	SysStr = 'bedtools intersect '+ \
			 '-wa -wb '+ \
			 '-a '+MergedFile+' '+ \
			 '-b '+' '.join(BedList)+' '+ \
			 '-names '+' '.join(SampleList)+' '+ \
			 '> '+IntersectFile
	os.system(SysStr)

	IntersectFile = '3.Macs/IntersectPeaks_Annotated.bed'
	SysStr = 'bedtools intersect '+ \
			 '-wa -wb '+ \
			 '-a '+AnnotatedMerged+' '+ \
			 '-b '+' '.join(BedList)+' '+ \
			 '-names '+' '.join(SampleList)+' '+ \
			 '> '+IntersectFile
	os.system(SysStr)


def MakeCountTable():

	CountDf = {};Header = ['Locus','Gene','Region']
	for ColName in Header:CountDf.setdefault(ColName,[])

	BedList = sorted(glob.glob('3.Macs/*_peaks.bed'))
	SampleList = [i.split('/')[-1].replace('_peaks.bed','') for i in BedList]
	for Sample in SampleList:
		Header.append(Sample)
		CountDf.setdefault(Sample,[])

	GeneList = [];RegionList = []
	LociList = [];PrevLocus = '';CountDic = {}
	IntersectFile = '3.Macs/IntersectPeaks_Annotated.bed'
	fin = open(IntersectFile,'r')
	for line in fin:

		Scaf,Start,End,Gene,Strand,Type,Sample = line.strip().split('\t')[:7]
		Locus = '%s:%s-%s' % (Scaf,Start,End)
		Pileup = line.strip().split('\t')[12]
		
		if not PrevLocus == Locus:
			LociList.append(Locus)
			GeneList.append(Gene)
			RegionList.append(Type)
			for IndSample in SampleList:
				CountDic.setdefault((IndSample,Locus),0)

		CountDic[(Sample,Locus)] = Pileup
		PrevLocus = Locus

	fin.close()

	for i in range(len(LociList)):
		CountDf[Header[0]].append(LociList[i])
		CountDf[Header[1]].append(GeneList[i])
		CountDf[Header[2]].append(RegionList[i])
		for Sample in SampleList:
			CountDf[Sample].append(CountDic[(Sample,LociList[i])])

	CountDf = pd.DataFrame(CountDf,columns=Header)
	CountDf.to_csv('3.Macs/CountDf.txt',sep='\t',index=False)


def MakeGeneCoDf(GffFile,GeneDfFile):

	GeneCoDf = {};Header = ['Scaf','Start','End','Gene','Strand','Region']
	for ColName in Header:GeneCoDf.setdefault(ColName,[])

	fin = open(GffFile,'r')
	for line in fin:

		if not '#' in line[0]:

			Type = line.strip().split('\t')[2]
			if Type in ['exon','intron']:

				GeneId = line.strip().split('Name=')[1].split(';')[0]
				Strand = line.strip().split('\t')[6]
				Scaf = line.strip().split('\t')[0]
				Start = int(line.strip().split('\t')[3])
				End = int(line.strip().split('\t')[4])

				Elements = [Scaf,Start,End,GeneId,Strand,Type]
				for i in range(len(Header)):
					GeneCoDf[Header[i]].append(Elements[i])

	fin.close()

	GeneCoDf = pd.DataFrame(GeneCoDf,columns=Header)
	GebeCoDf = GeneCoDf.drop_duplicates()
	GeneCoDf.to_csv(GeneDfFile,sep='\t',index=False,header=False)


def MakeUtrCoDf(GffFile,UtrDfFile):

	TransCoDic = {};CdsEndDic = {};CdsStartDic = {}
	fin = open(GffFile,'r')
	for line in fin:

		if not '#' in line[0]:

			if 'exon' == line.strip().split('\t')[2]:

				GeneId = line.strip().split('Name=')[1].split(';')[0]

				Strand = line.strip().split('\t')[6]
				ExonId = line.strip().split('ID=exon:')[1].split(';')[0]
				TransId = line.strip().split('ID=exon:')[1].split(';')[0].split(':')[0]

				Scaf = line.strip().split('\t')[0]
				Start = int(line.strip().split('\t')[3])
				End = int(line.strip().split('\t')[4])

				TransCoDic.setdefault((GeneId,TransId,Strand),[])
				TransCoDic[(GeneId,TransId,Strand)].append([Scaf,Start,End])

			if 'stop_codon' == line.strip().split('\t')[2]:

				GeneId = line.strip().split('Name=')[1].split(';')[0]
	
				TransId = line.strip().split('ID=stop_codon:')[1].split(';')[0]
				Strand = line.strip().split('\t')[6]

				if Strand == '+':CdsEnd = int(line.strip().split('\t')[4])
				elif Strand == '-':CdsEnd = int(line.strip().split('\t')[3])

				CdsEndDic[(GeneId,TransId)] = CdsEnd

			if 'start_codon' == line.strip().split('\t')[2]:

				GeneId = line.strip().split('Name=')[1].split(';')[0]
				TransId = line.strip().split('ID=start_codon:')[1].split(';')[0]
				Strand = line.strip().split('\t')[6]

				if Strand == '+':CdsStart = int(line.strip().split('\t')[3])
				elif Strand == '-':CdsStart = int(line.strip().split('\t')[4])

				CdsStartDic[(GeneId,TransId)] = CdsStart

	fin.close()

	UtrCoDic = {}
	for GeneId,TransId,Strand in sorted(TransCoDic.keys()):

		UtrCoDic.setdefault((GeneId,TransId,Strand),[])
		if (GeneId,TransId) in CdsEndDic:
			CdsEnd = CdsEndDic[(GeneId,TransId)]

			for Scaf,Start,End in TransCoDic[(GeneId,TransId,Strand)]:

				if Strand == '-':
					if Start < CdsEnd:
						if End < CdsEnd:
							UtrCoDic[(GeneId,TransId,Strand)].append([Scaf,Start,End,'3Utr'])
						elif End >= CdsEnd:
							UtrCoDic[(GeneId,TransId,Strand)].append([Scaf,Start,CdsEnd-1,'3Utr'])
					
				elif Strand == '+':
					if End > CdsEnd:
						if Start > CdsEnd:
							UtrCoDic[(GeneId,TransId,Strand)].append([Scaf,Start,End,'3Utr'])
						elif Start <= CdsEnd:
							UtrCoDic[(GeneId,TransId,Strand)].append([Scaf,CdsEnd+1,End,'3Utr'])

		if (GeneId,TransId) in CdsStartDic:
			CdsStart = CdsStartDic[(GeneId,TransId)]

			for Scaf,Start,End in TransCoDic[(GeneId,TransId,Strand)]:

				if Strand == '-':
					if End > CdsStart:
						if Start >= CdsStart:
							UtrCoDic[(GeneId,TransId,Strand)].append([Scaf,Start,End,'5Utr'])
						elif Start < CdsStart:
							UtrCoDic[(GeneId,TransId,Strand)].append([Scaf,CdsStart+1,End,'5Utr'])
				elif Strand == '+':
					if Start < CdsStart:
						if End <= CdsStart:
							UtrCoDic[(GeneId,TransId,Strand)].append([Scaf,Start,End,'5Utr'])
						elif End > CdsStart:
							UtrCoDic[(GeneId,TransId,Strand)].append([Scaf,Start,CdsStart-1,'5Utr'])

	UtrCoDf = {};Header = ['Scaf','Start','End','Gene','Strand','Region']
	for ColName in Header:UtrCoDf.setdefault(ColName,[])

	for GeneId,TransId,Strand in sorted(UtrCoDic.keys()):
		for Scaf,Start,End,Region in UtrCoDic[(GeneId,TransId,Strand)]:
			Elements = [Scaf,Start,End,GeneId,Strand,Region]
			for i in range(len(Elements)):UtrCoDf[Header[i]].append(Elements[i])

	UtrCoDf = pd.DataFrame(UtrCoDf,columns=Header)
	UtrCoDf = UtrCoDf.drop_duplicates()
	UtrCoDf.to_csv(UtrDfFile,sep='\t',index=False,header=False)


def MakeGenicDf(GffFile,GenicDfFile):

	GenicCoDf = {};Header = ['Scaf','Start','End','Gene','Strand','Region']
	for ColName in Header:GenicCoDf.setdefault(ColName,[])

	RegionLength = 10000

	fin = open(GffFile,'r')
	for line in fin:

		if not '#' in line[0]:

			GeneId = line.strip().split('Name=')[1].split(';')[0]
			Scaf = line.strip().split('\t')[0]
			Type = line.strip().split('\t')[2]
			Start = int(line.strip().split('\t')[3])
			End = int(line.strip().split('\t')[4])
			Strand = Strand = line.strip().split('\t')[6]

			if 'gene' == Type:

				if Strand == '+':

					if Start-1000 < 0:
						Promoter = [Scaf,1,Start,GeneId,Strand,'Promoter']
						Upstream = []
					elif Start-RegionLength-1000 < 0:
						Promoter = [Scaf,Start-1000,Start,GeneId,Strand,'Promoter']
						if Start-1001 < 0:
							Upstream = []
						else:
							Upstream = [Scaf,1,Start-1001,GeneId,Strand,'Upstream']
					else:
						Promoter = [Scaf,Start-1000,Start,GeneId,Strand,'Promoter']
						Upstream = [Scaf,Start-1000-RegionLength,Start-1001,GeneId,Strand,'Upstream']
					Downstream = [Scaf,End,End+RegionLength,GeneId,Strand,'Downstream']

				elif Strand == '-':

					Promoter = [Scaf,End,End+1000,GeneId,Strand,'Promoter']
					Upstream = [Scaf,End+1001,End+1000+RegionLength,GeneId,Strand,'Upstream']
					if Start-RegionLength < 0:
						Downstream = [Scaf,1,Start,GeneId,Strand,'Downstream']
					else:
						Downstream = [Scaf,Start-RegionLength,Start,GeneId,Strand,'Downstream']

				for i in range(len(Header)):

					GenicCoDf[Header[i]].append(Promoter[i])
					if not len(Upstream) == 0:
						GenicCoDf[Header[i]].append(Upstream[i])
					if not len(Downstream) == 0:
						GenicCoDf[Header[i]].append(Downstream[i])

	fin.close()

	GenicCoDf = pd.DataFrame(GenicCoDf,columns=Header)
	GenicCoDf = GenicCoDf.drop_duplicates()
	GenicCoDf.to_csv(GenicDfFile,sep='\t',index=False,header=False)


def RemoveDuplicate():

	BamList = sorted(glob.glob('2.Bwa/*Sorted.bam'))
	for BamFile in BamList:

		OutBam = BamFile.replace('Sorted.bam','_DupRmv.bam')
		OutMat = BamFile.replace('Sorted.bam','_DupMetrics.txt')
		SysStr = 'java -jar %spicard.jar '+ \
				 'MarkDuplicates '+ \
				 'REMOVE_DUPLICATES=true '+ \
				 'I=%s O=%s M=%s'
		os.system(SysStr %(ProgramDir,BamFile,OutBam,OutMat))

		os.system('samtools index %s' % OutBam)


def Filter_MtDna_MapQ():

	BamList = sorted(glob.glob('2.Bwa/*_DupRmv.bam'))
	for BamFile in BamList:

		MapQFiltered = BamFile.replace('_DupRmv.bam','_MapQFilter.bam')
		os.system('samtools view -b -q 10 %s > %s' % (BamFile,MapQFiltered))

		MtRemoved = BamFile.replace('_DupRmv.bam','_MapQFilter_MtRemoved_Ind.bam')
		Mt = BamFile.replace('_DupRmv.bam','_MapQFilter_Mt_Ind.bam')
		SysStr = 'samtools view -b '+ \
				 '-L MtLoci.txt '+ \
				 '-U %s %s > %s'
		os.system(SysStr % (MtRemoved,MapQFiltered,Mt))


def SampleReads_n_Merge():

	PBamList = sorted(glob.glob('2.Bwa/pec-pelvatac_*_MapQFilter_MtRemoved_Ind.bam'))
	SysStr = 'samtools merge 2.Bwa/PecPel_MapQFilter_MtRemoved.bam '+' '.join(PBamList)
	print(SysStr);os.system(SysStr)

	TBamList = sorted(glob.glob('2.Bwa/spinalatac_*_MapQFilter_MtRemoved_Ind.bam'))
	SysStr = 'samtools merge 2.Bwa/Spinal_MapQFilter_MtRemoved.bam '+' '.join(TBamList)
	print(SysStr);os.system(SysStr)


def FlagStats():

	MappedDic = {}
	BamList = sorted(glob.glob('2.Bwa/*_DupRmv.bam'))
	for BamFile in BamList:

		OutFile = BamFile+'.FlagStat'
		SysStr = ProgramDir+'sambamba-0.7.1-linux-static '+ \
				 'flagstat '+ \
				 '-t 6 '+ \
				 BamFile+ \
				 ' > '+OutFile

		if not OutFile.split('/')[-1] in os.listdir('2.Bwa/'):
			print(SysStr);os.system(SysStr)

		fin = open(OutFile,'r')
		for line in fin:

			if 'mapped (' in line:
				MappedReads = int(line.strip().split()[0])
				MappedDic[BamFile] = MappedReads

		fin.close()

	return MappedDic


def InsertSizeDist():

	BamList = sorted(glob.glob('2.Bwa/*_DupRmv.bam'))
	for BamFile in BamList:

		MetricsFile = BamFile.replace('.bam','_Metrics.txt')
		HistFile = BamFile.replace('.bam','_Hist.pdf')

		SysStr = 'java -jar %spicard.jar '+ \
				 'CollectInsertSizeMetrics I=%s O=%s H=%s'
		os.system(SysStr % (ProgramDir,BamFile,MetricsFile,HistFile))


def MakeGffFile(FileName):

	GffFile = '3.Macs/AtacPeaks.gff'

	fout = open(GffFile,'w')
	fin = open(FileName,'r')
	for line in fin:

		Scaf,Start,End = line.strip().split('\t')
		Locus = '%s:%s-%s' % (Scaf,Start,End)
		NewLine = '%s\tCAT\texon\t%s\t%s\t.\t+\t.\tgene_id "%s";\n'
		fout.write(NewLine % (Scaf,Start,End,Locus))

	fin.close();fout.close()


def MergePeaks():

	SampleList = []
	PeakFiles = sorted(glob.glob('3.Macs/*_peaks.xls'))
	for PeakFile in PeakFiles:

		SampleList.append('_'.join(PeakFile.split('/')[-1].split('_')[:3]))
		OutFile = PeakFile.replace('_peaks.xls','_peaks.bed')
		fout = open(OutFile,'w')
		fin = open(PeakFile,'r')
		for line in fin:

			if '#' in line[0]:

				pass

			elif 'chr' in line.split('\t')[0]:

				pass

			else:

				fout.write(line)

		fin.close()
		fout.close()

	AllPeakFile = '3.Macs/AllSamples_AllPeaks.bed'
	os.system('cat 3.Macs/*_peaks.bed > %s' % AllPeakFile)

	SortedFile = '3.Macs/AllSamples_AllPeaks_sorted.bed'
	SysStr = 'sort -k1,1 -k2,2n %s > %s' % (AllPeakFile,SortedFile)
	os.system(SysStr)

	ConsensusPeakFile = '3.Macs/ConsensusPeaks.bed'
	os.system('bedtools merge -i %s > %s' % (SortedFile,ConsensusPeakFile))

	PeakBedFiles = sorted(glob.glob('3.Macs/*_peaks.bed'))
	IntersectFile = '3.Macs/Intersect_perSample_Peaks.bed'
	SysStr = 'bedtools intersect '+ \
			 '-wa -wb '+ \
			 '-a '+ConsensusPeakFile+' '+ \
			 '-b '+' '.join(PeakBedFiles)+' '+ \
			 '-names '+' '.join(SampleList)+' '+ \
			 '> '+IntersectFile
	os.system(SysStr)

	fin = open(IntersectFile,'r');PeakDic = {}
	for line in fin:

		Scaf,Start,End,Sample = line.strip().split('\t')[:4]
		PeakDic.setdefault((Scaf,int(Start),int(End)),[])
		if not Sample in PeakDic[(Scaf,int(Start),int(End))]:
			PeakDic[(Scaf,int(Start),int(End))].append(Sample)

	fin.close()

	FilteredPeaks = '3.Macs/FilteredPeaks.bed'
	fout = open(FilteredPeaks,'w')
	for Scaf,Start,End in sorted(PeakDic.keys(),key=lambda x:(x[0],x[1],x[2])):
		if len(PeakDic[(Scaf,Start,End)]) >= 2:
			fout.write('%s\t%d\t%d\n' % (Scaf,Start,End))
	fout.close()

if __name__ == "__main__":

	Thread = 10
	ProgramDir = 'path/to/Programs/'
	ReferenceFile = '../Reference/LittleSkate_S1_RM.fa'
	GffFile = '../Reference/LEUER_HoxaEdited.gff'

	RawDataList = sorted(glob.glob('AtacSeqData/*_R1_001.fastq.gz'))

	## Original Flow ##

	#Trimmomatic()
	BwaIndex()
	Bwa()
	RemoveDuplicate()
	#InsertSizeDist()

	#Filter_MtDna_MapQ()

	#SampleReads_n_Merge()

	#MACS()
	#MergeBed()
