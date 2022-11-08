import os, glob

Ncores = 6
Bin = 10
#SmBinList = [100,200,400,1000,2000,10000]
#SmBinList = [10000,200,400,100]
SmBinList = [200]
BamFiles = sorted(glob.glob('2.Bwa/BracLumb_Merged.bam'))

for SmBin in SmBinList:

	for BamFile in BamFiles:

		BwOut = BamFile.replace('.bam','')+'_BinSize%d_SmoothLen%d.bw' % (Bin,SmBin)
		Cmd = 'bamCoverage '+ \
			  '-b %s '+ \
			  '--numberOfProcessors %d '+ \
			  '--effectiveGenomeSize 2652783500 '+ \
			  '--binSize %d '+ \
			  '--normalizeUsing CPM '+ \
			  '-o %s '+ \
			  '--smoothLength %d '+ \
			  '--outFileFormat bigwig'
		print(Cmd % (BamFile,Ncores,Bin,BwOut,SmBin))
		os.system(Cmd % (BamFile,Ncores,Bin,BwOut,SmBin))
