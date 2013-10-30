import re, sys, gc, argparse
from numpy import array

#Author: Nima C. Emami, UCSF BMI (Bioinformatics)
#License: github.com/nemami/PeakTracker/blob/master/LICENSE

class Genome:
	def __init__(self):
	# overestimate chromosome lengths for initializing all bins along chromosome
		self.all_chromsomes = {"chr1":250000000,"chr2":250000000,
							"chr3":200000000,"chr4":200000000,
							"chr5":185000000,"chr6":175000000,
							"chr7":165000000,"chr8":150000000,
							"chr9":145000000,"chr10":140000000,
							"chr11":140000000,"chr12":135000000,
							"chr13":120000000,"chr14":110000000,
							"chr15":105000000,"chr16":95000000,
							"chr17":85000000,"chr18":80000000,
							"chr19":60000000,"chr20":65000000,
							"chr21":50000000,"chr22":55000000,
							"chrX":160000000,"chrY":60000000}
		self.chrom, self.chromMap = "", {}
	def set_chrom(self, inputChrom, windowLength):
		self.chrom = inputChrom
		upperBoundLength = self.all_chromsomes[inputChrom]
		for i in range(upperBoundLength/windowLength):
			self.chromMap[i] = []
	def placePeak(self, new_peak, windowLength):
		for index in new_peak.get_chromIndex(windowLength):
			try:
				self.chromMap[index].append(new_peak)
			except:
				print 'Error: for ',new_peak.get_chrom(), \
					'peak indices out of bounds: ', \
					str(new_peak.get_chromStart()),str(new_peak.get_chromEnd(), \
					". Continuing.")
				continue
	def get_bins(self,count,signal,stdev):
		return [(key, peaks) for (key, peaks) in self.chromMap.iteritems() \
				if len(peaks) >= count and signalMean(peaks) > signal \
				and stdev*float(signalMean(peaks)) >= float(signalStdev(peaks))]

class track:
	def __init__(self):
		self.peaks = []
	def setPeaks(self,trackFile):
		for peakTuple in populatePeaks(trackFile):
			newPeak = peak()
			newPeak.setFields(peakTuple)
			self.peaks.append(newPeak)
	def getPeaks(self):
		return self.peaks

class peak:
	def __init__(self):
		# chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue
		self.chrom, self.chromStart, self.chromEnd, self.name, self.score, \
		self.strand, self.signalValue, self.pValue, self.qValue = "",0,0,"",0,"",0,0,0
	def setFields(self, fields):
		self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand, \
		self.signalValue, self.pValue, self.qValue = fields
	def get_chrom(self):
		return self.chrom
	def get_chromStart(self):
		return self.chromStart
	def get_chromEnd(self):
		return self.chromEnd
	def get_chromIndex(self, windowLength):
		floor = self.chromStart / windowLength
		ceiling = self.chromEnd / windowLength
		return range(floor,ceiling+1)
	def get_signalValue(self):
		return self.signalValue
	def inRange(self, locus):
		return locus >= self.chromStart and locus < self.chromEnd

def signalMean(peaks):
	signalValues = [peak.get_signalValue() for peak in peaks]
	numValues = array(signalValues)
	return float(numValues.mean(axis=0))

def signalStdev(peaks):
	signalValues = [peak.get_signalValue() for peak in peaks]
	numValues = array(signalValues)
	return float(numValues.std(axis=0))

def valueMean(values):
	numValues = array(values)
	return numValues.mean(axis=0)

def valueStdev(values):
	numValues = array(values)
	return numValues.std(axis=0)

def binBoundaries(peaks):
	firstPeak = peaks[0]
	min, max = firstPeak.get_chromStart(), firstPeak.get_chromEnd()
	for i in range(1,len(peaks)):
		peak = peaks[i]
		if peak.get_chromStart() < min:
			min = peak.get_chromStart()
		if peak.get_chromEnd() > max:
			max = peak.get_chromEnd()
	return min, max

def populatePeaks(inFile):
	allPeakLines = []
	regex = '^(chr[^\t]+)\t([^\t]+)+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t\n]+)'
	with open(inFile,"rU") as inFh:
		for line in inFh:
			match = re.search(regex,line)
			if match:
				allPeakLines.append( (match.group(1),int(match.group(2)),int(match.group(3)), \
					match.group(4),int(match.group(5)),match.group(6),float(match.group(7)), \
					float(match.group(8)),float(match.group(9))) )
	return allPeakLines

def main(inputFiles, outputFile, windowLength, inputCount, inputSignal, inputStdev):
	allTracks, totalBins, stddevPercent = [], 0, inputStdev/float(100)

	for inputFile in inputFiles:
		newTrack = track()
		newTrack.setPeaks(inputFile)
		allTracks.append(newTrack)
		print '* Processed input file: ',inputFile

	with open(outputFile,"w") as outFh:
		outFh.write("Chromosome\tWindow_Start\tWindow_End\tSignal_Mean\tSignal_StdDev\n")
	chromList = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',\
				 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',\
				 'chr20','chr21','chr22','chrX','chrY']

	for currentChrom in chromList:
		hg19 = Genome()
		hg19.set_chrom(currentChrom, windowLength)
		for currentTrack in allTracks:
			for peak_k in currentTrack.getPeaks():
				if peak_k.get_chrom() != currentChrom:
					continue
				hg19.placePeak(peak_k,windowLength)
		bins = hg19.get_bins(inputCount, inputSignal, stddevPercent)
		output, lastLocus = [ ], -1
		for locus, peaks in bins:
			binStart, binEnd = binBoundaries(peaks)
			signalValues = [peak.get_signalValue() for peak in peaks]
			if locus != lastLocus + 1:
				output.append(([locus,],signalValues))
			else:
				lastBinIndices, lastBinSignals = output.pop(-1)
				mergeSignals = lastBinSignals
				if len(signalValues) > len(lastBinSignals):
					mergeSignals = signalValues
				output.append((lastBinIndices+[locus,],mergeSignals))
			lastLocus = locus
		print "\t"+currentChrom+": found "+str(len(output))+" windows."
		totalBins += len(output)
		with open(outputFile,"a") as outFh:
			for outIndices, signals in output:
				binStart, binEnd = windowLength*outIndices[0], windowLength*outIndices[-1]
				outFh.write(currentChrom+"\t"+str(binStart)+"\t"+str(binEnd+windowLength)+"\t"+ \
							str(valueMean(signals))+"\t"+str(valueStdev(signals))+"\n")
		gc.collect()
	print 'For input criteria: \n\t', \
			'at least ',str(inputCount),' coinciding peaks per window,\n\t', \
			'at least ',str(inputSignal),' average signal value per window,\n\t', \
			'at most  ',str(inputStdev),'% signal value standard deviation, relative to average,\n\t', \
			'with window length of ',str(windowLength),',\n\t', \
			'returned ',str(totalBins),' total bins. Please see ',outputFile,' for further details.'

def errors(inputFiles, outputFile, windowLength, inputCount, inputSignal, inputStdev):
	if len(inputFiles) < 2:
		print "Error: user must input multiple .narrowpeak and/or .broadpeak files. Exiting."
		return True
	elif inputCount > inputFiles:
		print "Error: -count parameter must not exceed the number of input files. Exiting."
		return True
	return False

if __name__ == "__main__":
	parser = argparse.ArgumentParser( \
	description='''PeakTracker -- returns windows, or ranges of genomic loci, where 
						input [broad|narrow]peak files have overlap and exhibit enrichment.''')
	parser.add_argument('-input', help='Paths to input broadpeak and narrowpeak files.', \
							nargs='+', type = str, required=True)
	parser.add_argument('-output', help='''Name/path of output tab-delimited file of overlapping 
							windows.''', type = str, default="PeakTracker_output.txt")
	parser.add_argument('-window', type = int, help='''Length of windows to which peak signal values 
							will be mapped (Default = 250 bp).''', default=250)
	parser.add_argument('-count', help='''Threshold for number of peaks that must coincide in a single 
							returned bin (e.g. 2, "5") *Note: must be less than or equal 
							to the number of input files.''', type = int, required=True)
	parser.add_argument('-signal', help='''Threshold for average signal value in a returned bin 
							(e.g. 10, 95.6).''', type = float, required=True)
	parser.add_argument('-stdev', help='''Threshold for maximum standard deviation of values for 
							overlapping peaks in a given window, as a percent of average signal value 
							(-signal flag). Example: if invoking -signal 60, with -stdev 20, only peaks
							which (a) contain an average signal of at least 60, AND (b) contain a
							standard deviation of at most 15, will be returned. (Default = 25 %%).''', 
							type = float, default = 25.0)
	args = parser.parse_args()

	if errors(args.input, args.output, args.window, args.count, args.signal, args.stdev):
		sys.exit(1)

	main(args.input, args.output, args.window, args.count, args.signal, args.stdev)