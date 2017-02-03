import sys
import argparse

def main(args):

	input_file = args.input
	mapping_file = args.mapping
	media_name = args.control

	#Read in compound information
	compound_abundances = {}

	all_compounds = []
	i = 0
	f = open(input_file)
	compound_name_to_mz = {}
	for line in f.readlines():
		if i == 0:
			samples = line.strip().split(",")
			samples.pop(0)
			samples.pop(0)
			samples.pop()
			print len(samples)
			for sample in samples:
				compound_abundances[sample] = {}
		else:
			compounds = line.strip().split(",")
			compound_name = compounds[0]
			compound_name_to_mz[compound_name] = compounds[1]
			all_compounds.append(compound_name)
			j = 0
			for stat in compounds[2:]:
				compound_abundances[samples[j]][compound_name] = float(stat)
				j += 1
		i += 1

	f.close()
	#Get filenames of media

	f = open(mapping_file)
	i = 0
	media_files = []
	for line in f.readlines():
		cols = line.strip().split("\t")
		if i != 0:
			if media_name == cols[2].strip():
				print cols[0] + " is a media sample."
				media_files.append(cols[0].strip())
		i += 1

	f.close()

	f = open(input_file.split(".")[0] + '_filtered.txt', 'w+')
	f.write("Compound,Mass")
	for sample in sorted(compound_abundances.keys()):
		f.write("," + sample)
	f.write("\n")

	#for every peak, if it is not found in media, write the  out again
	found_compounds = []
	for compound in all_compounds:
		found_in_media = False

		if args.control:
			for sample in compound_abundances.keys():
				if sample in media_files:
					if compound_abundances[sample][compound] != 0:
						found_compounds.append(compound)
						found_in_media = True
						break

		#If subtracting controls
		if (not found_in_media) or (not args.control):
			
			#If not singleton
			abunds = 0
			for sample in sorted(compound_abundances.keys()):
				if compound_abundances[sample][compound] != 0:
					abunds += 1
			if (abunds > 1) or (not args.filter_singletons):

				#If minimum mz
				if (not args.min_mz) or (float(compound_name_to_mz[compound].split("+")[0]) >= float(args.min_mz)):
					f.write(compound + "," + compound_name_to_mz[compound])
					for sample in sorted(compound_abundances.keys()):
						f.write("," + str(compound_abundances[sample][compound]))
					f.write("\n")

	f.close()

if __name__ == '__main__':

	__author__ = "Alex Crits-Christoph"
	parser = argparse.ArgumentParser(description='Filters an MZstats compound table by removing compounds found in controls, small compounds, and singletons.')
	parser.add_argument('-i','--input', help='Path to input compound table',required=True)
	parser.add_argument('-m','--mapping', help='Path to input mapping file',required=True)
	parser.add_argument('-c','--control', help="Name of control samples in mapping file's third column - subtracts any compound found in these controls",required=False)
	parser.add_argument('-z','--min_mz', help="Minimum mz to retain in the compound table",required=False)
	parser.add_argument('--filter_singletons', dest='filter_singletons', action='store_true')

	args = parser.parse_args()

	main(args)