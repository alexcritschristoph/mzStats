import argparse
import sys
def main(args):
	#get mapping data
	group_files = []
	#If no min was suggested, min = length(group1)
	f = open(args.mapping)
	i = 0
	for line in f.readlines():
		cols = line.strip().split("\t")
		position = 0
		cols = line.split("\t")
		if cols[2].strip() == args.category:
			print cols[0] + " is in your group."
			group_files.append(cols[0].strip())
	f.close()

	#If min was not set, min = size of group
	if not args.min_count:
		minimum = len(group_files)
	else:
		minimum = int(args.min_count)

	#get compound data
	#Read in compound information
	compound_abundances = {}

	all_compounds = []
	i = 0
	f = open(args.input)
	compound_name_to_mz = {}
	for line in f.readlines():
		if i == 0:
			samples = line.strip().split(",")
			samples.pop(0)
			samples.pop(0)
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

	#For each comound: is it found in > n members of the group and NOT found in any other sample?
	f = open('compound_table_unique_' + args.category + '.txt', 'w+')
	f.write("Compound,Mass")
	for sample in sorted(compound_abundances.keys()):
		f.write("," + sample)
	f.write("\n")

	for compound in all_compounds:
		count_good = 0
		for sample in group_files:
			if compound_abundances[sample][compound] != 0:
				count_good += 1

		found = False
		for sample in compound_abundances.keys():
			if sample not in group_files:
				if compound_abundances[sample][compound] != 0:
					found = True
					break

		if not found and count_good >= minimum:
			f.write(compound + "," + compound_name_to_mz[compound])
			for sample in sorted(compound_abundances.keys()):
				f.write("," + str(compound_abundances[sample][compound]))
			f.write("\n")
	f.close()


if __name__ == '__main__':

	__author__ = "Alex Crits-Christoph"
	parser = argparse.ArgumentParser(description='Filters an MZstats compound table by finding compounds unique to one condition (e.g. one set of replicates) that are not found in other samples in the dataset.')
	parser.add_argument('-i','--input', help='Path to input compound table',required=True)
	parser.add_argument('-m','--mapping', help='Path to input mapping file',required=True)
	parser.add_argument('-c','--category', help="Name of group to find unique compounds for in the 3rd column of the mapping file.", required=True)
	parser.add_argument('-n','--min_count', help="", required=False)


	args = parser.parse_args()

	main(args)