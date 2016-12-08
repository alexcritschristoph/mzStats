from pyteomics import mzxml
from  sklearn.metrics.pairwise import cosine_similarity
import math
import numpy as np
import sys
from scipy.sparse import coo_matrix
from sys import getsizeof

def preprocess_sample(sample):
	print "Reading " + sample
	scans = []
	r = mzxml.read(sample)
	while True:
		try:
			scans.append(r.next())
		except:
			break

	print str(len(scans)) + " scans found in " + sample
	base_peaks = {}
	all_peaks = []
	for scan in scans:
		if scan['msLevel'] == '2':
			base_mz = scan['precursorMz'][0]['precursorMz']
			precursor_intensity = scan['precursorMz'][0]['precursorIntensity']
			intensities = scan['intensity array']
			mzs = scan['m/z array'] 
			num = scan['num']
			base_peaks[num] = {"base_mz":base_mz, "intensities":intensities, "mzs":mzs, "precursor_intensity": precursor_intensity}
			all_peaks = all_peaks + mzs.tolist()

	#Cluster peaks

	#Create cosine vector, intervals of 0.1;
	peak_min = math.floor(min(all_peaks))
	peak_max = math.ceil(max(all_peaks))
	
	#Get rid of really big variables...
	all_peaks = None
	scans = None
	r = None
	
	return peak_min, peak_max, base_peaks

def vectorize_peak(peak_min, peak_max, sample_data):	
	# peak_positions = np.arange(peak_min,peak_max,0.1, dtype=np.float32).tolist()
	vector_length = ( peak_max - peak_min ) * 10
	peak_vectors = {}
	print "Creating peak vectors..."
	for peak in sample_data:
		peak_vector = np.zeros((1,vector_length), dtype=np.float32)
	 	i = 0
		for p in sample_data[peak]['mzs']:
			pos = int((math.floor(p*10)/10 - peak_min) * 10)
			peak_vector[0,pos] = sample_data[peak]['intensities'][i]
			i += 1

	  	peak_vectors[peak] = peak_vector

	print len(peak_vectors)

	print "Finding unique peaks in sample..."
	# #Remove non-unique peaks; peaks that are most identical are grouped and the most intense peak from each group is kept.
	#Only compare peaks that have masses within ~5 DA of each other?
	similarities = []
	already_calculated = []
	peak_vectors_unique = []
	sims = []
	j = 0 
	for peak in peak_vectors:
		found = False
		i = 0
		j += 1
		for peak2 in peak_vectors_unique:
			if abs(sample_data[peak]['base_mz'] - sample_data[peak2]['base_mz']) <= 2:
				sim = cosine_similarity(peak_vectors[peak], peak_vectors[peak2])
				sims.append(sim[0][0])
				if sim[0][0] >= 0.97:	
					#Similar; take the peak with the higher intensity.
					if sample_data[peak]['precursor_intensity'] > sample_data[peak2]['precursor_intensity']:
						peak_vectors_unique[i] = peak			
					found = True
					break
			elif sample_data[peak]['base_mz'] - sample_data[peak2]['base_mz'] < -2:
				break

			i += 1
		#Not similar to already known peaks; add to uniques
		if not found:
			peak_vectors_unique.append(peak)
	
	#Create final data for this sample, return
	final_peaks = {}
	for peak in peak_vectors_unique:
		peak_data = sample_data[peak]
		peak_data['vector'] = coo_matrix(peak_vectors[peak]) #Store only as a COO matrix
		final_peaks[peak] = peak_data

	peak_vectors = None
	#Save peaks to disk in large X (samples) by n array; 
	return final_peaks


def compare_samples(samples_data):

	#Quit if only 1 sample...
	if len(samples_data) <= 1:
		print "There was only one sample processed! Please add more samples to your sample sheet..."
		sys.exit(1)
	#Cluster and count compounds across all samples
	compounds = []

	for sample in samples_data:
		sample_data = samples_data[sample]
		for peak in sample_data:
			sample_data[peak]['sample'] = sample
			found = False
			#Compare to already grouped compounds...
			i = 0
			for compound_group in compounds:
				mass_diff = sample_data[peak]['base_mz'] - compound_group[0]['base_mz']
				if mass_diff <= 2 and mass_diff >= -2:
					sim = cosine_similarity(sample_data[peak]['vector'].toarray(), compound_group[0]['vector'].toarray())
					if sim[0][0] >= 0.97: #Should this be slightly less stringent being between samples?
						compounds[i].append(sample_data[peak])
						found = True
				i += 1

			if not found:
				compounds.append([sample_data[peak]])

	print len(compounds)
	#Write data table to CSV

	line = "Compound,"
	for sample in samples_data:
		line = line + sample + ","
	line.rstrip(",")
	
	f = open('./compound_table.txt', 'a+')
	f.write(line + "\n")

	for compound_group in compounds:

		masses = []
		for compound in compound_group:
			masses.append(compound['base_mz'])

		line = str(np.mean(masses)) + "+-" + str(round(np.std(masses, ddof=0),4)) + ","
		for sample in samples_data:
			found = False
			for compound in compound_group:
				if sample == compound['sample']:
					found = True
					break
			if found:
				line = line + str(compound['precursor_intensity']) + ","
			else:
				line = line + "0,"
		line = line.rstrip(",")
		f.write(line + "\n")

	f.close()

def main():
	#Read in samples
	mapping_file = sys.argv[1]
	#get samples
	samples = []
	f = open(mapping_file)
	for line in f.readlines():
		samples.append(line.split("\t")[0])
	f.close()

	#Preprocess samples
	peak_min = 999999999
	peak_max = 0
	peak_data = {}
	for sample in samples:
		new_min, new_max, new_peak_data = preprocess_sample(sample)
		if new_min < peak_min:
			peak_min = new_min
		if new_max > peak_max:
			peak_max = new_max
		peak_data[sample] = new_peak_data

	print "***"
	print getsizeof(peak_data)
	#Calculate vectors for all samples
	for sample in peak_data:
		peak_data[sample] = vectorize_peak(peak_min, peak_max, peak_data[sample])

	print getsizeof(peak_data)
	print "Comparing samples..."
	
	#Compare samples
	compare_samples(peak_data)

if __name__ == '__main__':
	main()