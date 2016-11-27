from pyteomics import mzxml
from  sklearn.metrics.pairwise import cosine_similarity
import math

scans = []
r = mzxml.read('LG60FAVA_2_2_071816.mzXML')
while True:
	try:
		scans.append(r.next())
	except:
		print r
		break

print "done!"
print len(scans)
base_peaks = {}
peaks = []
all_peaks = []
for scan in scans:
	if scan['msLevel'] == '2':
		base_mz = scan['precursorMz'][0]['precursorMz']
		intensities = scan['intensity array']
		mzs = scan['m/z array'] 
		num = scan['num']
		base_peaks[num] = {"base_mz":base_mz, "intensities":intensities, "mzs":mzs}
		peaks.append(base_mz)
		all_peaks = all_peaks + intensities

#Cluster peaks
print len(peaks)

#Create cosine matrix, intervals of 0.1;
peak_min = math.floor(all_peaks)
peak_max = math.ceil(all_peaks)

#Calculate distance matrix
similarities = []
for peak in base_peaks:
	for peak2 in base_peaks:
		if peak != peak2:
			 sim = 	cosine_similarity(base_peaks[peak]['mzs'], base_peaks[peak2]['mzs'])
			 similarities.append(sim)
#Cluster based on 99% similarity

#