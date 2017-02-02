import argparse
from pyteomics import mzxml
import numpy as np

def main(sample):
	print "Reading " + sample
	
	scans = []
	r = mzxml.read(sample)
	while True:
		try:
			scans.append(r.next())
		except:
			break

	print str(len(scans)) + " scans found in " + sample

	min1 = 10000000000
	max1 = 0
	
	scans = {}

	for scan in scans:
		x = []
		y = []
		if scan['msLevel'] == '1':
			i = 0
			for point in scan['m/z array']:
				x.append(point)
				y.append(scan['intensity array'][i])
				i += 1
			scans[scan['num']] = [x,y]


	while True:
		input("Select a scan number from this file")
		
if __name__ == '__main__':
	__author__ = "Alex Crits-Christoph"
	parser = argparse.ArgumentParser(description='Plots 2D MS2 spectra and 3D plots of entire MS-MS runs.')
	parser.add_argument('input_file')
	args = parser.parse_args()

	main(args.input_file)
