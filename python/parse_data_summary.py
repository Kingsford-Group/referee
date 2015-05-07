import sys


data = {}
with open(sys.argv[1], "r") as f:
	for line in f:
		if "=== " in line:
			fname = line.strip().split()[-1]
			data[fname] = {"align": 0, "unalign": "XXX", "error": 0, "depth": "XXX"}
		elif "Saw" in line:
			# Saw 35677496 total alignments; of them 
			data[fname]["align"] = int(line.strip().split()[1])
		elif "Error rate: " in line:
			# Error rate: 2.85284%
			data[fname]["error"] = round( float( line.strip().split()[-1].strip("%") ), 2)
		elif "/mnt/scratch0/dfilippo/align" in line:
			# 5577669 /mnt/scratch0/dfilippo/align
			num = int( line.strip().split()[0] ) / 4.0
			if num % 4 == 0:
				data[fname]["unalign"] = num
			else:
				data[fname]["unalign"] = "XXX"

read_lens = {"SRR445718": 100, "SRR1294122": 101, "K562": 76, "P_aeruginosa": 51, "NA12878_S1": 101}
types = {"SRR445718": "RNA-seq", "SRR1294122": "RNA-seq", "K562": "RNA-seq", "P_aeruginosa": "RNA-seq", "NA12878_S1": "WGS"}

for fname, vals in data.iteritems():
	type="RNA-seq"
	rlen = read_lens[fname]
	print "{} & {} & {} & {} & {} ".format(fname, type, rlen, vals["align"],
		vals["unalign"]) #, vals["error"], vals["depth"])

