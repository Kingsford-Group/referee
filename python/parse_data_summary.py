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
			num = int( line.strip().split()[0] )
			if num % 4 == 0:
				data[fname]["unalign"] = num / 4
			else:
				if num % 3 == 0:
					data[fname]["unalign"] = num / 3

# TODO: read these from read_length.txt file
read_lens = {"SRR445718": 100, "SRR1294122": 101, "K562_cytosol_LID8465_TopHat_v2": 76, "P_aeruginosa_PAO1": 51, "NA12878_S1": 101}
types = {"SRR445718": "RNA-seq", "SRR1294122": "RNA-seq", "K562": "RNA-seq", "P_aeruginosa": "RNA-seq", "NA12878_S1": "WGS"}

f_out = "aligned_counts.txt"
f_al = open(f_out, "w")
for fname, vals in data.iteritems():
	type="RNA-seq"
	rlen = read_lens[fname]
	f_al.write( "{} {} {} {}\n".format(fname, vals["align"], vals["error"], vals["depth"] ) )
	print "{} & {} & {} & {} & {} ".format(fname, type, rlen, vals["align"],
		vals["unalign"]) #, vals["error"], vals["depth"])

f_al.close()
print "Wrote alignment counts to " + f_out
