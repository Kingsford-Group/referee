import sys
import os

# /mnt/scratch0/dfilippo/aligned/SRR445718_all.z\=0.5.sam.referee.seq2.log
def readLogForReadCount(fname, z):
	dir = "/mnt/scratch0/dfilippo/aligned/"
	full_name = dir + fname + "_all.z="
	edits = 1
	if z == 1:
	 	full_name += str(1) + ".sam.referee.seq2.log"
	else:
		full_name += str(z) + ".sam.referee.seq2.log"
	if not os.path.isfile(full_name):
		print "Could not find file", full_name
		return 1, edits

	print "Parsing", full_name
	f = open(full_name, "r")
	if not f:
		print "Could not open file"
		return 1, edits

	for line in f:
		if "Unique total reads" in line:
			n = int(line.strip().split()[-1])
		if "Total edit count:" in line:
			edits = int(line.strip().split()[-1])
			# print edits
	return n, edits

D = {}
zs = []
with open(sys.argv[1], "r") as f:
	for line in f:
		if "***" in line:
			_, fname, z = line.strip().split()
			# print z
			z = float(z.strip("z="))
			# print z
			zs.append(z)
			if not(fname in D):
				D[fname]={}
			if not(z in D[fname]):
				D[fname][z] = {}
				D[fname][z]["bytes"] = 0
	
			# get a number of aligned reads from the log
			reads, edits = readLogForReadCount(fname, z)
			D[fname][z]["N"] = reads
			# D[fname][z]["edits"] = edits 
			# print D[fname][z]

		# elif "offs.lz" in line:
		elif "edits:" in line:
			_, e = line.strip().split()
			e = int(e)
			D[fname][z]["edits"] = e
		elif ".lz" in line:
			size, f = line.strip().split()
			D[fname][z]["bytes"] += int(size)
		elif not("N" in D[fname]): # number of alignment w/in +- 8
			num_alignm, _ = line.strip().split()
			num_alignm = int(num_alignm)
			D[fname][z]["N"] = num_alignm * 8 + 2

readlen = 100

# print zs
zs = sorted(list(set(zs)) )
print " ".join( map(str, zs) )
for f, d in D.iteritems():
	# print f
	btbs = {}
	btbs_edits = {}
	for z, data in d.iteritems():
		# print z, data
		ratio = data["bytes"] * 8.0 / (data["N"] * readlen)
		btbs[z] =  ratio
		btbs_edits[z] = data["edits"] / float(data["N"])
	# print d
	# computes bit per byte	
	# print f, " ".join( map(str, [ d[ zs[i] ] for i in xrange(len(zs)) ]) )
	print f, " ".join( map(str, [ round( btbs[zs[i]], 2) if zs[i] in btbs else "--" for i in xrange(len(zs)) ]) )
	print f, "edits:", " ".join( map(str, [ round( btbs_edits[zs[i]], 2) if zs[i] in btbs_edits else "--" for i in xrange( len(zs) ) ] ) )
	print f, "seqsize:", " ".join( map(str, [ round(d[zs[i]]["bytes"] /1024.0/1024, 2) for i in xrange(len(zs))] ) )
