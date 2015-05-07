import sys
import os

################################
def getReadLengths(file_name):
	rl = {}
	with open(file_name, "r") as f:
		cnt = 0
		for line in f:
			part = line.strip()
			if cnt % 3 == 0: 
				name = part
				rl[name] = []
			elif cnt % 3 == 1: 
				rl[name].append( int(part) )
			elif cnt % 3 == 2: 
				rl[name].append( int(part) )
			cnt+=1
	return rl

################################
# in bytes
def getSize(name):
	if os.path.isfile(name):
		stats = os.stat(name)
		return stats.st_size
	else:
		print "[INFO] File not found: ", name
		return 0

################################
def getRefereeSize(file_name, bases):
	ref_size = 0
	ref_size += getSize(file_name + ".edits.lz")
	ref_size += getSize(file_name + ".has_edits.lz")
	ref_size += getSize(file_name + ".offs.lz")
	ref_size += getSize(file_name + ".left_clip.lz")
	ref_size += getSize(file_name + ".right_clip.lz")
	return "{} ({})".format( round(ref_size / 1024.0 / 1024, 2) , round(ref_size*8.0/bases, 2) )

################################
def parseDeezLog(log_name, read_len):
	if not os.path.isfile(log_name): # file does not exist
		print "No deez file for " + log_name
		return "---"

	# print "parsing", log_name
	with open(log_name, "r") as f:
		for line in f:
			if "Written" in line and "lines" in line:
				line_count = int( "".join( line.strip().split()[1].split(",") ) )
				
			elif "FIX" in line and "REP" in line and "seq" in line:
				parts = line.strip().split()
				fix = int(parts[1])
				rep = int(parts[3].strip("]seq:"))
				seq = int(parts[-1])
				assert fix + rep == seq
			elif "NUC" in line and "UNK"  in line:
				parts = line.strip().split()
				nuc = int(parts[1])
				unk = int(parts[3])
				op = int(parts[5])
				len = int(parts[7])
				loc = int(parts[9])
				sti = int(parts[11].strip("]edits:")) # sexually transmitted infections?..
				edits = int(parts[-1])
				assert nuc + unk + op + len + loc + sti == edits
			elif "IDX" in line and "readIDs" in line:
				parts = line.strip().split()
				idx = int(parts[1])
				s = int(parts[3].strip("]readIDs:"))
				readIDs = int(parts[-1])
				assert idx + s == readIDs
			elif "mapFlag" in line:
				mapFlag = int(line.strip().split()[-1])
			elif "mapQual" in line:
				mapQual = int(line.strip().split()[-1])
			elif "quals" in line:
				quals = int(line.strip().split()[-1])
			elif "pairedEnd" in line:
				pairedEnd = int(line.strip().split()[-1])
			elif "optfield" in line:
				optfield = int(line.strip().split()[-1])
	# print "Seq only:", seq + edits, "bytes", "bpb: ", (seq + edits) * 8/ (line_count * 100.0)

	return "{} ({})".format( round( (seq + edits)/1024.0 / 1024, 2), round( (seq+edits)*8.0 / (line_count * read_len),2 ) )
	# print "Total bases:", line_count, "* read_len"

################################
#
################################
def parseQuipLog(log_name):
	if not os.path.isfile(log_name): # file does not exist
		print "No quip file for " + log_name
		return "---"

	# print "parsing", log_name
	original_bases = 0
	compressed_bases = 0
	with open(log_name, "r") as f:
		for line in f:
			if "seq:" in line:
				parts = line.strip().split()
				compressed_bases += int(parts[1])
				original_bases += int(parts[3])
	# print "Seq only:", compressed_bases, "bytes", "bpb:", compressed_bases*8/float(original_bases)
	return "{} ({})".format( round(compressed_bases/1024.0/1024, 2), round( compressed_bases*8.0/original_bases, 2) )
	# print "Total bases:", original_bases, "bytes"

################################
# Main
################################
RLen = getReadLengths(sys.argv[1])

all_types = ["Referee", "Deez", "Quip"]
types_align = " ".join(["r" for i in range(len(all_types))])
print "\\begin{tabular}{l c r " + types_align + " r r}"
print "\\toprule"
print "File & $|r|$ & Total bases & ", " & ".join(all_types), " & Error rate & Depth \\\\"
print "\midrule"
for file_name in sys.argv[2:]:
	file_suff = file_name.split("/")[-1].split(".")[0]
	if file_suff in RLen:
		rl, aligned = RLen[file_suff]
	else:
		print "read info not available for " + file_suff
		continue
	super_name = file_name.split(file_suff)[0] + file_suff + ".sam"
	# print super_name
	referee = getRefereeSize(super_name, rl * aligned)
	deez = parseDeezLog(file_name + "stripped.deez.log", rl)
	quip = parseQuipLog(file_name + "stripped.quip.log")
	error_rate = "---"
	depth = "---"
	print "{} & {} & {} & {} & {} & {} & {} & {}".format(file_suff, rl, rl * aligned, referee, deez, quip, error_rate, depth)






