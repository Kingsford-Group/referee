import sys

inv_type = {
        "seq": ["offs", "left_clip", "right_clip", "edits", "has_edits"],
        "quals": ["quals"],
        "unaligned": ["unaligned"],
        "fields": ["opt", "flags"],
        "ids": ["ids"]
        }

data_type = {}
for k,v in inv_type.iteritems():
        for item in v:
                data_type[item] = k

#
#
#
def aggregate(items, method):
	f_out = open("all_files.sizes.aggregated", "w")
	for fname, data in items:
		f_out.write( "file, {}\n".format(fname) )
		f_out.write( "method, {}\n".format(method) )
		agg_s = {k: 0 for k in inv_type.keys() }
		for fsuffix, size in data[method].iteritems():
			t = fsuffix.split(".")[0]
			if not (t in data_type):
				mapped_type = "quals"
			else:
				mapped_type = data_type[t]
			agg_s[ mapped_type ] += size
		for t, s in agg_s.iteritems():
			f_out.write( "{},\t{}\n".format(t, s) )
	f_out.close()

#
#
#
data = {}

# file_type = sys.argv[1].split("/")[-1].split(".")[1]
with open(sys.argv[1], "r") as f:
	for line in f:
		if '#' == line[0]:
			file_type = line.strip().split()[1]
		else:
			parts = line.strip().split()
			file_size = int(parts[0])
			fname_parts = parts[-1].split("/")[-1].split(".")
			fprefix = fname_parts[0]
			if len(fname_parts) > 3:
				fsuffix = ".".join( fname_parts[2:-1] )
			else:
				fsuffix = fname_parts[0]
			if not(fprefix in data):
				data[fprefix] = {}
			if not(file_type in data[fprefix]):
				data[fprefix][file_type] = {}
			data[fprefix][file_type][fsuffix] = file_size

all_types = ["Referee", "Quip", "Deez", "BAM", "SAM"]
#for k,v in data.iteritems():
#	all_types += v.keys()
#all_types = list(set(all_types))
# print all_types
types_align = " ".join(["r" for i in range(len(all_types))])
print "\\begin{tabular}{l " + types_align + "}"
print "\\toprule"
print "File\t&", " & ".join(all_types), "\\\\"
print "\midrule"

items = sorted(data.items(), key=lambda x: sum( x[1]["SAM"].values() ) )

for fname, sizes in items:
	if "P_aer" in fname: fname = "P. aeruginosa"
	elif "SRR129" in fname: 
		fname = fname.split(".")[0]
	elif "SRR4457" in fname: fname = fname.split(".")[0]
	#elif "K562" in fname: fname = "Human RNA-seq II"
	elif "K562" in fname: fname = fname.split(".")[0].split("_")[0]
	# elif "SRR445718" in fname: fname = "Human RNA-seq III"
	#elif "NA128" in fname: fname = "Human genomic"
	elif "NA128" in fname: fname = fname.split(".")[0]
	
	if "_" in fname:
		fname = "\_".join( fname.split("_") )

	# print sizes

	ordered_sizes = []
	for file_type in all_types:
		if file_type in sizes:
			size = round( sum(sizes[file_type].values() ) / 1024.0 / 1024, 2)
		else:
			size = "xxx"
		ordered_sizes.append(size)
	
	# calculate % to the next best
	min_size = min(ordered_sizes)
	next_best = min_size * 10
	for i in xrange(len(ordered_sizes)):
		if ordered_sizes[i] > min_size and ordered_sizes[i] < next_best:
			next_best = ordered_sizes[i]
	if next_best == 0:
		percent = "xxx"
	else:
		percent = round( (1 - min_size / next_best) * 100, 1 )

	sizes_str = []
	for i in xrange(len(ordered_sizes)):
		s = ordered_sizes[i]
		if s == min_size:
			sizes_str.append( "\\textbf{" + str(s) + "} (" + str(percent) + "%)" )
		else:
			sizes_str.append( str(s) )
	# sizes_str = ["\\textbf{"+str(s)+"} (" +  if s == min_size else str(s) for s in ordered_sizes]
	print fname + "\t& " + " & ".join(sizes_str) + " \\\\"

print "\\bottomrule"

# aggregate data by type for referee
aggregate(items, "Referee")
