import sys
import os

fd_map = {}
fd_full = {}
orig_sizes = {}
orig_sizes_total = 0
orig_cluster_total = 0
comp_total = 0
comp_clusters_total = 0
with open(sys.argv[1], "r") as f:
	for line in f:
		if "Read" in line and "alignments" in line:
			_, alignments, _ = line.strip().split()[-3:]
			alignments = int(alignments)
		elif "Of them unaligned:" in line:
			unaligned = int(line.strip().split()[-1])
		elif "Assuming uniform read len of" in line:
			read_len = int(line.strip().split()[-2])
		elif "Opened " in line and "stream for compressed data" in line:
			parts = line.strip().split()
			fd = int(parts[-1].strip("(fd=").strip(")") )
			fname = parts[1].split("/")[-1]
			fd_map[fd] = fname
			fd_full[fd] = parts[1]
		elif "Closing fd" in line:
			fd_sizes = line.strip().split(";")
			for s in fd_sizes:
				parts = s.strip().split()
				if len(parts) == 1: continue
				fd = int(parts[1].strip(",").strip("Closing fd=") )
				bytes = int(parts[3])
				if os.path.isfile(fd_full[fd]):
					comp_bytes = os.path.getsize(fd_full[fd])
				else:
					comp_bytes = 0
				if comp_bytes != 0:
					ratio = float(bytes) / comp_bytes
				else:
					ratio = 1
				if "k=" in fd_full[fd] : # quality values
					orig_sizes_total += bytes
					comp_total += comp_bytes
					if not ("other" in fd_full[fd]):
						orig_cluster_total += bytes
						comp_clusters_total += comp_bytes				
					else:
						other_bytes = bytes
					print fd_map[fd], comp_bytes, bytes, ratio
# add membership vector
full_path = fd_full[12]
parts = full_path.split("/")
name_parts = parts[-1].split(".")
member_path = "/".join(parts[:-1]) + "/" + ".".join( name_parts[0:3] ) + ".membership"
print member_path
estimated_orig_bytes = 0
#if os.path.isfile(member_path):
#	estimated_orig_bytes = os.path.getsize(member_path)
#else:
#	print "Membership not included"
if os.path.isfile(member_path + ".lz"):
	comp_total += os.path.getsize(member_path + ".lz")
else:
	print "Compressed membership not included"

estimated_orig_bytes += (alignments - unaligned) * read_len
print "Original bytes:", estimated_orig_bytes
print "Total comp size:", comp_total
print "Estimated compression:", estimated_orig_bytes / float(comp_total)
print "Average compression for clusters:", (estimated_orig_bytes - other_bytes) / float(comp_clusters_total)
print "Or: ", float(comp_clusters_total) / (estimated_orig_bytes - other_bytes)
