import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

files = {}

# output = sys.argv[2]
# expected format:
# file, SRR123
# method, methodXXX
# quals, <size-in-bytes>
# seq, <size-in-bytes>
# ...
print "Reading sizes from", sys.argv[1]
fields = []
with open(sys.argv[1], "r") as f:
	for line in f:
		if "file" in line:
			# print line
			_, filename = line.strip().split(", ")
			if not(filename in files):
				files[filename] = {}
		elif "method" in line:
			_, method = line.strip().split(",")
			method = method.strip()
			files[filename][method] = {}
		else:
			data_type, bytes = line.strip().split(",")
			fields.append(data_type)
			files[filename][method][data_type] = float(bytes.strip() )

# print "Total sizes, bytes"
# methods = files.values()[0].keys()
# methods = list(set(methods))
# for method,v in data.iteritems():
# 	print method + ":", sum(v.values())
method = sys.argv[2]
print "Plotting for {} only".format(method)

width=0.2

fields = {f: [] for f in set(fields) }
colors = ["#1B9E77", "#E6AB02", "#A6761D", "#7570B3", "#D95F02"]

total_sizes = [(f, sum(v[method].values() ) ) for f,v in files.iteritems() ]
total_sizes = sorted(total_sizes, key=lambda x:x[1])
file_names, sizes = zip(*total_sizes)
print file_names
print sizes

totals = []
for file in file_names:
	data = files[file]
	S = 0
	for field,value in data[method].iteritems():
		fields[field].append(value)
		S += value
	totals.append(S)
	print file, 100 - ( data[method]["fields"] + data[method]["seq"] ) * 100.0 / float(S), "%"



X = range(len(file_names))
bottoms = np.zeros(len(file_names))
i = 0
for field, arr in fields.iteritems():
	Y = [arr[j] * 100.0 / totals[j] for j in xrange(len(arr))]
	plt.bar(X, Y, width, color=colors[i%len(colors)], bottom=bottoms, label=field, linewidth=0)
	bottoms += Y
	i+=1

# change fonts to Type 1
# they are type 3 by default
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

plt.ylabel('Size, % of total')
# plt.yscale("log")
# plt.yticks()
X = np.arange(len(files))
# plt.xticks(X + width/2.0, files.keys(), rotation=-20)
plt.xticks(X + width/2.0, file_names)
plt.xlim([-width, max(X) + 0.5])
plt.ylim([0,100])
plt.legend(loc=4)
# plt.legend()
fname = method + "-size-breakdown.pdf"
plt.savefig(fname)
print "Saved stacked bar chart to", fname
