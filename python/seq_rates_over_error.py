import sys
import matplotlib.pyplot as plt
# from matplotlib.font_manager import FontProperties

# data aggregated from logs by: python/parse_compr_rates_error_rate.py error_rates.txt
# error_rates.txt created by: referee_error_rates.sh 

data = {}
edits_data = {}
with open(sys.argv[1], "r") as f:
	i = 0
	for line in f:
		if "#" in line: continue
		if i == 0:
			error_rates = map(float, line.strip().split() )
		elif "edits" in line:
			parts = line.strip().split()
			fname = parts[0]
			edits = map(float, parts[2:])
			edits_data[fname] = edits
		else:
			parts = line.strip().split()
			fname = parts[0]
			comp_rates = map(float, parts[1:])
			data[fname] = comp_rates
		i+=1

print data, edits_data

# font = FontProperties()
# font.set_size('large')

# plot bit per base against error rate cutoff
# for k,v in data.iteritems():
	# plt.plot(v[:-1], "-o", label=k)
# plt.xlabel("max error rate per read")
# plt.xticks(range(len(comp_rates)), error_rates)

# plot bpb against edit count
# for fname, v in edits_data.iteritems():
# 	v = [d_v / float(10**6) for d_v in v]
# 	plt.plot(v, data[fname], "-o", label=fname)
# plt.xlabel("Total edit count, in millions")

fs = 16

# plot bpb against error rate
for fname, v in edits_data.iteritems():	
	plt.plot(v, data[fname], "-o", label=fname, markersize=8)
plt.xlabel("Error rate per alignment", fontsize=fs)

error_rate_at_05_srr445718 = 3.43
error_rate_at_05_srr1224129 = 2.68
plt.plot(error_rate_at_05_srr445718,	0.69, "b*", label="SRR445718 (Deez)", markersize=12)
plt.plot(error_rate_at_05_srr1224129, 	0.46, "g*", label="SRR1294122 (Deez)", markersize=12)
plt.plot(error_rate_at_05_srr445718, 	0.34, "bp", label="SRR445718 (Quip)", markersize=12)
plt.plot(error_rate_at_05_srr1224129, 	0.26, "gp", label="SRR1294122 (Quip)", markersize=12)

plt.ylabel("Bits per base", fontsize=fs)
plt.ylim([0, 0.75])

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
# plt.legend(loc=0)
# plt.xlim([-0.5, len(comp_rates) - 1.5])
plt.savefig("analyses/seq_rates_over_error.pdf")
plt.show()
