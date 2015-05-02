import sys
import os
import math

def createTableHeader(methods):
	#print "\\begin{table}[ht!]"
	#print "\\caption{XXX}"
	#print "\\label{tab:runtimes}"
	#print "\\centering"
	print "\\begin{tabular}{l c r r r }"
	print "\\toprule"
	print "File & {} \\\\ ".format( " & ".join(methods) )
	print "\\midrule"

def convertToTime(sec):
	if sec == None:
		t1_s = "XXX"
	else:
		h = int(sec) / 3600
		t1_s = "{}h".format(h) if h > 0 else ""
		m = (int(sec) - h * 3600) / 60
		t1_s += "{}m".format(m) if m > 0 else "00"
		sec = int(math.ceil(sec - h * 3600 - m * 60))
		t1_s += "{}s".format(sec)
	return t1_s

def createTableRow(file, t1, t2, t3):
	t1_s = convertToTime(t1)
	t2_s = convertToTime(t2)
	t3_s = convertToTime(t3)
	print "{} & {} & {} & {} \\\\".format(file, t1_s, t2_s, t3_s)

def createTableFooter():
	print "\\bottomrule"
	print "\\end{tabular}"
	print "\\end{table}"

def parseLogForTime(log_path):
	if not os.path.isfile(log_path):
		# print "Can't find log", log_path
		return None
	with open(log_path, "r") as f:
		for line in f:
			if "Elapsed (wall clock)" in line:
				t = line.strip().split()[-1]
				time_parts = t.split(":")
				sec = 0
				for p in time_parts[:-1]:
					sec = sec * 60 + int(p) * 60
				sec += float(time_parts[-1])
				return sec
	return None

FILES=["SRR445718", "SRR1294122", "P_aeruginosa_PAO1", "K562_cytosol_LID8465_TopHat_v2", "NA12878_S1"]
methods=["Referee","Deez","Quip -r"]
DIR = sys.argv[1]
createTableHeader(methods)
for file in FILES:
	ref_log = DIR + "/" + file + ".comp.log"
	deez_log = DIR + "/" + file + ".dz_comp.log"
	quip_log = DIR + "/" + file + ".quip.log"

	ref_time = parseLogForTime(ref_log)
	deez_time = parseLogForTime(deez_log)
	quip_time = parseLogForTime(quip_log)

	createTableRow(file, ref_time, deez_time, quip_time)

createTableFooter()
