import sys
import os
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
                        if "qual:" in line:
                                parts = line.strip().split()
                                compressed_bases += int(parts[1])
                                original_bases += int(parts[3])
        # print "Seq only:", compressed_bases, "bytes", "bpb:", compressed_bases*8/float(original_bases)
	return (compressed_bases, original_bases)
        # return "{} ({})".format( round(compressed_bases/1024.0/1024, 2), round( compressed_bases*8.0/original_bases, 2) )
        # print "Total bases:", original_bases, "bytes"

file_name = sys.argv[1]
if os.path.isfile(file_name):
        comp, orig = parseQuipLog(file_name)
	print round(comp / 1024.0 / 1024, 2)
