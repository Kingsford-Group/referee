# decompression tests in snakemake

data_dir="/data/referee/aligned"
file_prefix=["SRR445718.2m"]#, "NA12878_S1.1mil"]
# file_prefix="temp"
# file_prefix="P_aeruginosa_PAO1.1mil"
genome = data_dir + "/human_genome.fa"
# genome = data_dir + "/bact.fa"
# genomes = [hum, bacte]
diff_params = "--suppress-common-lines -i -w -y -B"

# seq_streams = ["has_edits", "edits", "left_clip", "offs", "right_clip", "flags"]

# recompile the tool if any of the essential compression/decompression
# sources have changed
rule compile:
	input:
		"include/RefereeCompress.hpp",
		"include/RefereeDecompress.hpp",
		"include/compress/Compressor.hpp",
		"include/compress/OutputBuffer.hpp",
		"include/decompress/Decompressor.hpp",
		"include/decompress/InputBuffer.hpp",
		"include/decompress/ClipStream.hpp",
		"include/decompress/OffsetsStream.hpp"
	output:
		"bin/referee"
	shell:
		'make'

##################################################################
# compress SAM file only taking care of the sequence files
##################################################################
rule compress_seq:
	input:
		"{FILE}.sam"

	output:
		# expand("{FILE}.sam.{TYPE}.lz", TYPE=seq_streams)
		"{FILE}.sam.has_edits.lz", 
		"{FILE}.sam.edits.lz", 
		"{FILE}.sam.left_clip.lz", 
		"{FILE}.sam.right_clip.lz", 
		"{FILE}.sam.offs.lz", 
		"{FILE}.sam.flags.lz"

	shell:
		'''
		bin/referee {input} -t 6 -r {genome}
		'''

##################################################################
# 
##################################################################
rule decompress_seq:
	input:
		"{FILE}.sam.has_edits.lz", 
		"{FILE}.sam.edits.lz", 
		"{FILE}.sam.left_clip.lz", 
		"{FILE}.sam.right_clip.lz", 
		"{FILE}.sam.offs.lz", 
		"{FILE}.sam.flags.lz"

	output:
		"{FILE}.sam.recovered"

	shell:
		# 'bin/referee -d {output} -t 6 -r {genome} > {output}.decomp.log 2>&1'
		'bin/referee -d {output} -t 6 -r {genome}'

##################################################################
# compute diff between the original and recovered sam files
##################################################################
rule compute_diffs:
	input:
		"{FILE}.sam",
		"{FILE}.sam.recovered"

	output:
		"{FILE}.sam.diffs"

	shell:
		# 'diff {diff_params} <(cut -f 10 {data_dir}/{file_prefix}.sam) {data_dir}/{file_prefix}.sam.recovered > {output}'
		'diff {diff_params} <(cut -f 1,3,4 {input[0]}) <(cut -f 1,3,4 {input[1]}) > {output}'

##################################################################
rule final:
	input:
		expand("{DIR}/{FILE}.sam.diffs", DIR=data_dir, FILE=file_prefix)