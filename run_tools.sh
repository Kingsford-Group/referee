# ./so_plzip.sh
# export DYLD_LIBRARY_PATH=plzip:$DYLD_LIBRARY_PATH
rm bin/rsupport
make rsupport
if [ ! -e bin/rsupport ]
	then
	echo "Compilation failed"
	exit
fi
time ./bin/rsupport depth /data/referee/aligned/SRR445718.sam 100
