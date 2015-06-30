# compile plzip as a shared library
# PLZIPDIR=$HOME/plzip_so
PLZIPDIR=plzip

ORIG=.
# compile into position-independent code
cd $PLZIPDIR

INCLUDE="-I /usr/local/include/ -I $PLZIPDIR -I$HOME/local/include/"
LIBS="-L /usr/local/lib/ -L /usr/lib/ -L$HOME/local/lib -llz -lpthread"
SRC="compress.cc dec_stream.cc dec_stdout.cc decompress.cc file_index.cc"
# echo "g++ -c -Wall -Werror -fPIC $INCLUDE $SRC $LIBS"
# g++ -c -Wall -Werror -fPIC $INCLUDE $SRC $LIBS

echo "g++ -shared -fPIC -o libplzip.so $SRC $INCLUDE $LIBS"
g++ -shared -fPIC -o libplzip.so $SRC $INCLUDE $LIBS

# make the lib available at runtime
export LD_LIBRARY_PATH=$PLZIPDIR:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=plzip/:$DYLD_LIBRARY_PATH

cd $ORIG
