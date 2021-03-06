CC=g++
CFLAGS=-std=c++11 -std=gnu++11 -O3 -w
# CCPARALL=-fopenmp -D_GLIBCXX_PARALLEL
#CFLAGSLIB=-shared -fPIC
#PLIBZ=
BIN=bin
LDFLAGS=-L /usr/local/lib/ -L plzip/
SRC=src/Referee.cpp
SRCSUPP=src/RefereeSupportTools.cpp
INCLUDE=-I include/ -I /usr/local/include/ -I plzip/ -I ~/include/ -I $(HOME)/lzlib/lzlib-1.6/
#LIBS=-lstaden-read -lplzip # -llz
LIBS=-lstaden-read -lpthread -lplzip
EXE=referee

# plzip so functionality
PLZIPDIR=plzip
PLZINCLUDE="-I /usr/local/include/ -I $PLZIPDIR -I$HOME/local/include/"
PLZLIBS="-L /usr/local/lib/ -L /usr/lib/ -L$HOME/local/lib -llz -lpthread"
PLZSRC="compress.cc dec_stream.cc dec_stdout.cc decompress.cc file_index.cc"

#SYSTEM=macos
SYSTEM=linux

# for mac: need to update DYLD_FALLBACK_LIBRARY_PATH to point to the directory containing TBB's *.dylib
# MAC
ifeq ($(SYSTEM),macos)
		# why do we need this include?
        # TBBINCL=-I $(HOME)/tbb/tbb43_20140724oss/include/ -I $(HOME)/tbb/tbb43_20140724oss/examples/common/utility/
        # TBBLIBS=-L $(HOME)/tbb/tbb43_20140724oss/build/macos_intel64_gcc_cc4.9.0_os10.9.5_release/ -ltbb
        # TBBINCL=-I $(HOME)/tbb/tbb43_20150424oss/include/
        # TBBLIBS=-L $(HOME)/tbb/tbb43_20150424oss/lib/ -ltbb
        # export DYLD_FALLBACK_LIBRARY_PATH=/Users/giantlynx/tbb/tbb43_20140724oss/build/macos_intel64_gcc_cc4.9.0_os10.9.5_release/
else
        # Linux
	# TBBLIBS=-L /opt/local/stow/intel_tbb_42_20130725oss/lib/ -ltbb
	# TBBINCL=-I /opt/local/stow/intel_tbb_42_20130725oss/include/
        # TBBLIBS=-L /opt/local/lib/ -ltbb
endif


all: $(EXE)

$(EXE):
	$(CC) $(CFLAGS) $(CCPARALL) $(LDFLAGS) -o $(BIN)/$@ $(SRC) $(INCLUDE) $(TBBINCL) $(LIBS) $(TBBLIBS)
clean:
	rm -f $(BIN)/$(EXE)
rsupport:
	$(CC) $(CFLAGS) $(CCPARALL) $(LDFLAGS) -o $(BIN)/$@ $(SRCSUPP) $(INCLUDE) $(LIBS) $(TBBLIBS)
plzipso:
	$(CC) -shared -fPIC -o libplzip.so $SRC $PLZINCLUDE $PLZLIBS
	export DYLD_LIBRARY_PATH=plzip/:$DYLD_LIBRARY_PATH