CC=g++
CFLAGS=-std=c++11 -O3 -w
#CFLAGSLIB=-shared -fPIC
#PLIBZ=
BIN=bin
LDFLAGS=-L /usr/local/lib/ -L plzip/ -L ~/lib/
SRC=src/Referee.cpp
INCLUDE=-I include/ -I /usr/local/include/ -I plzip/ -I ~/include/
#LIBS=-lstaden-read -lplzip # -llz
LIBS=-lstaden-read -lpthread -lplzip
EXE=referee

all: $(EXE)

#plibz: 
#	$(CC) $(CFLAGSLIB) -o lib$@.so compress.cc dec_stream.cc dec_stdout.cc decompress.cc file_index.cc -I /usr/local/include/ -I plzip -L /usr/local/lib/ -L /usr/lib/ -llz -lpthread

$(EXE):
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BIN)/$@ $(SRC) $(INCLUDE) $(LIBS)
clean:
	rm -f $(BIN)/$(EXE)

