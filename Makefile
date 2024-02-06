K?=31
L?=15
U?=40
LOG?=2
D?=0
T?=4
BATCH?=250000
OVERLAP?=2
COMPILE_TIME_PARAMETERS=-DKMER_SIZE=$(K) -DLOWER_KMER_FREQ=$(L) -DUPPER_KMER_FREQ=$(U) -DLOG_LEVEL=$(LOG) -DDEBUG=$(D) -DTHREAD_PER_TASK=$(T) -DMAX_SEND_BATCH=$(BATCH) -DOVERLAP_LEVEL=$(OVERLAP)
OPT=

ifeq ($(D), 1)
OPT+=-g -O2 -fsanitize=address -fno-omit-frame-pointer
else ifeq ($(D), 2)		# Debugging with MAP
OPT+=-g1 -O3 -fno-inline -fno-optimize-sibling-calls
else
OPT+=-O3
endif

FLAGS=$(OPT) $(COMPILE_TIME_PARAMETERS) -DTHREADED -fopenmp -std=c++17 -I./include -I./src

COMPILER=CC

OBJECTS=obj/logger.o \
		obj/dnaseq.o \
		obj/dnabuffer.o \
		obj/fastaindex.o \
		obj/hashfuncs.o \
		obj/kmerops.o \
		obj/memcheck.o \


all: ukmerc

install: ukmerc
	cp ukmerc $(HOME)/bin

ukmerc: obj/main.o $(OBJECTS)
	$(COMPILER) $(FLAGS) -o $@ $^ -lz

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(COMPILER) $(FLAGS) -c -o $@ $<

obj/main.o: src/main.cpp include/logger.hpp include/timer.hpp include/dnaseq.hpp include/dnabuffer.hpp include/fastaindex.hpp include/kmerops.hpp include/memcheck.hpp include/compiletime.h
obj/logger.o: src/logger.cpp include/logger.hpp
obj/dnaseq.o: src/dnaseq.cpp include/dnaseq.hpp
obj/dnabuffer.o: src/dnabuffer.cpp include/dnabuffer.hpp include/dnaseq.hpp
obj/fastaindex.o: src/fastaindex.cpp include/fastaindex.hpp include/dnaseq.hpp include/dnabuffer.hpp
obj/hashfuncs.o: src/hashfuncs.cpp include/hashfuncs.hpp
obj/kmerops.o: src/kmerops.cpp include/kmerops.hpp include/kmer.hpp include/dnaseq.hpp include/logger.hpp include/timer.hpp include/dnabuffer.hpp include/paradissort.hpp include/memcheck.hpp
obj/memcheck.o: src/memcheck.cpp include/memcheck.hpp

clean:
	rm -rf *.o obj/* ukmerc $(HOME)/bin/ukmerc