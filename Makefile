K?=31
M?=13
L?=15
U?=40
LOG?=2
D?=0
T?=4
T2?=16
SORT?=0
BATCH?=250000
COMPILE_TIME_PARAMETERS=-DKMER_SIZE=$(K) -DMINIMIZER_SIZE=$(M) -DLOWER_KMER_FREQ=$(L) -DUPPER_KMER_FREQ=$(U) -DLOG_LEVEL=$(LOG) -DDEBUG=$(D) -DTHREAD_PER_TASK=$(T) -DMAX_SEND_BATCH=$(BATCH) -DMAX_THREAD_MEMORY_BOUNDED=$(T2) -DSORT=$(SORT)
OPT=

# TODO: check if M is less than K

ifeq ($(D), 1)
OPT+=-g -O2 -fsanitize=address -fno-omit-frame-pointer
else ifeq ($(D), 2)		# Debugging with MAP
OPT+=-g1 -O3 -fno-inline -fno-optimize-sibling-calls
else
OPT+=-O3
endif

FLAGS=-pthread -m64 -mavx2 -DTHREADED -fopenmp -std=c++17 -I./include -I./src

COMPILER=CC

OBJECTS=obj/logger.o \
		obj/dnaseq.o \
		obj/dnabuffer.o \
		obj/fastaindex.o \
		obj/hashfuncs.o \
		obj/kmerops.o \
		obj/memcheck.o \
		raduls/sorting_network.o


all: ukmerc

install: ukmerc
	cp ukmerc $(HOME)/bin

ukmerc: obj/main.o $(OBJECTS)
	$(COMPILER) $(OPT) $(COMPILE_TIME_PARAMETERS) $(FLAGS) -o $@ $^ -lz

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(COMPILER) $(OPT) $(COMPILE_TIME_PARAMETERS) $(FLAGS) -c -o $@ $<

raduls/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(COMPILER) $(OPT) $(COMPILE_TIME_PARAMETERS) $(FLAGS) -c -o $@ $<

obj/main.o: src/main.cpp include/logger.hpp include/timer.hpp include/dnaseq.hpp include/dnabuffer.hpp include/fastaindex.hpp include/kmerops.hpp include/memcheck.hpp include/compiletime.h src/sorting_network.cpp include/raduls.h include/record.h include/small_sort.h include/sorting_network.h include/exceptions.h include/defs.h include/comp_and_swap.h
obj/logger.o: src/logger.cpp include/logger.hpp
obj/dnaseq.o: src/dnaseq.cpp include/dnaseq.hpp
obj/dnabuffer.o: src/dnabuffer.cpp include/dnabuffer.hpp include/dnaseq.hpp
obj/fastaindex.o: src/fastaindex.cpp include/fastaindex.hpp include/dnaseq.hpp include/dnabuffer.hpp
obj/hashfuncs.o: src/hashfuncs.cpp include/hashfuncs.hpp
obj/kmerops.o: src/kmerops.cpp include/kmerops.hpp include/kmer.hpp include/dnaseq.hpp include/logger.hpp include/timer.hpp include/dnabuffer.hpp include/paradissort.hpp include/memcheck.hpp src/sorting_network.cpp include/raduls.h include/record.h include/small_sort.h include/sorting_network.h include/exceptions.h include/defs.h include/comp_and_swap.h
obj/memcheck.o: src/memcheck.cpp include/memcheck.hpp
raduls/sorting_network.o: src/sorting_network.cpp include/raduls.h include/record.h include/small_sort.h include/sorting_network.h include/exceptions.h include/defs.h include/comp_and_swap.h

clean:
	rm -rf *.o obj/* raduls/* ukmerc $(HOME)/bin/ukmerc