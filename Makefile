K?=31
L?=15
U?=35
LOG?=2
D?=0
COMPILE_TIME_PARAMETERS=-DKMER_SIZE=$(K) -DLOWER_KMER_FREQ=$(L) -DUPPER_KMER_FREQ=$(U) -DLOG_LEVEL=$(LOG) -DDEBUG=$(D)
OPT=

ifeq ($(D), 1)
OPT+=-g -O2 -fsanitize=address -fno-omit-frame-pointer
else
OPT+=-O3
endif

FLAGS=$(OPT) $(COMPILE_TIME_PARAMETERS) -DTHREADED -Wall -std=c++17 -I./include -I./src

COMPILER=upcxx

OBJECTS=obj/logger.o \
		obj/timer.o \
		obj/dnaseq.o \
		obj/dnabuffer.o \
		obj/fastaindex.o \
		obj/hashfuncs.o \
		obj/kmer.o


all: ukmerc

install: ukmerc
	cp ukmerc $(HOME)/bin

ukmerc: obj/main.o $(OBJECTS)
	@echo upcxx -c -o $@ $^
	$(COMPILER) $(FLAGS) -o $@ $^

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	@echo upcxx $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	$(COMPILER) $(FLAGS) -c -o $@ $<

obj/main.o: src/main.cpp include/logger.hpp include/timer.hpp include/dnaseq.hpp include/dnabuffer.hpp include/fastaindex.hpp include/kmer.hpp include/hashfuncs.hpp include/compiletime.h
obj/logger.o: src/logger.cpp include/logger.hpp
obj/timer.o: src/timer.cpp include/timer.hpp
obj/dnaseq.o: src/dnaseq.cpp include/dnaseq.hpp
obj/dnabuffer.o: src/dnabuffer.cpp include/dnabuffer.hpp include/dnaseq.hpp
obj/fastaindex.o: src/fastaindex.cpp include/fastaindex.hpp include/dnaseq.hpp include/dnabuffer.hpp
obj/hashfuncs.o: src/hashfuncs.cpp include/hashfuncs.hpp
obj/kmer.o: src/kmer.cpp include/kmer.hpp include/dnaseq.hpp include/hashfuncs.hpp include/compiletime.h

clean:
	rm -rf *.o obj/* ukmerc $(HOME)/bin/ukmerc