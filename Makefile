# Kmer Size
K?=31
# Minimizer Size
M?=17
# Lower Kmer Frequency
L?=15
# Upper Kmer Frequency
U?=40
# Extension Information
EXT?=0

# Log Level (1, 2, 3, 4)
LOG?=2
# Debugging
D?=0

# Thread per Worker
T?=4
# Max Thread Memory Bounded
T2?=16
# Average Task per Worker
TPW?=3
# SORTING OPTION (0: runtime decision, 1: paradis, 2: raduls)
SORT?=0
# Max Send Batch
BATCH?=80000

# Dispatch Upper Coefficient
DISPATCH_UPPER = 1.5
# Dispatch Step
DISPATCH_STEP = 0.05
# Unbalanced Threshold
UNBALANCED_THRESHOLD = 2.3

PLAIN_CLASSIFIER=0

PLAIN_DISPATCHER=0

COMPILE_TIME_PARAMETERS=-DKMER_SIZE=$(K) -DMINIMIZER_SIZE=$(M) \
	-DLOWER_KMER_FREQ=$(L) -DUPPER_KMER_FREQ=$(U) -DLOG_LEVEL=$(LOG) -DDEBUG=$(D) \
	-DTHREAD_PER_WORKER=$(T) -DMAX_SEND_BATCH=$(BATCH) -DMAX_THREAD_MEMORY_BOUNDED=$(T2) \
	-DSORT=$(SORT) -DAVG_TASK_PER_WORKER=$(TPW) \
	-DDISPATCH_UPPER_COE=$(DISPATCH_UPPER) -DDISPATCH_STEP=$(DISPATCH_STEP) \
	-DUNBALANCED_RATIO=$(UNBALANCED_THRESHOLD) \
	-DPLAIN_CLASSIFIER=$(PLAIN_CLASSIFIER) -DPLAIN_DISPATCHER=$(PLAIN_DISPATCHER) \
	-DEXTENSION=$(EXT)
OPT=

# Check if M is less than K
ifneq ($(shell test $(M) -lt $(K) && echo 0 || echo 1), 0)
$(error ERROR: MINIMIZER_SIZE (M) must be less than KMER_SIZE (K))
endif

ifeq ($(D), 1)
OPT+=-g -O2 -fsanitize=address -fno-omit-frame-pointer
else ifeq ($(D), 2)
OPT+=-g1 -O3 -fno-inline -fno-optimize-sibling-calls
else
OPT+=-O3
endif

FLAGS=-pthread -m64 -mavx2 -DTHREADED -fopenmp -std=c++17 -I./include -I./src -I./dependency/Raduls -I./dependency/Paradis
LINK=-lm -fopenmp -O3 -mavx2 -fno-ipa-ra -fno-tree-vrp -fno-tree-pre  -std=c++17 -lpthread -DTHREADED

ifeq ($(NERSC_HOST), perlmutter)
COMPILER=CC
else
COMPILER=CXX
endif

LINKER=ld

OBJECTS=obj/logger.o \
		obj/dnaseq.o \
		obj/dnabuffer.o \
		obj/fastaindex.o \
		obj/hashfuncs.o \
		obj/kmerops.o \
		obj/memcheck.o \
		obj/hysortk.o



all: print lib

print:
	$(info ------ HySortK Compiletime Parameters ------ )
	$(info LOG_LEVEL: $(LOG), DEBUG: $(D))
	$(info KMER_SIZE: $(K), MINIMIZER_SIZE: $(M), EXTENSION: $(EXT))
	$(info LOWER_KMER_FREQ: $(L), UPPER_KMER_FREQ: $(U))
	$(info THREAD_PER_WORKER: $(T), AVG_TASK_PER_WORKER: $(TPW))
	$(info MAX_THREAD_MEMORY_BOUNDED: $(T2), MAX_SEND_BATCH: $(BATCH))
	$(info DISPATCH_UPPER_COE: $(DISPATCH_UPPER), DISPATCH_STEP: $(DISPATCH_STEP))
	$(info SORT: $(SORT) (0: Runtime decision, 1: PARADIS, 2: RADULS))
	$(info UNBALANCED_RATIO: $(UNBALANCED_THRESHOLD))
	$(info PLAIN_CLASSIFIER: $(PLAIN_CLASSIFIER), PLAIN_DISPATCHER: $(PLAIN_DISPATCHER))
	$(info ------------------------------------------- )

lib: $(OBJECTS)
	$(MAKE) -C dependency/Raduls
	$(LINKER) -r -o libhysortk.o $(OBJECTS) obj/sorting_network.o

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(COMPILER) $(OPT) $(COMPILE_TIME_PARAMETERS) $(FLAGS) -c -o $@ $<

standalone: all
	$(COMPILER) $(OPT) $(COMPILE_TIME_PARAMETERS) $(FLAGS) -c -o obj/standalone.o standalone/main.cpp
	$(COMPILER) $(OPT) $(COMPILE_TIME_PARAMETERS) $(FLAGS) -o hysortk obj/standalone.o libhysortk.o


obj/hysortk.o: src/hysortk.cpp include/logger.hpp include/timer.hpp include/dnaseq.hpp include/dnabuffer.hpp include/fastaindex.hpp include/kmerops.hpp include/memcheck.hpp include/compiletime.h 
obj/logger.o: src/logger.cpp include/logger.hpp
obj/dnaseq.o: src/dnaseq.cpp include/dnaseq.hpp
obj/dnabuffer.o: src/dnabuffer.cpp include/dnabuffer.hpp include/dnaseq.hpp
obj/fastaindex.o: src/fastaindex.cpp include/fastaindex.hpp include/dnaseq.hpp include/dnabuffer.hpp
obj/hashfuncs.o: src/hashfuncs.cpp include/hashfuncs.hpp
obj/kmerops.o: src/kmerops.cpp include/kmerops.hpp include/kmer.hpp include/dnaseq.hpp include/logger.hpp include/timer.hpp include/dnabuffer.hpp include/memcheck.hpp dependency/Paradis/paradissort.hpp
obj/memcheck.o: src/memcheck.cpp include/memcheck.hpp

clean:
	rm -rf *.o obj/* hysortk $(HOME)/bin/hysortk