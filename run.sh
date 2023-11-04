GASNET_BACKTRACE=1 upcxx-run -shared-heap 16G -n 2 -N 1  ./ukmerc /pscratch/sd/y/yfli03/ELBA_dataset/acinetobacter_baumannii/reads.fa

upcxx-run -shared-heap 16G -n 16 -N 2  ./ukmerc /pscratch/sd/y/yfli03/ELBA_dataset/acinetobacter_baumannii/reads.fa