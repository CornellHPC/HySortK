NAME=test_batchsize

mkdir $NAME

for size in 10000 20000 40000 80000 160000 320000 640000 1280000 
do
    make clean
    make K=31 M=17 L=2 U=50 LOG=2 BATCH=${size} -j8
    mv ukmerc $NAME/ukmerc_${size}
done