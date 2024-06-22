all: raduls
	
OBJ_DIR = ../../obj

CC 	= g++
CFLAGS	= -IRaduls -Wall -fopenmp -O3 -mavx2 -fno-ipa-ra -fno-tree-vrp -fno-tree-pre -m64 -std=c++17

$(OBJ_DIR)/sorting_network.o : sorting_network.cpp 
	$(CC) -O1 -m64 -std=c++17 -mavx2 -pthread -c $< -o $@

raduls: $(OBJ_DIR)/sorting_network.o

clean:
	-rm $(OBJ_DIR)/sorting_network.o

