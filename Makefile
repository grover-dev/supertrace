
SRC = ./src
INC = ./inc
OBJ = ./obj
CC = nvcc

# card code reference: https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/
CFLAGS = -g -std=c++11 -lcuda --gpu-architecture=compute_60 --gpu-code=sm_60

SOURCES := $(wildcard $(SRC)/*.cu) 
INCLUDES := $(wildcard $(SRC)/*.cuh)
OBJECTS := $(patsubst $(SRC)/%.cu, $(OBJ)/%.o, $(SOURCES))

MPIFLAGS = -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -lpthread -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi


# $(OBJ)/%.o: $(SRC)/%.cu
# 	$(CC)  -c $< -o $@ $(CFLAGS) -dc

# raytracer: $(OBJECTS) 
# 	$(CC) -o $@ $^  $(CFLAGS) -I$(INC)

all:
	nvcc $(CFLAGS) -c src/shapes.cu -dc -lcuda -o obj/shapes.o
	nvcc $(CFLAGS) $(MPIFLAGS) -c src/stl_raytracer.cu -dc -lcuda -o obj/stl_raytracer.o
	nvcc $(CFLAGS) obj/stl_raytracer.o obj/shapes.o -lm -lcudart  $(MPIFLAGS) -o raytracer -I ./src/






.PHONY: clean
# *~ core $(INCDIR)/*~
clean:
	rm -f $(OBJ)/*.o
	rm raytracer
	
