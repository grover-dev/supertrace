
SRC = ./src
INC = ./inc
OBJ = ./obj
CC = nvcc

# card code reference: https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/
CFLAGS = -O3 -std=c++11 -lcuda --gpu-architecture=compute_75 --gpu-code=sm_75

SOURCES := $(wildcard $(SRC)/*.cu) 
INCLUDES := $(wildcard $(SRC)/*.cuh)
OBJECTS := $(patsubst $(SRC)/%.cu, $(OBJ)/%.o, $(SOURCES))

$(OBJ)/%.o: $(SRC)/%.cu
	$(CC)  -c $< -o $@ $(CFLAGS) -dc

raytracer: $(OBJECTS) 
	$(CC) -o $@ $^  $(CFLAGS) -I$(INC)


.PHONY: clean
# *~ core $(INCDIR)/*~
clean:
	rm -f $(OBJ)/*.o
	rm raytracer
	
