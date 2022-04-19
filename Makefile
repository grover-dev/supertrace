
SRC = ./src
INC = ./inc
OBJ = ./obj
CC = nvcc

CFLAGS = -O3 -std=c++11  --gpu-architecture=compute_37 --gpu-code=sm_37

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
	
