
SRC = ./src
INC = ./inc
OBJ = ./obj
CC = nvcc

CFLAGS = -O3 -std=c++11 --gpu-architecture=compute_37 --gpu-code=sm_37

SOURCES := $(wildcard $(SRC)/*.cu) 
INCLUDES := $(wildcard $(INC)/*.hpp)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))

$(OBJ)/%.o: $(SRC)/%.c
	$(CC)  -c $< -o $@ $(CFLAGS)


raytracer: $(OBJECTS) 
	$(CC) -o $@ $^ -I$(INC) $(CFLAGS)


.PHONY: clean
# *~ core $(INCDIR)/*~
clean:
	rm -f $(OBJ)/*.o
	rm raytracer
	
