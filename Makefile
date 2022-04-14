
SRC = ./src
INC = ./inc
OBJ = ./obj
CC = g++

CFLAGS = -g

SOURCES := $(wildcard $(SRC)/*.cpp) 
INCLUDES := $(wildcard $(INC)/*.hpp)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))

$(OBJ)/%.o: $(SRC)/%.c
	$(CC)  -c $< -o $@ $(CFLAGS)


raytracer: $(OBJECTS) 
	$(CC) -o $@ $^ -I$(INC)


.PHONY: clean
# *~ core $(INCDIR)/*~
clean:
	rm -f $(OBJ)/*.o
	rm raytracer
	
