#compiler
CC = g++ #
#CC = icpc #intel compiler tends to be faster
#CC = nvcc # to compile cuda 

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
CFLAGS = -fopenmp -lhdf5 -lhdf5_cpp -lgsl -lgslcblas 
LIBS= -lgomp -lhdf5 -lhdf5_cpp -lgsl -lgslcblas

#source code directory
SRC = cpp

#build target
TARGET = part2s

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) -o $(TARGET) $(SRC)/$(TARGET).cpp
clean:
