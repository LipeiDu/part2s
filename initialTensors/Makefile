#compiler
#CC = g++ #
#CC = icpc #intel compiler tends to be faster
CC = nvcc # to compile cuda 

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
CFLAGS = -lhdf5 -lhdf5_cpp  
#LIBS= 

#source code directory
SRC = src

#build target
TARGET = ICpart2s

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cu
	$(CC) $(CFLAGS) $(LIBS) -o $(TARGET) $(SRC)/$(TARGET).cu
clean:
