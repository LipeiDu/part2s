#compiler
#CC = g++ #
CC = icpc #intel compiler tends to be faster
#CC = nvcc # to compile cuda 

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
CFLAGS = -lhdf5 -lhdf5_cpp  
#LIBS= 

#source code directory
SRC = src

#build target
TARGET = getEnergyUrQMD

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) -o $(TARGET) $(SRC)/$(TARGET).cpp
clean:
