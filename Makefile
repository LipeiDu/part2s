#compiler
#CC = g++ #
#CC = pgc++ #for OpenACC GPU acceleration
CC = icpc #intel compiler tends to be faster
#CC = nvcc # to compile cuda 

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
#CFLAGS = -acc #GPU acceleration with pgc++ 
#CFLAGS = -acc -ta=tesla:managed #also uses unified memory
CFLAGS = -O3  
#LIBS= 

#build target
SRC = src
TARGET = part2s

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) -o $(TARGET) $(SRC)/$(TARGET).cpp
clean:
