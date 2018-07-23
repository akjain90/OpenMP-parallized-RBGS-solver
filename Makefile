CXX = g++
CXXFLAGS = -std=c++0x -Wall -Wextra -ffast-math -march=native -Werror -Wshadow -O3 -fopenmp -DNDEBUG

INCLUDES =
LDFLAGS =
LIBS =

TARGET = rbgs
OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

$(TARGET).o: $(TARGET).cpp matrix.h matrix.cpp Timer.h Makefile 
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(TARGET).cpp

clean:
	@$(RM) -rf *.o $(TARGET)
