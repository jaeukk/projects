#=============================-
# Makefile to compile example codes
# for computing structure factors and 
# pair correlation functions
#=============================

#==============================
# directories
#==============================
IDIR = -I/usr/include/ -I/home/jaeukk/lib/gsl/include # include directory
ODIR = ./obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ./# directory for executables

TARGET = $(EXEdir)example.out# output file name

#==============================
#	Flags and libraries
#==============================
CC = g++
DEBUG = -g
LIBS = -lm -lgsl -lgslcblas

CFLAGS = -Wall -std=c++0x $(DEBUG) -I$(IDIR) -c #-O3 -fopenmp
LFLAGS = -Wall -std=c++0x $(DEBUG) $(LIBS) #-O3 -fopenmp

#
CORE_objs = GeometryVector.o PeriodicCellList.o RandomGenerator.o \
  etc.o #PairCorrelation.o StructureFactor.o
CORE_OBJS = $(patsubst %.o, $(ODIR)%.o, $(CORE_objs))
CORE_SRCS = $(patsubst %.o, ../%.cpp, $(CORE_OBJS))
CORE_HDRS = $(CORE_SRCS:.cpp=.h)

MAIN = main.cpp
MOBJ = $(patsubst %.cpp, $(ODIR)%.o, $(MAIN))



objs =main.o GeometryVector.o PeriodicCellList.o RandomGenerator.o \
  etc.o #PairCorrelation.o StructureFactor.o 
  
OBJ = $(patsubst %.o, $(ODIR)%.o, $(objs))
	 
DEPS = GeometryVector.h PeriodicCellList.h RandomGenerator.h \
  etc.h #PairCorrelation.h StructureFactor.h


all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LFLAGS)

$(OBJ): $(ODIR)%.o: $(notdir %.cpp) $(DEPS)
	echo $@ $<
	$(CC) -o $@ $< $(CFLAGS) 

clean:
	rm -f $(ODIR)*.o $(TARGET)
