#=============================-
# Makefile to compile a simple file format manager
# for converting a ConfigPack file to txt-format files, and vice versa.
#
# Author: Jaeuk Kim
# Date: September 2020
# Email: phy000.kim@gmail.com
#=============================

MACHINE = HOME
#MACHINE = DELLA

#==============================
# DELLA
#==============================
ifeq ($(MACHINE), DELLA)
#----- Directories ------------
IDIR = /usr/local/include/# include directory
ODIR = obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ./# directory for executables

#	Executable
TARGET = $(EXEdir)NR.out

#----- Flags and libraries -----
GSLLIB = -lgsl -lgslcblas
LIBS = -lm $(GSLLIB)
endif

#==============================
# MY LAPTOP
#==============================
ifeq ($(MACHINE), HOME)
$(info	$(MACHINE))
#----- Directories ------------
IDIR = -I/usr/include/ -I/home/jaeukk/lib/gsl/include# include directory
ODIR = ../obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ./# directory for executables

#	Executable
TARGET = $(EXEdir)1D-rsa.out#NR_xbootstrap.out# output file name

#----- Flags and libraries -----
GSLLIB = -L/home/jaeukk/lib/gsl/lib -lgsl -lgslcblas -lm
LIBS = -lm $(GSLLIB) 
endif

#==============================
#	Common flages
#==============================
CC = g++
DEBUG = -g
CFLAGS = -Wall -std=c++11 $(DEBUG) $(IDIR) -c -fopenmp -O3
LFLAGS = -Wall -std=c++11 $(DEBUG) $(LIBS) -fopenmp -O3

$(info	$(CFLAGS))
$(info	$(LFLAGS))

objs_core = GeometryVector.o PeriodicCellList.o RandomGenerator.o etc.o NumberVariance.o 
objs_netw = main.o NR_distribution.o
objs_misc = GenerateConfigs.o HardSphereFluids.o RandomSequentialAddition.o

DEPS = NR_distribution.h \
  ../GeometryVector.h ../PeriodicCellList.h ../RandomGenerator.h \
  ../etc.h ../NumberVariance.h \
  ../GenerateConfigs.h ../HardSphereFluids.h ../RandomSequentialAddition.h

OBJ_core = $(patsubst %.o, $(ODIR)%.o, $(objs_core))
OBJ_netw = $(patsubst %.o, $(ODIR)%.o, $(objs_netw))
OBJ_misc = $(patsubst %.o, $(ODIR)%.o, $(objs_misc))
OBJ = $(OBJ_core) $(OBJ_netw) $(OBJ_misc)

$(info	$(OBJ))

SRC_core = $(patsubst %.o, ../%.cpp, $(objs_core))
SRC_netw = $(patsubst %.o, %.cpp, $(objs_netw))
SRC_misc = $(patsubst %.o, ../%.cpp, $(objs_misc))

$(info	$(SRC_core))
$(info	$(SRC_netw))
$(info	$(SRC_misc))

#	echo 'core source files' $(SRC_core)
#	echo 'FFH source files' $(SRC_FFH)
#	$(SRC_core)
#	$(SRC_FFH)



all: $(OBJ) $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(PREPROC) -o $@ $(OBJ) $(LFLAGS)

# compilations...
$(OBJ_core): $(ODIR)%.o: $(addprefix ../, %.cpp) $(DEPS)
	echo $@ $<  
	$(CC) $(PREPROC) -o $@ $< $(CFLAGS)  

$(OBJ_misc): $(ODIR)%.o: $(addprefix ../, %.cpp) $(DEPS)
	echo $@ $<  
	$(CC) $(PREPROC) -o $@ $< $(CFLAGS)  

$(OBJ_netw): $(ODIR)%.o: $(notdir %.cpp) $(DEPS)
	echo $@ $<
	$(CC) $(PREPROC) -o $@ $< $(CFLAGS) 

#$(OBJ): $(ODIR)%.o: $(SRC) $(DEPS)
#	echo $@ $<
#	$(CC) -o $@ $< $(CFLAGS) 

clean:
	rm -f $(ODIR)*.o $(TARGET)

.PHONY: all
all:

