#=============================-
# Makefile to compile a simple program to manipulate QC configurations
#
# Author: Jaeuk Kim
# Date: October 2024
# Email: phy000.kim@gmail.com
#=============================

core = /home/jaeukk/new_codes/c++/cores/

MACHINE = HOME
#MACHINE = DELLA

#==============================
# DELLA
#==============================
ifeq ($(MACHINE), DELLA)
#echo $(MACHINE)
$(info	$(MACHINE))
#----- Directories ------------
IDIR = -I/usr/local/include/ -I/usr/include/ -I$(core) #-I/home/jaeukk/install/include # include directory
ODIR = obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ~/EXC/# directory for executables

#	Executable
TARGET = $(EXEdir)QC.out

#----- Flags and libraries -----
GSLLIB = -L/usr/lib64 -lgsl -lgslcblas
NLOPTLIB = -L/usr/lib64 -lnlopt #-L/home/jaeukk/install/lib64 -lnlopt
LIBS = -lm $(GSLLIB) $(NLOPTLIB)
PREPRO = -DDELLA
endif

#==============================
# MY LAPTOP
#==============================
ifeq ($(MACHINE), HOME)
#echo $(MACHINE)
$(info	$(MACHINE))
#----- Directories ------------
IDIR = -I/usr/include/ -I/home/jaeukk/lib/gsl/include -I/usr/local/include/ -I$(core) # include directory
ODIR = obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ./# directory for executables

#	Executable
TARGET = $(EXEdir)QC.out

#----- Flags and libraries -----
GSLLIB = -L/home/jaeukk/lib/gsl/lib -lgsl -lgslcblas
NLOPTLIB = -L/usr/local/lib -lnlopt
LIBS = -lm $(GSLLIB) $(NLOPTLIB)
PREPRO = 
endif

#==============================
#	Common flages
#==============================
CC = g++
DEBUG = -g
CFLAGS = -Wall -std=c++11 $(DEBUG) $(IDIR) -c -fopenmp #-O3
LFLAGS = -Wall -std=c++11 $(DEBUG) $(LIBS) -fopenmp #-O3

$(info	$(CFLAGS))
$(info	$(LFLAGS))

objs_core = GeometryVector.o RandomGenerator.o etc.o PeriodicCellList.o 
objs_main = main.o #ConfigDiff.o

DEPS = $(core)GeometryVector.h $(core)RandomGenerator.h $(core)etc.h $(core)PeriodicCellList.h ConfigDiff.h Perturbed.h

OBJ_core = $(patsubst %.o, $(ODIR)%.o, $(objs_core))
OBJ_main = $(patsubst %.o, $(ODIR)%.o, $(objs_main))
OBJ = $(OBJ_core)  $(OBJ_main) 

$(info	$(OBJ))

SRC_core = $(patsubst %.o, $(core)%.cpp, $(objs_core))
SRC_main = $(patsubst %.o, %.cpp, $(objs_main))

$(info	$(SRC_core))
$(info	$(SRC_main))
#$(info	$(SRC_pot))



all: $(OBJ) $(TARGET)

# compilations...
$(OBJ_core): $(ODIR)%.o: $(addprefix $(core), %.cpp) $(DEPS)
	echo $@ $<  
	$(CC) $(PREPRO) -o $@ $< $(CFLAGS)  

$(OBJ_main): $(ODIR)%.o: $(notdir %.cpp) $(DEPS)
	echo $@ $<
	$(CC) $(PREPRO) -o $@ $< $(CFLAGS) 

	echo "Start linking"

$(TARGET): $(OBJ)
	$(CC) $(PREPRO) -o $@ $(OBJ) $(LFLAGS)


# $(OBJ_pot): $(ODIR)%.o: $(addprefix ../, %.cpp) $(DEPS)
# 	echo $@ $<  
# 	$(CC) $(PREPRO) -o $@ $< $(CFLAGS)  


#$(OBJ): $(ODIR)%.o: $(SRC) $(DEPS)
#	echo $@ $<
#	$(CC) -o $@ $< $(CFLAGS) 

clean:
	rm -f $(ODIR)*.o $(TARGET)

.PHONY: all
all:

