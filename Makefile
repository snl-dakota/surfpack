## ----------------------------------------------------------
## Project: SURFPACK++
##
## File:        Makefile
## Author:      Mark Richards 
## Created:     February 2002
##
## Description:
## ++ A makefile for the SURFPACK++ library
## ----------------------------------------------------------

.SUFFIXES: .cpp .f .o

# CC=CC
# LINK=CC
 
# CPPFLAGS=-g -D__TESTING_MODE__
# LINKFLAGS = -G 

FC=g77
CC=g++
LINK=ar

#CPPFLAGS=-Wall -g -fPIC
#CPPFLAGS=-Wall -g -D__TESTING_MODE__ -pg -O3

CPPFLAGS=-Wall -g -D__TESTING_MODE__
LINKFLAGS = rcs 

DIR=$(SURFPACK)
INCLUDE=  -I$(DIR)/src -I$(DIR)/src/surfaces -I$(DIR)/src/surfaces/surfaceComponents -I$(DIR)/src/iterators -I$(DIR)/src/kriging

LDFLAGS= -L$(SURFPACK)/lib
LDLIBS= -lsurfpack -lkriging -lconmin -llapack -lblas -lm -lg2c 

CONMIN_SRCS = $(DIR)/src/conmin/conmin.f

CONMIN_OBJS = $(CONMIN_SRCS:.f=.o)

CONMIN = $(DIR)/lib/libconmin.a

KRIG_SRCS = $(DIR)/src/kriging/krig_model.f \
       $(DIR)/src/kriging/cpp_f77_conmin_cntl.f \
       $(DIR)/src/kriging/cpp_f77_dotcntl.f 
KRIG_OBJS = $(KRIG_SRCS:.f=.o)

.f.o:
	$(FC) -c $< -o $@

KRIG = $(DIR)/lib/libkriging.a 

#the last few files may be in a separate ParticleSwarm library
SURF_SRCS = $(DIR)/src/SurfPoint.cpp \
       $(DIR)/src/SurfData.cpp \
       $(DIR)/src/Surface.cpp \
       $(DIR)/src/iterators/AbstractSurfDataIterator.cpp \
       $(DIR)/src/iterators/SurfDataIterator.cpp \
       $(DIR)/src/iterators/SkipSurfDataIterator.cpp \
       $(DIR)/src/surfaces/KrigingSurface.cpp \
       $(DIR)/src/surfaces/PolynomialSurface.cpp 

SURF_OBJS = $(SURF_SRCS:.cpp=.o)

OTHER_SRCS = $(DIR)/interface/main.cpp \
	     $(DIR)/interface/main.cpp

OTHER_OBJS = $(OTHER_SRCS:.cpp=.o)

SURF = $(DIR)/lib/libsurfpack.a
NEWBIN = $(DIR)/bin/surfpack

.cpp.o:
	$(CC) -c $(CPPFLAGS) $(INCLUDE) $< -o $@


all: $(CONMIN) $(KRIG) $(SURF) $(NEWBIN)

$(CONMIN) : $(CONMIN_OBJS)
	$(LINK) $(LINKFLAGS) $@ $(CONMIN_OBJS)
	
$(KRIG): $(KRIG_OBJS)
	$(LINK) $(LINKFLAGS) $@ $(KRIG_OBJS) 

$(SURF): $(SURF_OBJS)
	$(LINK) $(LINKFLAGS)  $@ $(SURF_OBJS) 

$(NEWBIN): $(KRIG) $(SURF) $(OTHER_OBJS)
	$(CC) $(LDFLAGS) -o $@ $(DIR)/interface/main.o $(LDLIBS)

clean:
	rm $(CONMIN) $(KRIG) $(SURF) $(NEWBIN) $(SURF_OBJS) $(KRIG_OBJS) $(CONMIN_OBJS) $(OTHER_OBJS)
	   











#SRCS = $(DIR)/src/SurfPoint.cpp \
#       $(DIR)/src/SurfData.cpp \
#       $(DIR)/src/iterators/DataIterator.cpp \
#       $(DIR)/src/iterators/SkipDataIterator.cpp \
#       $(DIR)/src/iterators/SurfDataIterator.cpp \
#       $(DIR)/src/iterators/SkipSurfDataIterator.cpp \
#       $(DIR)/src/Surface.cpp \
#       $(DIR)/src/surfaces/KrigingSurface.cpp \
#       $(DIR)/src/surfaces/CubicSurface.cpp \
#       $(DIR)/src/surfaces/ANNSurface.cpp \
#       $(DIR)/src/surfaces/QuadSurface.cpp \
#       $(DIR)/src/surfaces/LinearSurfaceFactory.cpp \
#       $(DIR)/src/surfaces/LinearSurface.cpp \
#       $(DIR)/src/surfaces/surfaceComponents/MultiLayerPerceptron.cpp \
#       $(DIR)/src/swarm/OptimizationProblem.cpp \
#       $(DIR)/src/swarm/Particle.cpp \
#       $(DIR)/src/swarm/Swarm.cpp \
#       $(DIR)/src/swarm/NeuralNetwork.cpp \
#       $(DIR)/src/swarm/ANNWeights.cpp \
#       $(DIR)/src/swarm/Rastrigin.cpp 
