IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii
LIBS += -ltbb

CXX = g++
CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
#CXXFLAGS += -O3
CXXFLAGS += -std=c++11


PROGRAM = testexe

OBJS = main.o
OBJS += TwoInclusionLift_cell_mesh.o

DEPS = TwoInclusionLift.h

#link the object files into an executable
$(PROGRAM) : $(OBJS) #list all object files
	$(CXX) -o $@ $(CXXFLAGS) $(LFLAGS) $^ $(LIBS)

#make the object files '.o's
%.o : %.cpp $(DEPS)
	$(CXX) -o $@ $(CXXFLAGS) $(IFLAGS) -c $<

