#########################
## Library Directories ##
########################


####################
## Makefile Setup ##
####################

# Get the host-name if empty
ifeq ($(host-name),)
	host-name := $(shell hostname)
endif
# Get the kernel-name if empty
ifeq ($(kernel-name),)
	kernel-name := $(shell uname -s)
endif
# Get the arch-name if empty
ifeq ($(arch-name),)
	arch-name := $(shell uname -p)
endif

# Define the C++ compiler to use
CXX := /usr/local/opt/llvm/bin/clang++
CXX += -Xclang -fopenmp

# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPSFILE = $(DEPSDIR)/$(notdir $*.d)
DEPSFLAGS = -MD -MF $(DEPSFILE) #-MP

# Define any directories containing header files
#   To include directories use -Ipath/to/files
CPPFLAGS += -I/usr/local/include/
CPPFLAGS += -I. -I/usr/local/Cellar/sfml/2.5.1_1/include/
CPPFLAGS += -I/Users/sashayousefi/Desktop/cme212-sashayousefi/mtl4
CPPFLAGS += -I/usr/local/Cellar/boost/1.78.0/include/
CPPFLAGS += -I/usr/local/Cellar/libomp/13.0.1/include/
CPPFLAGS += -I/Users/sashayousefi/Desktop/cme212-sashayousefi/thrust

# Define cxx compile flags
CXXFLAGS := -std=c++14 -O3 -W -Wall -Wextra
CXXFLAGS += -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP

# Define any directories containing libraries
#   To include directories use -Lpath/to/files
LDFLAGS += -L/usr/local/lib/ 
LDFLAGS += -L/usr/local/Cellar/sfml/2.5.1_1/lib/
LDFLAGS += -L/usr/local/Cellar/boost/1.78.0/lib/
LDFLAGS += -L/usr/local/Cellar/libomp/13.0.1/lib -lomp

# Define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX
ifeq ($(kernel-name), Linux)
  LDLIBS += -lsfml-graphics -lsfml-window -lsfml-system -lX11 -lGL -lpthread
endif
ifeq ($(kernel-name), Darwin)
  LDLIBS += -L/usr/local/opt/sfml/lib
  LDLIBS += -L/usr/local/opt/libomp/lib
  LDLIBS += -lsfml-graphics -lsfml-window -lsfml-system -framework OpenGL
endif

####################
## DEFAULT Rules ##
####################

# Suffix replacement rules
#   $^: the name of the prereqs of the rule
#   $<: the name of the first prereq of the rule
#   $@: the name of the target of the rule

# 'make' - default rule
all: $(EXEC)

# Default rule for creating an exec of $(EXEC) from a .o file
#$(EXEC):
% : %.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Default rule for creating a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -o $@ -c $<

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	-$(RM) *.o *.dSYM  $(EXEC) gtest.a gtest_main.a
	$(RM) -r $(DEPSDIR)

# Define rules that do not actually generate the corresponding file
.PHONY: clean all

# Include the dependency files
-include $(wildcard $(DEPSDIR)/*.d)

########################
## googletest Options ##
########################

# Points to the root of Google Test, relative to where this file is.
# MODIFY THIS FOR CME212
GTEST_DIR = /Users/sashayousefi/Desktop/cme212-sashayousefi/googletest/googletest

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# Flags passed to the C++ compiler for tests
GTEST_CXXFLAGS += -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP

# Define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX
ifeq ($(kernel-name), Linux)
  GTEST_LDLIBS += -lsfml-graphics -lsfml-window -lsfml-system -lX11 -lGL -lpthread
endif
ifeq ($(kernel-name), Darwin)
  GTEST_LDLIBS += -L/usr/local/opt/sfml/lib
  GTEST_LDLIBS += -L/usr/local/opt/libomp/lib
  GTEST_LDLIBS += -lsfml-graphics -lsfml-window -lsfml-system -framework OpenGL
endif

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^


################
## Test Rules ##
################

#HW0

gtest_hw0.o : gtest_hw0.cpp Graph.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) -c $<

gtest_hw0 : gtest_hw0.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $^ -o $@ $(GTEST_LDLIBS)

#HW1

gtest_hw1.o : gtest_hw1.cpp Graph.hpp shortest_path.hpp subgraph.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) -c $<

gtest_hw1 : gtest_hw1.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $^ -o $@ $(GTEST_LDLIBS)

#HW2

gtest_hw2_nodes.o : gtest_hw2_nodes.cpp Graph.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) -c $<

gtest_hw2_nodes : gtest_hw2_nodes.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $^ -o $@ $(GTEST_LDLIBS)

gtest_hw2_edges.o : gtest_hw2_edges.cpp Graph.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) -c $<

gtest_hw2_edges : gtest_hw2_edges.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $^ -o $@ $(GTEST_LDLIBS)

gtest_hw2_mass_spring.o : gtest_hw2_mass_spring.cpp Graph.hpp mass_spring.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) -c $<

gtest_hw2_mass_spring : gtest_hw2_mass_spring.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $^ -o $@ $(GTEST_LDLIBS)

#HW3

gtest_hw3.o : gtest_hw3.cpp mtl_test.cpp Graph.hpp IdentityMatrix.hpp GraphSymmetricMatrix.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $(INCLUDES) -c $<

gtest_hw3 : gtest_hw3.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $^ -o $@ $(GTEST_LDLIBS)

#HW4

gtest_hw4_testomp.o : gtest_hw4_testomp.cpp Graph.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $(INCLUDES) -c $< 

gtest_hw4_testomp : gtest_hw4_testomp.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $^ -o $@ $(GTEST_LDLIBS) $(LDFLAGS)

gtest_hw4.o : gtest_hw4.cpp Graph.hpp SpaceSearcher.hpp MortonCoder.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $(INCLUDES) -c $<

gtest_hw4 : gtest_hw4.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(GTEST_CXXFLAGS) $^ -o $@ $(GTEST_LDLIBS) $(LDFLAGS)
