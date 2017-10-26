# =====================================================================
# RadMC build system
# =====================================================================
#
# If a Makefile.in exists in this directory, then use it.
#
-include Makefile.in


#
# Any macros that are omitted receive these default values:
AR       ?= ar rcu
RANLIB   ?= ranlib
CXX      ?= c++
CXXFLAGS ?= -std=c++11 -Wall -O3


# Build macros
# =====================================================================
PYSRC    := src/PythonWrapper.cpp
PYOBJ    := $(PYSRC:%.cpp=%.o)
PYDEP    := $(PYSRC:%.cpp=%.d)
SRC      := $(filter-out $(PYSRC), $(wildcard src/*.cpp))
OBJ      := $(SRC:%.cpp=%.o)
DEP      := $(SRC:%.cpp=%.d)
CXXFLAGS += -MMD -MP


# Build rules
# =====================================================================
#
default: radmc radmc.so

radmc: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

src/PythonWrapper.o: src/PythonWrapper.cpp	
	$(CXX) $^ -o $@ -c $(CXXFLAGS) -DRADMC_PYTHON $(shell python3.6-config --cflags)

radmc.so : $(OBJ) $(PYOBJ)
	$(CXX) $^ -o $@ -shared $(shell python3.6-config --ldflags)

clean:
	$(RM) $(OBJ) $(DEP) $(PYOBJ) $(PYDEP) radmc radmc.so

show:
	@echo $(SRC)

-include $(DEP)
