# =====================================================================
# RadMC build system
# =====================================================================
#
# If a Makefile.in exists in this directory, then use it.
#
-include Makefile.in


#
# Any macros that are omitted receive these default values:
LUA_ARCH ?= generic
AR       ?= ar rcu
RANLIB   ?= ranlib
CXX      ?= c++
CXXFLAGS ?= -std=c++11 -Wall -O3


# Build macros
# =====================================================================
SRC      := $(wildcard src/*.cpp)
OBJ      := $(SRC:%.cpp=%.o)
DEP      := $(SRC:%.cpp=%.d)
CXXFLAGS += -MMD -MP


# Build rules
# =====================================================================
#
radmc: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

radmc.so : src/PythonWrapper.cpp $(OBJ)
	$(CXX) $^ -o $@ -shared -std=c++11 -DRADMC_PYTHON $(shell python-config --cflags --ldflags)

clean:
	$(RM) $(OBJ) $(DEP) radmc radmc.so

-include $(DEP)
