SRC := $(filter-out src/main.cpp src/test.cpp, $(wildcard src/*.cpp))
HDR := $(wildcard src/*.hpp)
OBJ := $(SRC:.cpp=.o)
CFLAGS = -std=c++11 -Wall -O3

%.o : %.cpp $(HDR)
	$(CXX) $(CFLAGS) -o $@ -c $<

radmc : src/main.o $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

test : src/test.o $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

show :
	@echo $(SRC)
	@echo $(OBJ)

clean :
	$(RM) $(OBJ) src/main.o src/test.o radmc test
