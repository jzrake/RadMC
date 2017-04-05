

SRC := $(filter-out src/main.cpp src/test.cpp, $(wildcard src/*.cpp))
HDR := $(wildcard src/*.hpp)
OBJ := $(SRC:.cpp=.o)
CFLAGS = -Wall -O0 -g


%.o : %.cpp $(HDR)
	$(CXX) $(CFLAGS) -o $@ -c -std=c++11 $<

radmc : src/main.o $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

test : src/test.o $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

show :
	@echo $(SRC)
	@echo $(OBJ)

clean :
	$(RM) $(OBJ) src/main.o src/test.o radmc test
