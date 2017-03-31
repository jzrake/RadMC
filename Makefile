

SRC := $(filter-out main.cpp test.cpp, $(wildcard *.cpp))
HDR := $(wildcard *.hpp)
OBJ := $(SRC:.cpp=.o)
CFLAGS = -Wall -O3


%.o : %.cpp $(HDR)
	$(CXX) $(CFLAGS) -o $@ -c -std=c++11 $<

main : main.o $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

test : test.o $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

show :
	@echo $(SRC)
	@echo $(OBJ)

clean :
	$(RM) $(OBJ) main.o test.o main test
