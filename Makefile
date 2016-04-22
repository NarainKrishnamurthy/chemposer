## Variables

SRC := src
SOURCES := $(wildcard src/*.cpp)
OBJS    := $(patsubst src/%,build/%,$(SOURCES:.cpp=.o))
CC  = clang-omp
CFLAGS = -Wall -Wextra -O2 -std=c++11 -fopenmp
EXEC_NAME = chemposer


## Targets

.PHONY: dirs clean

all: $(EXEC_NAME)

$(EXEC_NAME): dirs $(OBJS)
	$(CC) $(CFLAGS) -o $(EXEC_NAME) $(OBJS) -lc++

build/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
			
dirs:
	mkdir -p build

clean:
	rm -rf $(EXEC_NAME) build



 
