## Variables

SRC := src
SOURCES := $(wildcard src/*.cpp)
OBJS    := $(patsubst src/%,build/%,$(SOURCES:.cpp=.o))
CC  = g++
CFLAGS = -Wall -Wextra -O2 -std=c++11
EXEC_NAME = chemposer


## Targets

all: $(EXEC_NAME)

$(EXEC_NAME): $(SOURCES) $(OBJS)
	$(CC) $(CFLAGS) -o $(EXEC_NAME) $(OBJS)

build/%.o: $(SRC)/%.cpp build
	$(CC) $(CFLAGS) -c $< -o $@

build:
	mkdir -p $@

clean:
	rm -rf $(EXEC_NAME) build



 
