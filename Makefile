
SRC := src
SOURCES := $(wildcard src/*.cpp)
NVSOURCES := $(wildcard src/*.cu)
OBJS    := $(patsubst src/%,build/%,$(SOURCES:.cpp=.o))
NVOBJS    := $(patsubst src/%,nvbuild/%,$(NVSOURCES:.cu=.o))
NCC = nvcc
GCC = g++
NVFLAGS = -std=c++11 -O2
CFLAGS = -Wall -Wextra -O2 -std=c++11
EXEC_NAME = chemposer


## Targets

.PHONY: dirs clean

all: $(EXEC_NAME)

$(EXEC_NAME): dirs $(OBJS) $(NVOBJS)
	$(NCC) $(NVFLAGS) -o $(EXEC_NAME) $(OBJS) $(NVOBJS)

build/%.o: $(SRC)/%.cpp
	$(GCC) $(CFLAGS) -c $< -o $@

nvbuild/%.o: $(SRC)/%.cu
	$(NCC) $(NVFLAGS) -c $< -o $@


dirs:
	mkdir -p build nvbuild

clean:
	rm -rf $(EXEC_NAME) build nvbuild

