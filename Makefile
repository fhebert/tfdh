
# Makefile for a simple, single-exec, C++ build
# * expects all source files (cpp,h) in ./src
# * produces object and dependency files in ./build
# * produces executable in ./bin

EXECUTABLE := tfdh

CXX := clang++
CPPFLAGS :=
CXXFLAGS := -O3 -Wall -Wextra -std=c++11 -march=native
LIBS := -lm -lgsl

SRCS := $(wildcard src/*.cpp)
OBJS := $(subst src,build,$(SRCS:.cpp=.o))
DEPS := $(subst src,build,$(SRCS:.cpp=.d))

.PHONY: all
all: bin/$(EXECUTABLE)

bin/$(EXECUTABLE): $(OBJS) | bin
	@ echo "  CXXLD     $(EXECUTABLE)"
	@ $(CXX) $(CXXFLAGS) $(LIBS) $(OBJS) -o bin/$(EXECUTABLE)

build/%.o: src/%.cpp | build
	@ echo "  CXX       $*.cpp"
	@ $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c src/$*.cpp -o build/$*.o -MMD

build:
	@ mkdir -p build

bin:
	@ mkdir -p bin

.PHONY: clean
clean:
	@ $(RM) $(OBJS) $(DEPS) bin/$(EXECUTABLE)

.PHONY: immaculate
immaculate: clean
	@ $(RM) -r build bin

-include $(DEPS)

