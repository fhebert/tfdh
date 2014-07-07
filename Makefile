
CXX = clang++
CPPFLAGS = -DDEBUG
CXXFLAGS = -O3 -Wall -Wextra -std=c++11 -march=native

SRCS = $(wildcard src/*.cpp)
OBJS = $(subst src,build,$(SRCS:.cpp=.o))
DEPS = $(subst src,build,$(SRCS:.cpp=.d))


all: exec

exec: $(OBJS) | bin
	@ echo "  CXXLD     tfdh"
	@ $(CXX) $(CXXFLAGS) $(OBJS) -o bin/tfdh -lgsl

build/%.o: src/%.cpp | build
	@ echo "  CXX       $*.cpp"
	@ $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c src/$*.cpp -o build/$*.o
	@ $(CXX) $(CPPFLAGS) src/$*.cpp -MM -MF build/$*.d.preprocess
	@ sed -e 's|.*:|build/$*.o:|' < build/$*.d.preprocess > build/$*.d
	@ $(RM) build/$*.d.preprocess

build:
	@ mkdir -p build

bin:
	@ mkdir -p bin

clean:
	@ $(RM) $(OBJS) $(DEPS) bin/tfdh

immaculate: clean
	@ $(RM) -r build bin

-include $(DEPS)

