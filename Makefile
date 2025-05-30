executables := inte_p inte_r
headers := $(wildcard src/inc/*)

CXX ?= clang++
CXXFLAGS = -std=c++20 -O3 -march=native -Isrc/inc


all: $(executables)

$(executables): % : src/%.cpp $(headers) build
	$(CXX) $(CXXFLAGS) $< -o build/$@

build:
	mkdir -p build
