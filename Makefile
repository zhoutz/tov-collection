executables := inte_p inte_r inte_h
headers := $(wildcard src/inc/*)

CXX = clang++
CXXFLAGS = -std=c++20 -Isrc/inc -O3 -march=native -ffast-math


all: $(executables)

$(executables): % : src/%.cpp $(headers) build
	$(CXX) $(CXXFLAGS) $< -o build/$@

build:
	mkdir -p build
