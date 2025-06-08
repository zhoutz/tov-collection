executables := inte_p inte_r inte_h inte_xi inte_uvy_p
headers := $(wildcard src/inc/*)

ifneq (command line,$(origin CXX))
  CXX := clang++
endif

CXXFLAGS := -std=c++20 -Isrc/inc -O3 -march=native -ffast-math

all: $(executables)

$(executables): % : src/%.cpp $(headers) build
	$(CXX) $(CXXFLAGS) $< -o build/$@

build:
	mkdir -p build
