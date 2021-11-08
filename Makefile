
test.exe: main.cpp graph.h  vec2.h generate_graph.h network.h fire.h fire2.h ConfigFile.h minimize.h
	g++ main.cpp -std=c++11 -Wall -Werror  -o test.exe  -O3

.PHONY : clean
clean:
	rm test.exe

