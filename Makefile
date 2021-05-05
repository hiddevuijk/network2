
test.exe: main.cpp graph.h  vec2.h generate_graph.h
	g++ main.cpp -std=c++11 -Wall -Werror  -o test.exe  -g


