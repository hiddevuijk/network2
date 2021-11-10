
test.exe: main.cpp network.h
	g++ main.cpp -std=c++11 -Wall -Werror  -o test.exe -lgsl -g


.PHONY : clean
clean:
	rm test.exe

