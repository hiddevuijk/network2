
test.exe: main.cpp main2.cpp network.h minimize.h minimize_gsl.h fire2.h get_gamma_list.h
	g++ main.cpp -std=c++11 -Wall -Werror  -o test.exe -lgsl -O3



.PHONY : clean
clean:
	rm test.exe

