build: compile link

all: build run

compile:
	g++ -c main.cpp

link:
	g++ -o main.exe main.o -std=c++11 -lsuperlu -lm -lblas -lstdc++ -larmadillo

run:
	./main.exe 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.2 0 1e10 1 n
	./main.exe 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.2 0 1e10 2 n
	./main.exe 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.2 0 1e10 3 n
	python3 plotter.py
