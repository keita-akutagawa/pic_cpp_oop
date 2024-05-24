CC = g++
CFLAGS = -std=c++20 -Wall -pedantic-errors -fopenmp -O3 -mtune=native -march=native
SRCDIR = ../lib_pic1d_cpp
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(SRCDIR)/%.o, $(SOURCES))
MAIN = main_two_stream_electron.o
EXECUTABLE = program

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS) $(MAIN)
	$(CC) $(CFLAGS) $(OBJS) $(MAIN) -o $(EXECUTABLE)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(MAIN): main_two_stream_electron.cpp
	$(CC) $(CFLAGS) -c $< -o $@

run: $(EXECUTABLE)
	./$(EXECUTABLE)

clean:
	rm -f $(OBJS) $(MAIN) $(EXECUTABLE)
