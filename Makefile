# CXX=clang++
CXX=g++
CXXFLAGS=-O3 -DNDEBUG -std=c++11 -I/usr/include -g
tests=test1.exe test2.exe test3.exe
.PHONY: clean run all

all: $(tests)

$(tests): %.exe: %.cc simultaneousDiagonalization.h \
	simultaneousDiagonalization.o
	$(CXX) -o $@ $(CXXFLAGS) $< simultaneousDiagonalization.o	

simultaneousDiagonalization.o: simultaneousDiagonalization.h

.cxx.o: 
	$(CXX) -c $(CXXFLAGS) $<
.cc.o: 
	$(CXX) -c $(CXXFLAGS) $<


clean:
	rm -rf *.o *.exe *.txt

run: $(tests)
	./test1.exe | tee test1.txt
	./test2.exe | tee test2.txt
	./test3.exe | tee test3.txt

