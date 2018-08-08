COMPILER=g++-8
HEADERS=Mesh.hpp
FILE=main.cpp
EXE=./a.out
REMOVER=./remover.sh

all:
	$(COMPILER) $(HEADERS) $(FILE)
	$(EXE)
	$(REMOVER)