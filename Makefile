# FLAGS=-g -Wall # debug
# FLAGS=-O2 -g -pg -Wall # profile
FLAGS=-O3 -Wall # release

main:
	g++ $(FLAGS) -std=c++0x -o afkm main.cpp vfkm/Vector.cpp vfkm/PolygonalPath.cpp vfkm/Vector2D.cc vfkm/Grid.cpp vfkm/Optimizer.cpp vfkm/ConstraintMatrix.cpp attributetrajectory.cpp util.cpp randomsampler.cpp -I.

clean:
	rm afkm
	rm *.o
