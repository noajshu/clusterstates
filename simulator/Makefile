DEPS=blossom *.cpp

all: blossom simulate
# g++  libs/blossom5-v2.05.src/MinCost/*.o libs/blossom5-v2.05.src/GEOM/*.o libs/blossom5-v2.05.src/*.o -o cubic cubic.cpp


blossom:
	cd libs/blossom5-v2.05.src/ && module load gcc/10.1.0 || true && make


simulate: $(DEPS)
	module load gcc/10.1.0 || true && \
		g++ -std=c++11 -I libs/cxxopts/ -I libs/blossom5-v2.05.src/ \
		libs/blossom5-v2.05.src/PMduals.o libs/blossom5-v2.05.src/PMexpand.o  \
		libs/blossom5-v2.05.src/PMinit.o libs/blossom5-v2.05.src/PMinterface.o \
		libs/blossom5-v2.05.src/PMmain.o libs/blossom5-v2.05.src/PMrepair.o \
		libs/blossom5-v2.05.src/PMshrink.o libs/blossom5-v2.05.src/misc.o \
		libs/blossom5-v2.05.src/MinCost/MinCost.o libs/blossom5-v2.05.src/GEOM/GPMinit.o \
		libs/blossom5-v2.05.src/GEOM/GPMinterface.o libs/blossom5-v2.05.src/GEOM/GPMkdtree.o \
		libs/blossom5-v2.05.src/GEOM/GPMmain.o -O3 \
		functions.cpp noise_models.cpp simulate.cpp -o simulate
