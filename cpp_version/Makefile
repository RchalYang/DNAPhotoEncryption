all: dna chaotic_system reduction transition rotate dna_io_reduce reduction_io_reduce rotate_io_reduce
	nvcc -o out1 test_opencv.cpp ./lib/reduce.o ./lib/rotate.o ./lib/dna.o ./lib/transition.o ./lib/chaotic_system.o -std=c++11 `pkg-config --cflags --libs opencv` --gpu-architecture compute_60
	nvcc -o out2 reduce_IO.cu ./lib_io_reduce/reduce_io_reduce.o ./lib_io_reduce/rotate_io_reduce.o ./lib_io_reduce/dna_io_reduce.o ./lib/transition.o ./lib/chaotic_system.o -std=c++11 `pkg-config --cflags --libs opencv` --gpu-architecture compute_60
	g++ -std=c++11 -o outS SequencialVersion.cpp ./lib/chaotic_system.o ./lib/transition.o `pkg-config --cflags --libs opencv`

chaotic_system:
	g++ -std=c++11 -o ./lib/chaotic_system.o -c ./lib/chaotic_system.cpp 

reduction: dna
	nvcc -std=c++11 -o ./lib/reduce.o -c ./lib/reduce.cu --gpu-architecture compute_60

dna_io_reduce:
	nvcc  -o ./lib_io_reduce/dna_io_reduce.o -c ./lib_io_reduce/dna_io_reduce.cu -std=c++11 --gpu-architecture compute_60
	
rotate_io_reduce:
	nvcc  -o ./lib_io_reduce/rotate_io_reduce.o -c ./lib_io_reduce/rotate_io_reduce.cu -std=c++11 --gpu-architecture compute_60

reduction_io_reduce: dna_io_reduce
	nvcc -std=c++11 -o ./lib_io_reduce/reduce_io_reduce.o -c ./lib_io_reduce/reduce_io_reduce.cu --gpu-architecture compute_60

dna:
	nvcc  -o ./lib/dna.o -c ./lib/dna.cu -std=c++11 --gpu-architecture compute_60
	
rotate:
	nvcc  -o ./lib/rotate.o -c ./lib/rotate.cu -std=c++11 --gpu-architecture compute_60


transition:
	g++ -std=c++11 -o ./lib/transition.o -c ./lib/transition.cpp

clean:
	rm ./lib/*.o
	rm ./out
	# rm .