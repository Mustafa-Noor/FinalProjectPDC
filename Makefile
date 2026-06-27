CC := gcc
MPICC := mpicc
CFLAGS := -std=c99 -Wall -Wextra -pedantic -O0 -g

SEQ_SRC := src/sequential/main_seq.c
SEQ_BIN := bin/main_seq
PTH_SRC := src/pthreads/main_pth.c
PTH_BIN := bin/main_pth
OMP_SRC := src/openmp/main_omp.c
OMP_BIN := bin/main_omp
MPI_SRC := src/mpi/main_mpi.c
MPI_BIN := bin/main_mpi
GEN_SRC := data/input_generator.c
GEN_BIN := bin/input_generator

.PHONY: all sequential pthreads openmp mpi input_generator clean

all: sequential pthreads openmp mpi input_generator

sequential: $(SEQ_BIN)

pthreads: $(PTH_BIN)

openmp: $(OMP_BIN)

mpi: $(MPI_BIN)

input_generator: $(GEN_BIN)

$(SEQ_BIN): $(SEQ_SRC)
	mkdir -p bin
	$(CC) $(CFLAGS) $(SEQ_SRC) -o $(SEQ_BIN)

$(PTH_BIN): $(PTH_SRC)
	mkdir -p bin
	$(CC) $(CFLAGS) $(PTH_SRC) -o $(PTH_BIN) -pthread

$(OMP_BIN): $(OMP_SRC)
	mkdir -p bin
	$(CC) $(CFLAGS) $(OMP_SRC) -o $(OMP_BIN) -fopenmp

$(MPI_BIN): $(MPI_SRC)
	mkdir -p bin
	$(MPICC) $(CFLAGS) $(MPI_SRC) -o $(MPI_BIN)

$(GEN_BIN): $(GEN_SRC)
	mkdir -p bin
	$(CC) $(CFLAGS) $(GEN_SRC) -o $(GEN_BIN)

clean:
	rm -rf bin
