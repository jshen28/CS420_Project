#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv){

	int rank, procs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	int *data, *wdata;
	
	data = (int*) malloc(sizeof(int) * 3);

	int i, j;

	for(i = 0; i < 3; i++)
		data[i] = rank * (i+1);
	
	if(rank == 0){
		wdata = (int*) malloc(sizeof(int) * procs * 3);
	} 
	
	MPI_Gather(data, 3, MPI_INT, wdata, 3, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank == 0){
		for(i = 0; i < procs; i++){
			for(j = 0; j < 3; j++){
				printf("%d ", wdata[3*i+j]);
			}
			printf("\n");
		}
	}
	MPI_Finalize();
	return 0;
}