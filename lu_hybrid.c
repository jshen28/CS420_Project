#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mpi.h"

// For now use 2-dimension matrix, try array to see if more efficient
typedef double **mat;

// Fill the matrix with 0
void initZero(mat x, int n) {
    int i,j;
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			x[i][j] = 0;
}

// Initialize matrix and fill with 0
mat initNew(int n) {
	mat x = malloc(sizeof(double*) * n);
	x[0] = malloc(sizeof(double) * n * n);

    int i;
	for (i = 0; i < n; i++)
		x[i] = x[0] + n * i;
	initZero(x,n);

	return x;
}

// Clean up
void matDel(mat x) {
	free(x[0]);
	free(x);
}

void free2DArr(double** m, int row, int col){
	int i;
	for(i = 0; i < row; i++){
		free(m[i]);
	}
	free(m);
}

// Copy matrix to proper type
mat matCopy(void *s, int n) {
	mat x = initNew(n);
	int i,j;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			x[i][j] = ((double(*)[n])s)[i][j];
	return x;
}

// Actual Decompsition
// Now, matrix L should have all 1 on diagonal entries. And U is the 
// upper triangular matrix of A after algorithm.
int LUDecomp(mat A, mat L, int n) {
	// Declarations
	int i, j, k, p;
	double *ptr, *ptrRow, *ptrCol;

    // Standard Gaussian Elimination
    for (k = 0; k < n - 1; k++) {
    	// Pivoting necessary (just return failure)
    	if (A[k][k] == 0) return -1;

    	L[k][k] = 1.0;
    	for (i = k + 1; i < n; i++)
    		L[i][k] = A[i][k] / A[k][k];              // Compute current column of L
    		
    	for (j = k + 1; j < n; j++) {
    		for (i = k + 1; i < n; i++) {
    			A[i][j] -= L[i][k] * A[k][j];         // Update submatrix of A
    		}
    	}
    }

    // Success
    return 1;
}

void randMatrixGene(double*** m, int row, int col, int seed){
	srand((unsigned) time(NULL) + (unsigned) (seed * seed));
	int i, j;
	for (i = 0; i < row; ++i){
		for (j = 0; j < col; ++j){
			(*m)[i][j] =((double)rand()/(double)RAND_MAX);
		}
	}
}

void initZero2(double** L, int row, int col){
	int i, j;
	for(i = 0; i < row; i++){
		for(j = 0; j < col; j++){
			L[i][j] = 0;
		}
	}
}

void initL(double*** L, int row, int col, int rank){
	int i, ii, j;
	for(i = 0; i < row; i++){
		ii = i + rank * row;
		for(j = 0; j < col; j++){
			if( ii == j){
				(*L)[i][j] = 1.0;
			}
			else{
				(*L)[i][j] = 0.0;
			}
		}
	}
	
}

double init2DArray(double*** m, int row, int col){
	*m = (double**) malloc(sizeof(double*) * row);
	int i;
	for(i = 0; i < row; i++){
		(*m)[i] = (double*) malloc(sizeof(double) * col);
	}
}

void MD2SD(double*** data, double** ldata, int row, int col){
	int i, j, count;
	double a;
	count = 0;
	for(i = 0; i < row; i++){
		for(j = 0; j < col; j++){
			a = (*data)[i][j];		
			(*ldata)[count++] = a;
		}
	}
}

void showData(double* data, int size){
	int i = 0;
	for(i = 0; i < size; i++){
		printf("%.4lf ", data[i]);
	}
	printf("\n");
}


void matrixPro(double *a, double *b, int row, int col){

	double *c;
	c = (double*) malloc(sizeof(double) * row * row);

	int i, j, k, l;
	for(k = 0; k< row; k++){
		for(l = 0; l < row; l++){
			for(i = 0; i < col; i++){
				c[k*col+l] += a[k*col+i] * b[i*row+l];
			}
		}
	}
	printf("Product of L and U\n");
	for(i = 0; i < row; i++){
		for(j = 0; j< row; j++){
			printf("%lf ", c[i*row+j]);
		}
		printf("\n");
	}

	free(c);
}

// Simple row decomposition version;
int LUDecomp_RD(int* argc, char*** argv){

	
	// Declarations
	int i, j, k, p, rank, procs, n, size;
	double **l, **data, *ptrRow;
	double *ldata, *ll, *wdata;

	// Initialization
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);


	if(*argc <= 1){
		if(rank == 0)
			printf("Not enough input you stupid ass!\n");
		MPI_Finalize();
		return -1;
	}
	size = atoi((*argv)[1]);
	
	// Initialized a random matrix
	n = size / procs; // suppose divisible;

	init2DArray(&data, n, size);
	randMatrixGene(&data, n, size, rank);

	init2DArray(&l, n, size);
	initL(&l, n, size, rank);
	
    wdata = (double*) malloc(sizeof(double) * size * size);
    ldata = (double*) malloc(sizeof(double) * n * size);
   	MD2SD(&data, &ldata, n, size);
   	
   	MPI_Gather(ldata, n*size, MPI_DOUBLE, wdata, n*size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   	
	// Parallelized version
	for(k = 0; k < size - 1; k++){
		ptrRow = (double*) malloc(sizeof(double) * (size - k) );
		
		if( k/n == rank){
			p = k%n;
			for(i = k; i < size; i++){
				ptrRow[i-k] = data[p][i]; // spread this information 
			}
		}
		MPI_Bcast(ptrRow, size-k, MPI_DOUBLE, k/n, MPI_COMM_WORLD); // collective call;
		if( fabs((ptrRow[0])) < 0.0001 ){
			MPI_Finalize();
			return -1;
		}
		
		if(k >= n * (rank+1)) continue;

		int cc = (k/n == rank)? (k%n)+1 : 0;

		#pragma omp parallel for
		for(i = cc; i < n; i++){
			l[i][k] = data[i][k] / ptrRow[0];

			for(j = k; j < size; j++){
				data[i][j] -= l[i][k] * ptrRow[j-k];
			}
		}

		free(ptrRow);
	}
	
	// send data back to process 0;
	// MPI_Gather()
	double *U, *L;
	if(rank == 0){
		U = (double*) malloc(sizeof(double) * size * size);
		L = (double*) malloc(sizeof(double) * size * size);
	}

	
	MD2SD(&data, &ldata, n, size);
	ll = (double*) malloc(sizeof(double) * n * size);
	MD2SD(&l, &ll, n, size);

	MPI_Gather(ldata, size*n, MPI_DOUBLE, U, size*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(ll, size*n, MPI_DOUBLE, L, size*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//free2DArr(data, n, size);
	//free2DArr(l, n, size);
	
	if(rank == 0){
		printf("Whole data:\n");
		for(i = 0; i < size; i++){
			for(j = 0; j< size; j++){
				printf("%lf ", wdata[i*size+j]);
			}
			printf("\n");
		}
		printf("U:\n");
		for(i = 0; i < size; i++){
			for(j = 0; j< size; j++){
				printf("%lf ", U[i*size+j]);
			}
			printf("\n");
		}
		
		printf("L:\n");
		for(i = 0; i < size; i++){
			for(j = 0; j< size; j++){
				printf("%lf ", L[i*size+j]);
			}
			printf("\n");
		}

		matrixPro(L, U, size, size);

		free(U); free(L);
	}

	free(ll); free(ldata);
	
	MPI_Finalize();
	
	return 1;
}

void printMatrix(mat x, int n) {
    int i,j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%8.4g", x[i][j]);
			printf(j < n - 1 ? "  " : i == n - 1 ? "\n\n" : "\n");
		}
	}
}

int main(int argc, char** argv) {


	LUDecomp_RD(&argc, &argv);
	/*
	double A1[][3] = {{8,2,9}, {4,9,4}, {6,7,9}};

    mat A = matCopy(A1,3);

	mat L1 = initNew(3);
	int res = LUDecomp(A,L1,3);
	if (res == 1) {
		printMatrix(A,3);
		printMatrix(L1,3);
	}
	else
        printf("Need pivoting.");

    matDel(A);
    matDel(L1);
	*/
	return 0;
}
