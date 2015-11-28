#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
	double *ptr, *ptrRow, *ptrCol, temp;

    // Standard Gaussian Elimination
    for (k = 0; k < n - 1; k++) {
    	// Pivoting necessary (just return failure)
    	if (A[k][k] == 0) return -1;

    	L[k][k] = 1.0;
    	for (i = k + 1; i < n; i++)
    		l[i][k] = A[i][k] / A[k][k];              // Compute current column of L
    		
    	for (j = k + 1; j < n; j++) {
    		for (i = k + 1; i < n; i++) {
    			A[i][j] -= l[i][k] * A[k][j];         // Update submatrix of A
    		}
    	}
    }

    // Success
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

int main() {
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

	return 0;
}
