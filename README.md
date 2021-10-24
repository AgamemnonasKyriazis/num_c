# NUM_C LINEAR ALGEBRA LIB
This is a library that implements linear algebra functions in C.

## Description
  This is a premature low-tested version of what I hope to be the final result.
This library contains the most common matrix operations found in linear algebra.
Also there are some semi-optimized implementations of cache heavy functions such as matrix multiplication and transposition.
These methods are designed to use 16KB of cache memory maximum in order to be able to run in slow or cache-starved systems.
I will try to implement and optimize the library in the future but for now I think it gets the job done for anyone that would like to use it.

## About
  This library uses a datatype called matrix_t. This datatype contains both the data and the metadata of the matrix (rows and collumns).

### List of implemented functions

* matrix_t *allocatem(int rows, int collumns);
Creates a new matrix_t by allocating memory and setting the size of rows and collumns

* void releasem(matrix_t *m_t);
Calls free iteratively to return the allocated memory back to the system

* matrix_t *copy(matrix_t *m_t1);
Copies the data of matrix m_t1 to a new matrix that it returns

* matrix_t *get_random(int rows, int collumns);
Returns a matrix with set dimensions but with random elements
* matrix_t *get_random_norm(int rows, int collumns);
Returns a matrix with set dimensions but with random elements with range from -1 to 1
* void init();
Is needed once before requesting a random matrix and never again!

* void printm(matrix_t *m_t);
Prints the matrix m_t1

* matrix_t *zeros(int rows, int collumns);
Returns a matrix whose all elements are zeroes
* matrix_t *ones(int rows, int collumns);
Returns a matrix whose all elements are ones
* matrix_t *identity(int rows, int collumns);
Returns the identity matrix with set dimensions

* matrix_t *transpose(matrix_t *m_t);
* matrix_t *transpose2(matrix_t *m_t);
Both traspose a matrix, transpose2 is unoptimized

* matrix_t *multiply(matrix_t *m_t1, matrix_t *m_t2);
* matrix_t *multiply2(matrix_t *m_t1, matrix_t *m_t2);
Both multiply two matrices but multiply2 is unoptimized
* matrix_t *multiplys(matrix_t *m_t1, double scalar);
Multiplies a matrix with a scalar

* matrix_t *add(matrix_t *m_t1, matrix_t *m_t2);
* matrix_t *add2(matrix_t *m_t1, matrix_t *m_t2);
Both add two matrices but add2 is unoptimized
* matrix_t *adds(matrix_t *m_t1, double scalar);
Adds a scalar to a matrix

* matrix_t *sub(matrix_t *m_t1, matrix_t *m_t2);
* matrix_t *sub2(matrix_t *m_t1, matrix_t *m_t2);
Both subtract two matrices but sub2 is unoptimized
* matrix_t *subs(matrix_t *m_t1, double scalar);
Subtracts a scalar from a matrix

* double det(matrix_t *m_t1);
Calculates the determinant of a matrix !unoptimized!

* matrix_t *inverse(matrix_t *m_t1);
Calculates the inverse of a matrix !unoptimized!
