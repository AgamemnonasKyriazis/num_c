#ifndef __NUM_C__
#define __NUM_C__
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
/*
**  formula BN = floor[ sqrt( CACHE_SIZE / (N * sizeof(double) ) ) ]
*/
#define B2 45 //Max N for 2 NxN matrices
#define B3 36 //Max N for 3 NxN matrices

#define MAX_RAND  32000

struct matrix_s {
    double **matrix;
    int rows;
    int collumns;
};
typedef struct matrix_s matrix_t;

void init();
matrix_t *allocatem(int rows, int collumns);
void releasem(matrix_t *m_t);
matrix_t *copy(matrix_t *m_t1);
matrix_t *get_random(int rows, int collumns);
matrix_t *get_random_norm(int rows, int collumns);
void printm(matrix_t *m_t);

matrix_t *zeros(int rows, int collumns);
matrix_t *ones(int rows, int collumns);
matrix_t *identity(int rows, int collumns);

matrix_t *transpose(matrix_t *m_t);
matrix_t *transpose2(matrix_t *m_t);

matrix_t *multiply(matrix_t *m_t1, matrix_t *m_t2);
matrix_t *multiply2(matrix_t *m_t1, matrix_t *m_t2);
matrix_t *multiplys(matrix_t *m_t1, double scalar);

matrix_t *add(matrix_t *m_t1, matrix_t *m_t2);
matrix_t *add2(matrix_t *m_t1, matrix_t *m_t2);
matrix_t *adds(matrix_t *m_t1, double scalar);

matrix_t *sub(matrix_t *m_t1, matrix_t *m_t2);
matrix_t *sub2(matrix_t *m_t1, matrix_t *m_t2);
matrix_t *subs(matrix_t *m_t1, double scalar);

double det(matrix_t *m_t1);

matrix_t *inverse(matrix_t *m_t1);

#endif