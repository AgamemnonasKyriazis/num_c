#include "num_c.h"

void init() {
    srand (time(NULL));
}

matrix_t *allocatem(int rows, int collumns) {
    int i;
    matrix_t *m_t = (matrix_t *)malloc(sizeof(matrix_t));
    m_t->matrix = (double **)malloc(rows*sizeof(double *));
    m_t->rows = rows;
    m_t->collumns = collumns;
    for(i = 0; i < rows; i++)
        m_t->matrix[i] = (double *)malloc(collumns*sizeof(double));
    return m_t;
}

void releasem(matrix_t *m_t) {
    int i = m_t->rows;
    while(i) {
        free(m_t->matrix[--i]);
    }
    free(m_t);
}

matrix_t *copy(matrix_t *m_t1) {
    matrix_t *c = NULL;
    int i, j;
    c = allocatem(m_t1->rows, m_t1->collumns);
    for(i = 0; i < m_t1->rows; i = i+1) {
        for(j = 0; j < m_t1->collumns; j = j+1) {
            c->matrix[i][j] = m_t1->matrix[i][j];
        }
    }
    return c;
}

matrix_t *get_random(int rows, int collumns) {
    matrix_t *c = NULL;
    int i, j;
    c = allocatem(rows, collumns);
    for(i = 0; i < rows; i = i+1) {
        for(j = 0; j < collumns; j = j+1) {
            c->matrix[i][j] = (double)((rand()%MAX_RAND)/sin(rand()));
        }
    }
    return c;
}

matrix_t *get_random_norm(int rows, int collumns) {
    matrix_t *c = NULL;
    int i, j;
    c = allocatem(rows, collumns);
    for(i = 0; i < rows; i = i+1) {
        for(j = 0; j < collumns; j = j+1) {
            c->matrix[i][j] = (double)(sin(rand()));
        }
    }
    return c;
}

void printm(matrix_t *m_t) {
    int i, j;
    for(i = 0; i < m_t->rows; i = i+1) {
        for(j = 0; j < m_t->collumns; j = j+1) {
            printf("|%f|", m_t->matrix[i][j]);
        }
        printf("\n");
    }
}

matrix_t *zeros(int rows, int collumns) {
    matrix_t *m_t = NULL;
    int i, j;
    m_t = allocatem(rows, collumns);
    for(i = 0; i < m_t->rows; i = i+1) {
        for(j = 0; j < m_t->collumns; j = j+1) {
            m_t->matrix[i][j] = 0;
        }
    }
    return m_t;
}

matrix_t *ones(int rows, int collumns) {
    matrix_t *m_t = NULL;
    int i, j;
    m_t = allocatem(rows, collumns);
    for(i = 0; i < m_t->rows; i = i+1) {
        for(j = 0; j < m_t->collumns; j = j+1) {
            m_t->matrix[i][j] = 1;
        }
    }
    return m_t;
}

matrix_t *identity(int rows, int collumns) {
    matrix_t *m_t = NULL;
    int i, j;
    if(rows != collumns) {
        exit(1);
    }
    m_t = allocatem(rows, collumns);
    for(i = 0; i < m_t->rows; i = i+1) {
        for(j = 0; j < m_t->collumns; j = j+1) {
            m_t->matrix[i][j] = (i == j)? 1 : 0;
        }
    }
    return m_t;
}

matrix_t *transpose(matrix_t *m_t) {
    int i, j,  m, n;
    matrix_t *c = allocatem(m_t->collumns, m_t->rows);
    int Bi, Bj;
    int rm1 = m_t->rows, cm1 = m_t->collumns;


    Bi = (B2 < cm1)? B2 : cm1;
    Bj = (B2 < rm1)? B2 : rm1;

    for(i = 0; i < rm1; i = i+Bi) {
        Bi = (Bi < rm1-i)? Bi : rm1-i;
        for(j = 0; j < cm1; j = j+Bj) {
            Bj = (Bj < cm1-j)? Bj : cm1-j;
            for(m = 0; m < Bi; m++) {
                for(n = 0; n < Bj; n++) {
                    c->matrix[j+n][i+m] = m_t->matrix[i+m][j+n];
                }
            }
        }
    }
    return c;
}

matrix_t *transpose2(matrix_t *m_t) {
    int i, j;
    matrix_t *c = allocatem(m_t->collumns, m_t->rows);
    int rm1 = m_t->rows, cm1 = m_t->collumns;

    for(i = 0; i < rm1; i = i+1) {
        for(j = 0; j < cm1; j = j+1) {
            c->matrix[j][i] = m_t->matrix[i][j];
        }
    }
    return c; 
}

matrix_t *multiply(matrix_t *m_t1, matrix_t *m_t2) {
    matrix_t *c = NULL;
    double sum;
    int i, j, k, m, n, p;
    int Bi, Bj, Bk;
    int rm1 = m_t1->rows, cm1 = m_t1->collumns;
    int rm2 = m_t2->rows, cm2 = m_t2->collumns;

    if(cm1 != rm2)
        exit(1);

    c = allocatem(rm1, cm2);
    Bi = (B3 < rm1)? B3 : rm1;
    Bj = (B3 < cm2)? B3 : cm2;
    Bk = (B3 < cm1)? B3 : cm1;

    for(i = 0; i < rm1; i = i+Bi) {
        Bi = (Bi < rm1-i)? Bi : rm1-i;
        for(j = 0; j < cm2; j = j+Bj) {
            Bj = (Bj < cm2-j)? Bj : cm2-j;
            for(k = 0; k < cm1; k = k+Bk) {
                Bk = (Bk < cm1-k)? Bk : cm1-k;
                for(m = 0; m < Bi; m++) {
                    for(n = 0; n < Bj; n++) {
                        sum = (double)0;
                        for(p = 0; p < Bk; p++) {
                            sum += m_t1->matrix[i+m][k+p]*m_t2->matrix[k+p][j+n];
                        }
                        c->matrix[i+m][j+n] = sum;
                    }
                }
            }
        }
    }
    return c;
}

matrix_t *multiply2(matrix_t *m_t1, matrix_t *m_t2) {
    matrix_t *c = NULL;
    double sum;
    int i, j, k;
    int rm1 = m_t1->rows, cm1 = m_t1->collumns;
    int rm2 = m_t2->rows, cm2 = m_t2->collumns;

    if(cm1 != rm2)
        exit(1);

    c = allocatem(rm1, cm2);

    for(i = 0; i < rm1; i++) {
        for(j = 0; j < cm2; j++) {
           sum = (double)0;
           for(k = 0; k < cm1; k++) {
               sum += m_t1->matrix[i][k]*m_t2->matrix[k][j];
           }
           c->matrix[i][j] = sum;
        }
    }
    return c;
}

matrix_t *multiplys(matrix_t *m_t1, double scalar) {
    matrix_t *c = NULL;
    int i, j;
    int rows = m_t1->rows;
    int collumns = m_t1->collumns;

    c = allocatem(rows, collumns);

    for(i = 0; i < rows; i = i+1) {
        for(j = 0; j < collumns; j = j+1) {
            c->matrix[i][j] = m_t1->matrix[i][j] * scalar;
        }
    }
    return c;
}

matrix_t *add(matrix_t *m_t1, matrix_t *m_t2) {
    matrix_t *c = NULL;
    int i, j, m, n;
    int Bi, Bj;
    int rows = m_t1->rows;
    int collumns = m_t1->collumns;

    if((m_t1->rows != m_t2->rows)&&(m_t1->collumns != m_t2->collumns))
        exit(1);

    c = allocatem(rows, collumns);
    Bi = (B3 < rows)? B3 : rows;
    Bj = (B3 < collumns)? B3 : collumns;

    for(i = 0; i < rows; i = i+Bi) {
        Bi = (Bi < rows-i)? Bi : rows-i;
        for(j = 0; j < collumns; j = j+Bj) {
            Bj = (Bj < collumns-j)? Bj : collumns-j;
            for(m = 0; m < Bi; m++) {
                for(n = 0; n < Bj; n++) {
                    c->matrix[i+m][j+n] = m_t1->matrix[i+m][j+n] + m_t2->matrix[i+m][j+n];
                }
            }
        }
    }
    return c;
}

matrix_t *add2(matrix_t *m_t1, matrix_t *m_t2) {
    matrix_t *c = NULL;
    int i, j;
    int rows = m_t1->rows;
    int collumns = m_t1->collumns;

    if((m_t1->rows != m_t2->rows)&&(m_t1->collumns != m_t2->collumns))
        exit(1);

    c = allocatem(rows, collumns);

    for(i = 0; i < rows; i = i+1) {
        for(j = 0; j < collumns; j = j+1) {
            c->matrix[i][j] = m_t1->matrix[i][j] + m_t2->matrix[i][j];
        }
    }
    return c;
}

matrix_t *adds(matrix_t *m_t1, double scalar) {
    matrix_t *c = NULL;
    int i, j;
    int rows = m_t1->rows;
    int collumns = m_t1->collumns;

    c = allocatem(rows, collumns);

    for(i = 0; i < rows; i = i+1) {
        for(j = 0; j < collumns; j = j+1) {
            c->matrix[i][j] = m_t1->matrix[i][j] + scalar;
        }
    }
    return c;
}

matrix_t *sub(matrix_t *m_t1, matrix_t *m_t2) {
    matrix_t *c = NULL;
    int i, j, m, n;
    int Bi, Bj;
    int rows = m_t1->rows;
    int collumns = m_t1->collumns;

    if((m_t1->rows != m_t2->rows)&&(m_t1->collumns != m_t2->collumns))
        exit(1);

    c = allocatem(rows, collumns);
    Bi = (B3 < rows)? B3 : rows;
    Bj = (B3 < collumns)? B3 : collumns;

    for(i = 0; i < rows; i = i+Bi) {
        Bi = (Bi < rows-i)? Bi : rows-i;
        for(j = 0; j < collumns; j = j+Bj) {
            Bj = (Bj < collumns-j)? Bj : collumns-j;
            for(m = 0; m < Bi; m++) {
                for(n = 0; n < Bj; n++) {
                    c->matrix[i+m][j+n] = m_t1->matrix[i+m][j+n] - m_t2->matrix[i+m][j+n];
                }
            }
        }
    }
    return c;
}

matrix_t *sub2(matrix_t *m_t1, matrix_t *m_t2) {
    matrix_t *c = NULL;
    int i, j;
    int rows = m_t1->rows;
    int collumns = m_t1->collumns;

    if((m_t1->rows != m_t2->rows)&&(m_t1->collumns != m_t2->collumns))
        exit(1);

    c = allocatem(rows, collumns);

    for(i = 0; i < rows; i = i+1) {
        for(j = 0; j < collumns; j = j+1) {
            c->matrix[i][j] = m_t1->matrix[i][j] - m_t2->matrix[i][j];
        }
    }
    return c;
}

matrix_t *subs(matrix_t *m_t1, double scalar) {
    matrix_t *c = NULL;
    int i, j;
    int rows = m_t1->rows;
    int collumns = m_t1->collumns;

    c = allocatem(rows, collumns);

    for(i = 0; i < rows; i = i+1) {
        for(j = 0; j < collumns; j = j+1) {
            c->matrix[i][j] = m_t1->matrix[i][j] - scalar;
        }
    }
    return c;    
}

double det(matrix_t *m_t1) {
    double res = 0;
    int i, m, n;
    if(m_t1->rows == 2 && m_t1->collumns == 2) {
        res = (m_t1->matrix[0][0]*m_t1->matrix[1][1]) - (m_t1->matrix[1][0]*m_t1->matrix[0][1]);
        return res;
    }

    for(i = 0; i < m_t1->collumns; i = i+1) {
        matrix_t *c = allocatem(m_t1->rows-1, m_t1->collumns-1);
        for(m = 1; m < m_t1->rows; m = m+1) {
            int k = 0;
            for(n = 0; n < m_t1->collumns; n = n+1) {
                if(n == i)
                    continue;
                c->matrix[m-1][k] = m_t1->matrix[m][n];
                k++;
            }
        }
        res = res + (pow(-1, i) * m_t1->matrix[0][i] * det(c));
        releasem(c);
    }
    return res;
}
//
void eliminate(double *r1, double *r2, int col, int r_len, int target, double *ri1, double *ri2) {
    double fact = (double)(r2[col]-target) / r1[col];
    for(int i = 0; i < r_len; i = i+1) {
        r2[i] -= (double)fact * r1[i];
        ri2[i] -= (double)fact * ri1[i];
    }
}
//
matrix_t *inverse(matrix_t *m_t1) {
    int i, j, n = m_t1->rows;
    double max_e;
    matrix_t *m_alpha = copy(m_t1); 
    matrix_t *m_identity = identity(m_t1->rows, m_t1->collumns);

    if(det(m_alpha) == 0)
        exit(1);

    for(i = 0; i < m_alpha->rows; i = i+1) {
        for(j = i+1; j < m_alpha->collumns; j = j+1)
            eliminate(m_alpha->matrix[i], m_alpha->matrix[j], i, m_alpha->collumns, 0, m_identity->matrix[i], m_identity->matrix[j]);
    }

    for(i = m_alpha->rows-1; i >= 0; i = i-1) {
        for(j = i-1; j >= 0; j = j-1)
            eliminate(m_alpha->matrix[i], m_alpha->matrix[j], i, m_alpha->collumns, 0, m_identity->matrix[i], m_identity->matrix[j]);
    }

    for(i = 0; i < m_alpha->collumns; i = i+1) {
        eliminate(m_alpha->matrix[i], m_alpha->matrix[i], i, m_alpha->collumns, 1, m_identity->matrix[i], m_identity->matrix[i]);
    }

    return m_identity;
}