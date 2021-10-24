#include "num_c.h"

int main() {
    matrix_t *m_t1 = identity(4, 4);
    printm(m_t1);
    matrix_t *m_t2 = ones(4, 4);
    printf("\n");

    matrix_t *m_t3 = sub(m_t1, m_t2);
    matrix_t *m_t5 = adds(m_t3, 7.0f);
    printm(m_t3);
    printf("\n");
    double temp = 0;
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            m_t3->matrix[i][j] = ++temp;
            if(i == j)
                m_t3->matrix[i][j] = 1;
        }
    }
    printm(m_t3);
    printf("\n");
    double b = det(m_t3);
    printf("DET(m_t3) = %f\n", b);

    matrix_t *alpha = inverse(m_t3);
    printm(alpha);
    printf("\n");
    matrix_t *m_t4 = multiply(alpha, m_t3);
    printm(m_t4);
    releasem(m_t1);
    releasem(m_t2);
    releasem(m_t3);
    releasem(m_t5);
    releasem(m_t4);
    releasem(alpha);
}