extern double A[N][M];
extern double B[M][N];
extern double C[N][N];

/**
 * The basic dense matrix multiply.
 * Use this as a stub for your other implementations.
 */
void basic_MM()
{
    int i, j, k;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            C[i][j] = 0;
            for (k = 0; k < M; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
