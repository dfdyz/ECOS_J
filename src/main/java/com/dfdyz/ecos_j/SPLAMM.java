package com.dfdyz.ecos_j;

public class SPLAMM {

    public static void splaCumsum(int[] p, int[] w, int m) {
        int cumsum = 0;
        for (int i = 0; i < m; i++) {
            p[i] = cumsum;
            cumsum += w[i];
            w[i] = p[i];
        }
    }

    public static void pinv(int n, int[] p, int[] pinv) {
        for (int i = 0; i < n; i++) { pinv[p[i]] = i; }
    }

    public static SpMat transposeSparseMatrix(SpMat M, int[] MtoMt) {
        int j, i, k, q;
        int[] w;

        SpMat A = newSparseMatrix(M.n, M.m, M.nnz);
        if (M.nnz == 0) return A;

        w = new int[M.m];
        for (i = 0; i < M.m; i++) { w[i] = 0; }
        for (k = 0; k < M.nnz; k++) { w[M.ir[k]]++; }

        splaCumsum(A.jc, w, M.m);

        for (j = 0; j < M.n; j++) {
            for (k = M.jc[j]; k < M.jc[j + 1]; k++) {
                q = w[M.ir[k]]++;
                A.ir[q] = j;
                A.pr[q] = M.pr[k];
                MtoMt[k] = q;
            }
        }

        return A;
    }

    public static SpMat newSparseMatrix(int m, int n, int nnz) {
        int[] jc = new int[n + 1];
        int[] ir = new int[nnz];
        double[] pr = new double[nnz];
        jc[n] = nnz;
        return ecosCreateSparseMatrix(m, n, nnz, jc, ir, pr);
    }

    public static SpMat ecosCreateSparseMatrix(int m, int n, int nnz, int[] jc, int[] ir, double[] pr) {
        SpMat M = new SpMat();
        M.m = m;
        M.n = n;
        M.nnz = nnz;
        M.jc = jc;
        M.ir = ir;
        M.pr = pr;
        if (M.jc != null) M.jc[n] = nnz;
        return M;
    }

    public static void freeSparseMatrix(SpMat M) {
        M.ir = null;
        M.jc = null;
        M.pr = null;
    }

    public static void permuteSparseSymmetricMatrix(SpMat A, int[] pinv, SpMat C, int[] PK) {
        int i, i2, j, j2, k, q;
        int[] w;

        w = new int[A.n];
        for (j = 0; j < A.n; j++) { w[j] = 0; }
        for (j = 0; j < A.n; j++) {
            j2 = pinv[j];
            for (k = A.jc[j]; k < A.jc[j + 1]; k++) {
                i = A.ir[k];
                if (i > j) continue;
                i2 = pinv[i];
                w[(i2 > j2 ? i2 : j2)]++;
            }
        }

        splaCumsum(C.jc, w, A.n);

        for (j = 0; j < A.n; j++) {
            j2 = pinv[j];
            for (k = A.jc[j]; k < A.jc[j + 1]; k++) {
                i = A.ir[k];
                if (i > j) continue;
                i2 = pinv[i];
                q = w[(i2 > j2 ? i2 : j2)]++;
                C.ir[q] = (i2 < j2 ? i2 : j2);
                C.pr[q] = A.pr[k];
                if (PK != null) PK[k] = q;
            }
        }
    }

    public static SpMat copySparseMatrix(SpMat A) {
        int i;
        SpMat B = newSparseMatrix(A.m, A.n, A.nnz);
        for (i = 0; i <= A.n; i++) { B.jc[i] = A.jc[i]; }
        for (i = 0; i < A.nnz; i++) { B.ir[i] = A.ir[i]; }
        for (i = 0; i < A.nnz; i++) { B.pr[i] = A.pr[i]; }
        return B;
    }

    public static void printDenseMatrix(double[] M, int dim1, int dim2, String name) {
        if (GlblOpts.PRINTLEVEL > 2) {
            GlblOpts.printText("%s = \n\t", name);
            for (int i = 0; i < dim1; i++) {
                for (int j = 0; j < dim2; j++) {
                    if (j < dim2 - 1)
                        GlblOpts.printText("% 14.12e,  ", M[i * dim2 + j]);
                    else
                        GlblOpts.printText("% 14.12e;  ", M[i * dim2 + j]);
                }
                if (i < dim1 - 1) {
                    GlblOpts.printText("\n\t");
                }
            }
            GlblOpts.printText("\n");
        }
    }

    public static void printDenseMatrixI(int[] M, int dim1, int dim2, String name) {
        if (GlblOpts.PRINTLEVEL > 2) {
            GlblOpts.printText("%s = \n\t", name);
            for (int i = 0; i < dim1; i++) {
                for (int j = 0; j < dim2; j++) {
                    if (j < dim2 - 1)
                        GlblOpts.printText("%d,  ", M[i * dim2 + j]);
                    else
                        GlblOpts.printText("%d;  ", M[i * dim2 + j]);
                }
                if (i < dim1 - 1) {
                    GlblOpts.printText("\n\t");
                }
            }
            GlblOpts.printText("\n");
        }
    }

    public static void printSparseMatrix(SpMat M) {
        if (GlblOpts.PRINTLEVEL > 2) {
            int k = 0;
            for (int j = 0; j < M.n; j++) {
                int row_strt = M.jc[j];
                int row_stop = M.jc[j + 1];
                if (row_strt == row_stop) continue;
                else {
                    for (int i = row_strt; i < row_stop; i++) {
                        GlblOpts.printText("\t(%3d,%3d) = %g\n", M.ir[i] + 1, j + 1, M.pr[k++]);
                    }
                }
            }
        }
    }

    public static void dumpSparseMatrix(SpMat M, String fn) {
        if (GlblOpts.PRINTLEVEL > 2) {
            int k = 0;
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < M.n; j++) {
                int row_strt = M.jc[j];
                int row_stop = M.jc[j + 1];
                if (row_strt == row_stop) continue;
                else {
                    for (int i = row_strt; i < row_stop; i++) {
                        sb.append(String.format("%d\t%d\t%20.18e\n", M.ir[i] + 1, j + 1, M.pr[k++]));
                    }
                }
            }
            sb.append(String.format("%d\t%d\t%20.18e\n", M.m, M.n, 0.0));
            System.out.print(sb.toString());
        }
    }
}
