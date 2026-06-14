package com.dfdyz.ecos_j;

public class SPLA {

    public static void sparseMV(SpMat A, double[] x, double[] y, int a, int newVector) {
        int i, j;
        if (newVector > 0) {
            for (i = 0; i < A.m; i++) { y[i] = 0; }
        }
        if (A.nnz == 0) return;
        if (a > 0) {
            for (j = 0; j < A.n; j++) {
                for (i = A.jc[j]; i < A.jc[j + 1]; i++) {
                    y[A.ir[i]] += A.pr[i] * x[j];
                }
            }
        } else {
            for (j = 0; j < A.n; j++) {
                for (i = A.jc[j]; i < A.jc[j + 1]; i++) {
                    y[A.ir[i]] -= A.pr[i] * x[j];
                }
            }
        }
    }

    public static void sparseMtVm(SpMat A, double[] x, double[] y, int newVector, int skipDiagonal) {
        int i, j, k;
        if (newVector > 0) {
            for (j = 0; j < A.n; j++) { y[j] = 0; }
        }
        if (A.nnz == 0) return;
        if (skipDiagonal != 0) {
            for (j = 0; j < A.n; j++) {
                for (k = A.jc[j]; k < A.jc[j + 1]; k++) {
                    i = A.ir[k];
                    y[j] -= (i == j ? 0 : A.pr[k] * x[i]);
                }
            }
        } else {
            for (j = 0; j < A.n; j++) {
                for (k = A.jc[j]; k < A.jc[j + 1]; k++) {
                    y[j] -= A.pr[k] * x[A.ir[k]];
                }
            }
        }
    }

    public static void vadd(int n, double[] x, double[] y) {
        for (int i = 0; i < n; i++) { y[i] += x[i]; }
    }

    public static void vsubscale(int n, double a, double[] x, double[] y) {
        for (int i = 0; i < n; i++) { y[i] -= a * x[i]; }
    }

    public static double norm2(double[] v, int n) {
        double normsquare = 0;
        for (int i = 0; i < n; i++) { normsquare += v[i] * v[i]; }
        return Math.sqrt(normsquare);
    }

    public static double norminf(double[] v, int n) {
        double norm = 0;
        double mv;
        for (int i = 0; i < n; i++) {
            if (v[i] > norm) { norm = v[i]; }
            mv = -v[i];
            if (mv > norm) { norm = mv; }
        }
        return norm;
    }

    public static double eddot(int n, double[] x, double[] y) {
        double z = 0;
        for (int i = 0; i < n; i++) { z += x[i] * y[i]; }
        return z;
    }

    public static double eddot(int n, double[] x, double[] y, int xOffset, int yOffset) {
        double z = 0;
        for (int i = 0; i < n; i++) { z += x[xOffset + i] * y[yOffset + i]; }
        return z;
    }

    public static double norm2(double[] v, int offset, int n) {
        double normsquare = 0;
        for (int i = 0; i < n; i++) { normsquare += v[offset + i] * v[offset + i]; }
        return Math.sqrt(normsquare);
    }
}
