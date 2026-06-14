package com.dfdyz.ecos_j;

public class LDL {

    public static void ldl_symbolic2(int n, int[] Ap, int[] Ai, int[] Lp, int[] Parent, int[] Lnz, int[] Flag) {
        int i, k, p, p2;

        for (k = 0; k < n; k++) {
            Parent[k] = -1;
            Flag[k] = k;
            Lnz[k] = 0;
            p2 = Ap[k + 1];
            for (p = Ap[k]; p < p2; p++) {
                i = Ai[p];
                for (; Flag[i] != k; i = Parent[i]) {
                    if (Parent[i] == -1) Parent[i] = k;
                    Lnz[i]++;
                    Flag[i] = k;
                }
            }
        }

        Lp[0] = 0;
        for (k = 0; k < n; k++) {
            Lp[k + 1] = Lp[k] + Lnz[k];
        }
    }

    public static int ldl_numeric2(int n, int[] Ap, int[] Ai, double[] Ax,
                                   int[] Lp, int[] Parent, int[] Sign, double eps, double delta,
                                   int[] Lnz, int[] Li, double[] Lx,
                                   double[] D, double[] Y, int[] Pattern, int[] Flag,
                                   double[] t1, double[] t2) {
        double yi, l_ki;
        int i, k, p, p2, len, top;
        Timer clock = new Timer();

        for (k = 0; k < n; k++) {
            if (GlblOpts.PROFILING > 1) TimerOps.tic(clock);

            Y[k] = 0.0;
            top = n;
            Flag[k] = k;
            Lnz[k] = 0;
            p2 = Ap[k + 1];
            for (p = Ap[k]; p < p2; p++) {
                i = Ai[p];
                Y[i] = Ax[p];
                for (len = 0; Flag[i] != k; i = Parent[i]) {
                    Pattern[len++] = i;
                    Flag[i] = k;
                }
                while (len > 0) Pattern[--top] = Pattern[--len];
            }

            if (GlblOpts.PROFILING > 1) {
                t1[0] += TimerOps.toc(clock);
                TimerOps.tic(clock);
            }

            D[k] = Y[k];
            Y[k] = 0.0;
            for (; top < n; top++) {
                i = Pattern[top];
                yi = Y[i];
                Y[i] = 0.0;
                p2 = Lp[i] + Lnz[i];
                for (p = Lp[i]; p < p2; p++) {
                    Y[Li[p]] -= Lx[p] * yi;
                }
                l_ki = yi / D[i];
                D[k] -= l_ki * yi;
                Li[p] = k;
                Lx[p] = l_ki;
                Lnz[i]++;
            }

            D[k] = (Sign[k] * D[k] <= eps ? Sign[k] * delta : D[k]);

            if (GlblOpts.PROFILING > 1) t2[0] += TimerOps.toc(clock);
        }
        return n;
    }

    public static int ldl_numeric2(int n, int[] Ap, int[] Ai, double[] Ax,
                                   int[] Lp, int[] Parent, int[] Sign, double eps, double delta,
                                   int[] Lnz, int[] Li, double[] Lx,
                                   double[] D, double[] Y, int[] Pattern, int[] Flag) {
        double yi, l_ki;
        int i, k, p, p2, len, top;

        for (k = 0; k < n; k++) {
            Y[k] = 0.0;
            top = n;
            Flag[k] = k;
            Lnz[k] = 0;
            p2 = Ap[k + 1];
            for (p = Ap[k]; p < p2; p++) {
                i = Ai[p];
                Y[i] = Ax[p];
                for (len = 0; Flag[i] != k; i = Parent[i]) {
                    Pattern[len++] = i;
                    Flag[i] = k;
                }
                while (len > 0) Pattern[--top] = Pattern[--len];
            }

            D[k] = Y[k];
            Y[k] = 0.0;
            for (; top < n; top++) {
                i = Pattern[top];
                yi = Y[i];
                Y[i] = 0.0;
                p2 = Lp[i] + Lnz[i];
                for (p = Lp[i]; p < p2; p++) {
                    Y[Li[p]] -= Lx[p] * yi;
                }
                l_ki = yi / D[i];
                D[k] -= l_ki * yi;
                Li[p] = k;
                Lx[p] = l_ki;
                Lnz[i]++;
            }

            D[k] = (Sign[k] * D[k] <= eps ? Sign[k] * delta : D[k]);
        }
        return n;
    }

    public static void ldl_lsolve(int n, double[] X, int[] Lp, int[] Li, double[] Lx) {
        int j, p, p2;
        for (j = 0; j < n; j++) {
            p2 = Lp[j + 1];
            for (p = Lp[j]; p < p2; p++) {
                X[Li[p]] -= Lx[p] * X[j];
            }
        }
    }

    public static void ldl_lsolve2(int n, double[] b, int[] Lp, int[] Li, double[] Lx, double[] x) {
        int j, p, p2;
        for (j = 0; j < n; j++) { x[j] = b[j]; }
        for (j = 0; j < n; j++) {
            p2 = Lp[j + 1];
            for (p = Lp[j]; p < p2; p++) {
                x[Li[p]] -= Lx[p] * x[j];
            }
        }
    }

    public static void ldl_dsolve(int n, double[] X, double[] D) {
        for (int j = 0; j < n; j++) { X[j] /= D[j]; }
    }

    public static void ldl_ltsolve(int n, double[] X, int[] Lp, int[] Li, double[] Lx) {
        int j, p, p2;
        for (j = n - 1; j >= 0; j--) {
            p2 = Lp[j + 1];
            for (p = Lp[j]; p < p2; p++) {
                X[j] -= Lx[p] * X[Li[p]];
            }
        }
    }
}
