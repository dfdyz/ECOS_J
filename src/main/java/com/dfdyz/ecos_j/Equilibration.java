package com.dfdyz.ecos_j;

public class Equilibration {

    static void max_rows(double[] E, SpMat mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = mat.jc[i]; j < mat.jc[i + 1]; j++) {
                int row = mat.ir[j];
                E[row] = Math.max(GlblOpts.fabs(mat.pr[j]), E[row]);
            }
        }
    }

    static void max_cols(double[] E, SpMat mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = mat.jc[i]; j < mat.jc[i + 1]; j++) {
                E[i] = Math.max(GlblOpts.fabs(mat.pr[j]), E[i]);
            }
        }
    }

    static void sum_sq_rows(double[] E, SpMat mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = mat.jc[i]; j < mat.jc[i + 1]; j++) {
                int row = mat.ir[j];
                E[row] += (mat.pr[j] * mat.pr[j]);
            }
        }
    }

    static void sum_sq_cols(double[] E, SpMat mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = mat.jc[i]; j < mat.jc[i + 1]; j++) {
                E[i] += (mat.pr[j] * mat.pr[j]);
            }
        }
    }

    static void equilibrate_rows(double[] E, SpMat mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = mat.jc[i]; j < mat.jc[i + 1]; j++) {
                int row = mat.ir[j];
                mat.pr[j] /= E[row];
            }
        }
    }

    static void equilibrate_cols(double[] E, SpMat mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = mat.jc[i]; j < mat.jc[i + 1]; j++) {
                mat.pr[j] /= E[i];
            }
        }
    }

    static void restore(double[] D, double[] E, SpMat mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = mat.jc[i]; j < mat.jc[i + 1]; j++) {
                int row = mat.ir[j];
                mat.pr[j] *= (D[row] * E[i]);
            }
        }
    }

    static void use_ruiz_equilibration(PWork w) {
        int i, j, ind, iter;
        int num_cols = (w.A != null) ? w.A.n : w.G.n;
        int num_A_rows = (w.A != null) ? w.A.m : 0;
        int num_G_rows = w.G.m;
        double[] xtmp = new double[num_cols];
        double[] Atmp = new double[num_A_rows];
        double[] Gtmp = new double[num_G_rows];
        double total;

        for (i = 0; i < num_cols; i++) { w.xequil[i] = 1.0; }
        for (i = 0; i < num_A_rows; i++) { w.Aequil[i] = 1.0; }
        for (i = 0; i < num_G_rows; i++) { w.Gequil[i] = 1.0; }

        for (iter = 0; iter < GlblOpts.EQUIL_ITERS; iter++) {
            for (i = 0; i < num_cols; i++) { xtmp[i] = 0.0; }
            for (i = 0; i < num_A_rows; i++) { Atmp[i] = 0.0; }
            for (i = 0; i < num_G_rows; i++) { Gtmp[i] = 0.0; }

            if (w.A != null) max_cols(xtmp, w.A);
            if (num_G_rows > 0) max_cols(xtmp, w.G);

            if (w.A != null) max_rows(Atmp, w.A);
            if (num_G_rows > 0) max_rows(Gtmp, w.G);

            ind = w.C.lpc.p;
            for (i = 0; i < w.C.nsoc; i++) {
                total = 0.0;
                for (j = 0; j < w.C.soc[i].p; j++) { total += Gtmp[ind + j]; }
                for (j = 0; j < w.C.soc[i].p; j++) { Gtmp[ind + j] = total; }
                ind += w.C.soc[i].p;
            }

            for (i = 0; i < w.C.nexc; i++) {
                total = 0.0;
                for (j = 0; j < 3; j++) { total += Gtmp[ind + j]; }
                for (j = 0; j < 3; j++) { Gtmp[ind + j] = total; }
                ind += 3;
            }

            for (i = 0; i < num_cols; i++) {
                xtmp[i] = GlblOpts.fabs(xtmp[i]) < 1e-6 ? 1.0 : Math.sqrt(xtmp[i]);
            }
            for (i = 0; i < num_A_rows; i++) {
                Atmp[i] = GlblOpts.fabs(Atmp[i]) < 1e-6 ? 1.0 : Math.sqrt(Atmp[i]);
            }
            for (i = 0; i < num_G_rows; i++) {
                Gtmp[i] = GlblOpts.fabs(Gtmp[i]) < 1e-6 ? 1.0 : Math.sqrt(Gtmp[i]);
            }

            if (w.A != null) equilibrate_rows(Atmp, w.A);
            if (num_G_rows > 0) equilibrate_rows(Gtmp, w.G);

            if (w.A != null) equilibrate_cols(xtmp, w.A);
            if (num_G_rows > 0) equilibrate_cols(xtmp, w.G);

            for (i = 0; i < num_cols; i++) { w.xequil[i] *= xtmp[i]; }
            for (i = 0; i < num_A_rows; i++) { w.Aequil[i] *= Atmp[i]; }
            for (i = 0; i < num_G_rows; i++) { w.Gequil[i] *= Gtmp[i]; }
        }

        for (i = 0; i < num_A_rows; i++) { w.b[i] /= w.Aequil[i]; }
        for (i = 0; i < num_G_rows; i++) { w.h[i] /= w.Gequil[i]; }
    }

    public static void set_equilibration(PWork w) {
        use_ruiz_equilibration(w);
    }

    public static void unset_equilibration(PWork w) {
        int i;
        int num_A_rows = (w.A != null) ? w.A.m : 0;
        int num_G_rows = w.G.m;

        if (w.A != null) restore(w.Aequil, w.xequil, w.A);
        if (num_G_rows > 0) restore(w.Gequil, w.xequil, w.G);

        for (i = 0; i < num_A_rows; i++) { w.b[i] *= w.Aequil[i]; }
        for (i = 0; i < num_G_rows; i++) { w.h[i] *= w.Gequil[i]; }
    }
}
