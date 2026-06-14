package com.dfdyz.ecos_j;

public class ExpConeOps {

    public static double socres(double[] u, int p) {
        double res = u[0] * u[0];
        for (int i = 1; i < p; i++) { res -= u[i] * u[i]; }
        return res;
    }

    public static void evalExpHessian(double[] w, double[] v, double mu) {
        double x = w[0];
        double y = w[1];
        double z = w[2];
        double l = Math.log(-y / x);
        double r = -x * l - x + z;
        v[0] = mu * ((r * r - x * r + l * l * x * x) / (r * x * x * r));
        v[1] = mu * ((z - x) / (r * r * y));
        v[2] = mu * ((r * r - x * r + x * x) / (r * r * y * y));
        v[3] = mu * (-l / (r * r));
        v[4] = mu * (-x / (r * r * y));
        v[5] = mu * (1 / (r * r));
    }

    public static void evalExpGradient(double[] w, double[] g) {
        double x = w[0];
        double y = w[1];
        double z = w[2];
        double l = Math.log(-y / x);
        double r = -x * l - x + z;
        g[0] = (l * x - r) / (r * x);
        g[1] = (x - r) / (y * r);
        g[2] = -1 / r;
    }

    public static double evalBarrierValue(double[] siter, double[] ziter, int fc, int nexc) {
        double l, u, v, w, x, y, z, o;
        double primal_barrier = 0.0;
        double dual_barrier = 0.0;

        int j;
        double[] ziterOffset = new double[ziter.length];
        double[] siterOffset = new double[siter.length];
        System.arraycopy(ziter, 0, ziterOffset, 0, ziter.length);
        System.arraycopy(siter, 0, siterOffset, 0, siter.length);
        int zIdx = fc;
        int sIdx = fc;

        for (j = 0; j < nexc; j++) {
            u = ziterOffset[zIdx];
            v = ziterOffset[zIdx + 1];
            w = ziterOffset[zIdx + 2];
            x = siterOffset[sIdx];
            y = siterOffset[sIdx + 1];
            z = siterOffset[sIdx + 2];

            l = Math.log(-v / u);
            dual_barrier += -Math.log(w - u - u * l) - Math.log(-u) - Math.log(v);

            o = WrightOmega.wrightOmega(1 - x / z - Math.log(z) + Math.log(y));
            o = (o - 1) * (o - 1) / o;
            primal_barrier += -Math.log(o) - 2 * Math.log(z) - Math.log(y) - 3;

            zIdx += 3;
            sIdx += 3;
        }
        return primal_barrier + dual_barrier;
    }

    public static void scaleToAddExpcone(double[] y, double[] x, ExpCone[] expc, int nexc, int fc) {
        int l;
        int xIdx = fc;
        int yIdx = fc;

        for (l = 0; l < nexc; l++) {
            y[yIdx] += expc[l].v[0] * x[xIdx] + expc[l].v[1] * x[xIdx + 1] + expc[l].v[3] * x[xIdx + 2];
            y[yIdx + 1] += expc[l].v[1] * x[xIdx] + expc[l].v[2] * x[xIdx + 1] + expc[l].v[4] * x[xIdx + 2];
            y[yIdx + 2] += expc[l].v[3] * x[xIdx] + expc[l].v[4] * x[xIdx + 1] + expc[l].v[5] * x[xIdx + 2];

            xIdx += 3;
            yIdx += 3;
        }
    }

    public static int evalExpPrimalFeas(double[] s, int fc, int nexc) {
        double x1, x2, x3, tmp1, psi;
        for (int j = 0; j < nexc; j++) {
            int idx = fc + 3 * j;
            x1 = s[idx];
            x2 = s[idx + 1];
            x3 = s[idx + 2];
            tmp1 = Math.log(x2 / x3);
            psi = x3 * tmp1 - x1;
            if (psi < 0 || x2 < 0 || x3 < 0) {
                return 0;
            }
        }
        return 1;
    }

    public static int evalExpDualFeas(double[] z, int fc, int nexc) {
        double x1, x2, x3, tmp1, psi;
        for (int j = 0; j < nexc; j++) {
            int idx = fc + 3 * j;
            x1 = z[idx];
            x2 = z[idx + 1];
            x3 = z[idx + 2];
            tmp1 = Math.log(-x2 / x1);
            psi = -x1 - x1 * tmp1 + x3;
            if (0 < x1 || x2 < 0 || psi < 0) {
                return 0;
            }
        }
        return 1;
    }
}
