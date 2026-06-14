package com.dfdyz.ecos_j;

public class ConeOps {

    public static void bring2cone(Cone C, double[] r, double[] s) {
        double alpha = -GlblOpts.GAMMA;
        double cres, r1square;
        int i, l, j;

        i = 0;
        for (int ii = 0; ii < C.lpc.p; ii++) {
            if (r[i] <= 0 && -r[i] > alpha) { alpha = -r[i]; }
            i++;
        }

        for (l = 0; l < C.nsoc; l++) {
            cres = r[i]; i++;
            r1square = 0;
            for (j = 1; j < C.soc[l].p; j++) { r1square += r[i] * r[i]; i++; }
            cres -= Math.sqrt(r1square);
            if (cres <= 0 && -cres > alpha) { alpha = -cres; }
        }

        alpha += 1.0;
        i = 0;

        for (int ii = 0; ii < C.lpc.p; ii++) {
            s[i] = r[i] + alpha;
            i++;
        }

        for (l = 0; l < C.nsoc; l++) {
            s[i] = r[i] + alpha; i++;
            for (j = 1; j < C.soc[l].p; j++) { s[i] = r[i]; i++; }
        }
    }

    public static void unitInitialization(Cone C, double[] s, double[] z, double scaling) {
        int i = 0, l, j;

        for (int ii = 0; ii < C.lpc.p; ii++) {
            s[i] = scaling;
            z[i] = scaling;
            i++;
        }

        for (l = 0; l < C.nsoc; l++) {
            s[i] = scaling; z[i] = scaling; i++;
            for (j = 1; j < C.soc[l].p; j++) { s[i] = 0.0; z[i] = 0.0; i++; }
        }

        for (l = 0; l < C.nexc; l++) {
            s[i] = scaling * (-1.051383945322714);
            s[i + 1] = scaling * (1.258967884768947);
            s[i + 2] = scaling * (0.556409619469370);
            z[i] = scaling * (-1.051383945322714);
            z[i + 1] = scaling * (1.258967884768947);
            z[i + 2] = scaling * (0.556409619469370);
            i = i + 3;
        }
    }

    public static int updateScalings(Cone C, double[] s, double[] z, double[] lambda, double mu) {
        int i, l, k, p;
        double sres, zres, snorm, znorm, gamma, one_over_2gamma;
        double[] sk;
        double[] zk;
        double a, c, d, w, temp, divisor;
        double u0, u0_square, u1, v1, d1, c2byu02_d, c2byu02, c_square;

        i = 0;
        for (int ii = 0; ii < C.lpc.p; ii++) {
            C.lpc.v[i] = GlblOpts.safeDivPos(s[i], z[i]);
            C.lpc.w[i] = Math.sqrt(C.lpc.v[i]);
            i++;
        }

        k = C.lpc.p;
        for (l = 0; l < C.nsoc; l++) {
            sk = s;
            zk = z;
            int skOffset = k;
            int zkOffset = k;
            p = C.soc[l].p;

            sres = socres(s, p, skOffset);
            zres = socres(z, p, zkOffset);
            if (sres <= 0 || zres <= 0) { return GlblOpts.OUTSIDE_CONE; }

            snorm = Math.sqrt(sres);
            znorm = Math.sqrt(zres);
            for (i = 0; i < p; i++) { C.soc[l].skbar[i] = GlblOpts.safeDivPos(s[skOffset + i], snorm); }
            for (i = 0; i < p; i++) { C.soc[l].zkbar[i] = GlblOpts.safeDivPos(z[zkOffset + i], znorm); }
            C.soc[l].eta_square = GlblOpts.safeDivPos(snorm, znorm);
            C.soc[l].eta = Math.sqrt(C.soc[l].eta_square);

            gamma = 1.0;
            for (i = 0; i < p; i++) { gamma += C.soc[l].skbar[i] * C.soc[l].zkbar[i]; }
            gamma = Math.sqrt(0.5 * gamma);
            one_over_2gamma = GlblOpts.safeDivPos(0.5, gamma);
            a = one_over_2gamma * (C.soc[l].skbar[0] + C.soc[l].zkbar[0]);
            w = 0;
            for (i = 1; i < p; i++) {
                C.soc[l].q[i - 1] = one_over_2gamma * (C.soc[l].skbar[i] - C.soc[l].zkbar[i]);
                w += C.soc[l].q[i - 1] * C.soc[l].q[i - 1];
            }
            C.soc[l].w = w;
            C.soc[l].a = a;

            temp = 1.0 + a;
            c = 1.0 + a + GlblOpts.safeDivPos(w, temp);
            divisor = temp * temp;
            d = 1 + GlblOpts.safeDivPos(2, temp) + GlblOpts.safeDivPos(w, divisor);

            c_square = c * c;
            divisor = 1.0 + w * d;
            d1 = 0.5 * (a * a + w * (1.0 - GlblOpts.safeDivPos(c_square, divisor)));
            if (d1 < 0) { d1 = 0; }
            u0_square = a * a + w - d1;
            u0 = Math.sqrt(u0_square);
            c2byu02 = GlblOpts.safeDivPos((c * c), u0_square);
            c2byu02_d = c2byu02 - d;
            if (c2byu02_d <= 0) { return GlblOpts.OUTSIDE_CONE; }
            v1 = Math.sqrt(c2byu02_d);
            u1 = Math.sqrt(c2byu02);
            C.soc[l].d1 = d1;
            C.soc[l].u0 = u0;
            C.soc[l].u1 = u1;
            C.soc[l].v1 = v1;

            k += C.soc[l].p;
        }

        k = C.fexv;
        for (l = 0; l < C.nexc; l++) {
            evalExpHessian(z, k, C.expc[l].v, mu);
            evalExpGradient(z, k, C.expc[l].g);
            k += 3;
        }

        scale(z, C, lambda);

        return GlblOpts.INSIDE_CONE;
    }

    public static double evalSymmetricBarrierValue(double[] siter, double[] ziter, double tauIter, double kapIter, Cone C, double D) {
        double barrier = 0.0;
        int j, k = 0, l;
        double normAccumS = 0.0;
        double normAccumZ = 0.0;
        int socDim;

        for (k = 0; k < C.lpc.p; k++)
            barrier -= (siter[k] <= 0 || ziter[k] <= 0) ? GlblOpts.ECOS_INFINITY : (Math.log(siter[k]) + Math.log(ziter[k]));

        barrier -= (tauIter <= 0 || kapIter <= 0) ? GlblOpts.ECOS_INFINITY : (Math.log(tauIter) + Math.log(kapIter));

        for (l = 0; l < C.nsoc; l++) {
            socDim = C.soc[l].p;
            normAccumS = 0.0;
            normAccumZ = 0.0;
            normAccumS = siter[k] * siter[k];
            normAccumZ = ziter[k] * ziter[k];
            k++;
            for (j = 1; j < socDim; j++) {
                normAccumS -= siter[k] * siter[k];
                normAccumZ -= ziter[k] * ziter[k];
                k++;
            }
            barrier -= normAccumS <= 0.0 ? GlblOpts.ECOS_INFINITY : 0.5 * Math.log(normAccumS);
            barrier -= normAccumZ <= 0.0 ? GlblOpts.ECOS_INFINITY : 0.5 * Math.log(normAccumZ);
        }
        return barrier - D - 1;
    }

    public static void scale(double[] z, Cone C, double[] lambda) {
        int i, j, l, cone_start;
        double zeta, factor;

        for (i = 0; i < C.lpc.p; i++) { lambda[i] = C.lpc.w[i] * z[i]; }

        cone_start = C.lpc.p;
        for (l = 0; l < C.nsoc; l++) {
            zeta = 0;
            for (i = 1; i < C.soc[l].p; i++) { zeta += C.soc[l].q[i - 1] * z[cone_start + i]; }

            factor = z[cone_start] + GlblOpts.safeDivPos(zeta, (1 + C.soc[l].a));

            lambda[cone_start] = C.soc[l].eta * (C.soc[l].a * z[cone_start] + zeta);
            for (i = 1; i < C.soc[l].p; i++) {
                j = cone_start + i;
                lambda[j] = C.soc[l].eta * (z[j] + factor * C.soc[l].q[i - 1]);
            }

            cone_start += C.soc[l].p;
        }
    }

    public static void scale2add(double[] x, double[] y, Cone C) {
        int i, l, cone_start, conesize, conesize_m1;
        double[] x1arr, y1arr;
        double eta_square;
        double[] q;
        double d1, u0, u1, v1;
        double v1x3_plus_u1x4;
        double qtx2;

        for (i = 0; i < C.lpc.p; i++) { y[i] += C.lpc.v[i] * x[i]; }

        cone_start = C.lpc.p;

        for (l = 0; l < C.nsoc; l++) {
            SOCone soc = C.soc[l];
            conesize = soc.p;
            conesize_m1 = conesize - 1;
            eta_square = soc.eta_square;
            d1 = soc.d1;
            u0 = soc.u0;
            u1 = soc.u1;
            v1 = soc.v1;
            q = soc.q;

            int xOff = cone_start;
            int yOff = cone_start;

            y[yOff] += eta_square * (d1 * x[xOff] + u0 * x[xOff + conesize + 1]);

            v1x3_plus_u1x4 = v1 * x[xOff + conesize] + u1 * x[xOff + conesize + 1];
            qtx2 = 0;
            for (i = 0; i < conesize_m1; i++) {
                y[yOff + 1 + i] += eta_square * (x[xOff + 1 + i] + v1x3_plus_u1x4 * q[i]);
                qtx2 += q[i] * x[xOff + 1 + i];
            }

            y[yOff + conesize] += eta_square * (v1 * qtx2 + x[xOff + conesize]);
            y[yOff + conesize + 1] += eta_square * (u0 * x[xOff] + u1 * qtx2 - x[xOff + conesize + 1]);

            cone_start += conesize + 2;
        }

        ExpConeOps.scaleToAddExpcone(y, x, C.expc, C.nexc, cone_start);
    }

    public static void unscale(double[] lambda, Cone C, double[] z) {
        int i, j, l, cone_start;
        double zeta, factor;

        for (i = 0; i < C.lpc.p; i++) { z[i] = GlblOpts.safeDivPos(lambda[i], C.lpc.w[i]); }

        cone_start = C.lpc.p;
        for (l = 0; l < C.nsoc; l++) {
            zeta = 0;
            for (i = 1; i < C.soc[l].p; i++) { zeta += C.soc[l].q[i - 1] * lambda[cone_start + i]; }

            factor = -lambda[cone_start] + GlblOpts.safeDivPos(zeta, (1 + C.soc[l].a));

            z[cone_start] = GlblOpts.safeDivPos((C.soc[l].a * lambda[cone_start] - zeta), C.soc[l].eta);
            for (i = 1; i < C.soc[l].p; i++) {
                j = cone_start + i;
                z[j] = GlblOpts.safeDivPos((lambda[j] + factor * C.soc[l].q[i - 1]), C.soc[l].eta);
            }

            cone_start += C.soc[l].p;
        }
    }

    public static double conicProduct(double[] u, double[] v, Cone C, double[] w) {
        int i, j, k = 0, cone_start, conesize;
        double u0, v0, mu;

        mu = 0;

        for (i = 0; i < C.lpc.p; i++) {
            w[k] = u[i] * v[i];
            mu += w[k] < 0 ? -w[k] : w[k];
            k++;
        }

        cone_start = C.lpc.p;
        for (i = 0; i < C.nsoc; i++) {
            conesize = C.soc[i].p;
            u0 = u[cone_start];
            v0 = v[cone_start];
            w[k] = SPLA.eddot(conesize, u, v, cone_start, cone_start);
            mu += w[k] < 0 ? -w[k] : w[k];
            k++;
            for (j = 1; j < conesize; j++) { w[k++] = u0 * v[cone_start + j] + v0 * u[cone_start + j]; }
            cone_start += conesize;
        }

        return mu;
    }

    public static void conicDivision(double[] u, double[] w, Cone C, double[] v) {
        int i, j, k, cone_start, conesize;
        double rho, zeta, u0, w0, factor, temp;

        for (i = 0; i < C.lpc.p; i++) { v[i] = GlblOpts.safeDivPos(w[i], u[i]); }

        cone_start = C.lpc.p;
        for (i = 0; i < C.nsoc; i++) {
            conesize = C.soc[i].p;
            u0 = u[cone_start];
            w0 = w[cone_start];
            rho = u0 * u0;
            zeta = 0;
            for (j = 1; j < conesize; j++) {
                k = cone_start + j;
                rho -= u[k] * u[k];
                zeta += u[k] * w[k];
            }
            temp = GlblOpts.safeDivPos(zeta, u0) - w0;
            factor = GlblOpts.safeDivPos(temp, rho);
            temp = u0 * w0 - zeta;
            v[cone_start] = GlblOpts.safeDivPos(temp, rho);
            for (j = 1; j < conesize; j++) {
                k = cone_start + j;
                v[cone_start + j] = factor * u[k] + GlblOpts.safeDivPos(w[k], u0);
            }
            cone_start += C.soc[i].p;
        }
    }

    public static void getSOCDetails(SOCone soc, int[] conesizeOut, double[] eta_squareOut, double[] d1Out, double[] u0Out, double[] u1Out, double[] v1Out, double[][] qOut) {
        conesizeOut[0] = soc.p;
        eta_squareOut[0] = soc.eta_square;
        d1Out[0] = soc.d1;
        u0Out[0] = soc.u0;
        u1Out[0] = soc.u1;
        v1Out[0] = soc.v1;
        qOut[0] = soc.q;
    }

    public static void unstretch(int n, int p, Cone C, int[] Pinv, double[] Px, double[] dx, double[] dy, double[] dz) {
        int i, j = 0, k = 0, l;

        for (i = 0; i < n; i++) { dx[i] = Px[Pinv[k++]]; }
        for (i = 0; i < p; i++) { dy[i] = Px[Pinv[k++]]; }
        for (i = 0; i < C.lpc.p; i++) { dz[j++] = Px[Pinv[k++]]; }
        for (l = 0; l < C.nsoc; l++) {
            for (i = 0; i < C.soc[l].p; i++) { dz[j++] = Px[Pinv[k++]]; }
            k += 2;
        }
        for (l = 0; l < C.nexc; l++) {
            for (i = 0; i < 3; i++) { dz[j++] = Px[Pinv[k++]]; }
        }
    }

    // Helper: socres with offset
    static double socres(double[] u, int p, int offset) {
        double res = u[offset] * u[offset];
        for (int i = 1; i < p; i++) { res -= u[offset + i] * u[offset + i]; }
        return res;
    }

    // Helper: eddot with offset
    static double eddot(int n, double[] x, double[] y, int xOffset, int yOffset) {
        double z = 0;
        for (int i = 0; i < n; i++) { z += x[xOffset + i] * y[yOffset + i]; }
        return z;
    }

    // Helper: scale with offset
    static void scale(double[] z, int zOffset, Cone C, double[] lambda) {
        int i, j, l, cone_start;
        double zeta, factor;

        for (i = 0; i < C.lpc.p; i++) { lambda[i] = C.lpc.w[i] * z[zOffset + i]; }

        cone_start = C.lpc.p;
        for (l = 0; l < C.nsoc; l++) {
            zeta = 0;
            for (i = 1; i < C.soc[l].p; i++) { zeta += C.soc[l].q[i - 1] * z[zOffset + cone_start + i]; }

            factor = z[zOffset + cone_start] + GlblOpts.safeDivPos(zeta, (1 + C.soc[l].a));

            lambda[cone_start] = C.soc[l].eta * (C.soc[l].a * z[zOffset + cone_start] + zeta);
            for (i = 1; i < C.soc[l].p; i++) {
                j = cone_start + i;
                lambda[j] = C.soc[l].eta * (z[zOffset + j] + factor * C.soc[l].q[i - 1]);
            }

            cone_start += C.soc[l].p;
        }
    }

    // Helper: evalExpHessian with offset
    static void evalExpHessian(double[] w, int wOffset, double[] v, double mu) {
        double x = w[wOffset];
        double y = w[wOffset + 1];
        double z_ = w[wOffset + 2];
        double l = Math.log(-y / x);
        double r = -x * l - x + z_;
        v[0] = mu * ((r * r - x * r + l * l * x * x) / (r * x * x * r));
        v[1] = mu * ((z_ - x) / (r * r * y));
        v[2] = mu * ((r * r - x * r + x * x) / (r * r * y * y));
        v[3] = mu * (-l / (r * r));
        v[4] = mu * (-x / (r * r * y));
        v[5] = mu * (1 / (r * r));
    }

    // Helper: evalExpGradient with offset
    static void evalExpGradient(double[] w, int wOffset, double[] g) {
        double x = w[wOffset];
        double y = w[wOffset + 1];
        double z_ = w[wOffset + 2];
        double l = Math.log(-y / x);
        double r = -x * l - x + z_;
        g[0] = (l * x - r) / (r * x);
        g[1] = (x - r) / (y * r);
        g[2] = -1 / r;
    }
}
