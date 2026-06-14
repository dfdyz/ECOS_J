package com.dfdyz.ecos_j;

public class KKTOps {

    public static int kkt_factor(KKT KKT, double eps, double delta) {
        int nd = LDL.ldl_numeric2(
            KKT.PKPt.n,
            KKT.PKPt.jc,
            KKT.PKPt.ir,
            KKT.PKPt.pr,
            KKT.L.jc,
            KKT.Parent,
            KKT.Sign,
            eps,
            delta,
            KKT.Lnz,
            KKT.L.ir,
            KKT.L.pr,
            KKT.D,
            KKT.work1,
            KKT.Pattern,
            KKT.Flag
        );
        return nd == KKT.PKPt.n ? GlblOpts.KKT_OK : GlblOpts.KKT_PROBLEM;
    }

    public static int kkt_factor(KKT KKT, double eps, double delta, double[] t1, double[] t2) {
        int nd = LDL.ldl_numeric2(
            KKT.PKPt.n,
            KKT.PKPt.jc,
            KKT.PKPt.ir,
            KKT.PKPt.pr,
            KKT.L.jc,
            KKT.Parent,
            KKT.Sign,
            eps,
            delta,
            KKT.Lnz,
            KKT.L.ir,
            KKT.L.pr,
            KKT.D,
            KKT.work1,
            KKT.Pattern,
            KKT.Flag,
            t1,
            t2
        );
        return nd == KKT.PKPt.n ? GlblOpts.KKT_OK : GlblOpts.KKT_PROBLEM;
    }

    public static int kkt_solve(KKT KKT, SpMat A, SpMat G, double[] Pb,
                                 double[] dx, double[] dy, double[] dz,
                                 int n, int p, int m, Cone C, int isinit, int nitref, int staticreg) {
        int MTILDE = m + 2 * C.nsoc;

        int i, k, l, j, kk, kItRef;
        int[] Pinv = KKT.Pinv;
        double[] Px = KKT.work1;
        double[] dPx = KKT.work2;
        double[] e = KKT.work3;
        double[] Pe = KKT.work4;
        double[] truez = KKT.work5;
        double[] Gdx = KKT.work6;
        double[] ex = new double[n + p + MTILDE];
        System.arraycopy(e, 0, ex, 0, n + p + MTILDE);
        double[] ey = new double[p];
        double[] ez = new double[MTILDE];
        for (int ii = 0; ii < p; ii++) ey[ii] = ex[n + ii];
        double bnorm = 1.0 + SPLA.norminf(Pb, n + p + MTILDE);
        double nex = 0;
        double ney = 0;
        double nez = 0;
        double nerr;
        double nerr_prev = GlblOpts.ECOS_NAN;
        double error_threshold = bnorm * GlblOpts.LINSYSACC;
        int nK = KKT.PKPt.n;

        LDL.ldl_lsolve2(nK, Pb, KKT.L.jc, KKT.L.ir, KKT.L.pr, Px);
        LDL.ldl_dsolve(nK, Px, KKT.D);
        LDL.ldl_ltsolve(nK, Px, KKT.L.jc, KKT.L.ir, KKT.L.pr);

        for (kItRef = 0; kItRef <= nitref; kItRef++) {
            ConeOps.unstretch(n, p, C, Pinv, Px, dx, dy, dz);

            k = 0; j = 0;
            if (staticreg > 0) {
                for (i = 0; i < n; i++) { ex[i] = Pb[Pinv[k++]] - GlblOpts.DELTASTAT * dx[i]; }
            } else {
                for (i = 0; i < n; i++) { ex[i] = Pb[Pinv[k++]]; }
            }
            if (A != null) SPLA.sparseMtVm(A, dy, ex, 0, 0);
            SPLA.sparseMtVm(G, dz, ex, 0, 0);
            nex = SPLA.norminf(ex, n);

            if (p > 0) {
                if (staticreg > 0) {
                    for (i = 0; i < p; i++) { ey[i] = Pb[Pinv[k++]] + GlblOpts.DELTASTAT * dy[i]; }
                } else {
                    for (i = 0; i < p; i++) { ey[i] = Pb[Pinv[k++]]; }
                }
                SPLA.sparseMV(A, dx, ey, -1, 0);
                ney = SPLA.norminf(ey, p);
            }

            kk = 0; j = 0;
            SPLA.sparseMV(G, dx, Gdx, 1, 1);
            if (staticreg > 0) {
                int dzoffset = 0;
                for (i = 0; i < C.lpc.p; i++) {
                    ez[kk++] = Pb[Pinv[k++]] - Gdx[j++] + GlblOpts.DELTASTAT * dz[dzoffset++];
                }
                for (l = 0; l < C.nsoc; l++) {
                    for (i = 0; i < C.soc[l].p; i++) {
                        ez[kk++] = i < (C.soc[l].p - 1) ?
                            Pb[Pinv[k++]] - Gdx[j++] + GlblOpts.DELTASTAT * dz[dzoffset++] :
                            Pb[Pinv[k++]] - Gdx[j++] - GlblOpts.DELTASTAT * dz[dzoffset++];
                    }
                    ez[kk] = 0;
                    ez[kk + 1] = 0;
                    k += 2;
                    kk += 2;
                }
                for (l = 0; l < C.nexc; l++) {
                    for (i = 0; i < 3; i++) {
                        ez[kk++] = Pb[Pinv[k++]] - Gdx[j++] + GlblOpts.DELTASTAT * dz[dzoffset++];
                    }
                }
            } else {
                for (i = 0; i < C.lpc.p; i++) {
                    ez[kk++] = Pb[Pinv[k++]] - Gdx[j++];
                }
                for (l = 0; l < C.nsoc; l++) {
                    for (i = 0; i < C.soc[l].p; i++) {
                        ez[kk++] = Pb[Pinv[k++]] - Gdx[j++];
                    }
                    ez[kk] = 0;
                    ez[kk + 1] = 0;
                    k += 2;
                    kk += 2;
                }
                for (l = 0; l < C.nexc; l++) {
                    for (i = 0; i < 3; i++) {
                        ez[kk++] = Pb[Pinv[k++]] - Gdx[j++];
                    }
                }
            }
            for (i = 0; i < MTILDE; i++) { truez[i] = Px[Pinv[n + p + i]]; }
            if (isinit == 0) {
                ConeOps.scale2add(truez, ez, C);
            } else {
                SPLA.vadd(MTILDE, truez, ez);
            }
            nez = SPLA.norminf(ez, MTILDE);

            nerr = Math.max(nex, nez);
            if (p > 0) { nerr = Math.max(nerr, ney); }

            if (kItRef > 0 && nerr > nerr_prev) {
                for (i = 0; i < nK; i++) { Px[i] -= dPx[i]; }
                kItRef--;
                break;
            }

            if (kItRef == nitref || (nerr < error_threshold) || (kItRef > 0 && nerr_prev < GlblOpts.IRERRFACT * nerr)) {
                break;
            }
            nerr_prev = nerr;

            for (i = 0; i < nK; i++) { Pe[Pinv[i]] = e[i]; }

            LDL.ldl_lsolve2(nK, Pe, KKT.L.jc, KKT.L.ir, KKT.L.pr, dPx);
            LDL.ldl_dsolve(nK, dPx, KKT.D);
            LDL.ldl_ltsolve(nK, dPx, KKT.L.jc, KKT.L.ir, KKT.L.pr);

            for (i = 0; i < nK; i++) { Px[i] += dPx[i]; }
        }

        ConeOps.unstretch(n, p, C, Pinv, Px, dx, dy, dz);

        return kItRef;
    }

    public static void kkt_update(SpMat PKP, int[] P, Cone C) {
        int i, j, k, conesize, conesize_m1;
        double eta_square, d1, u0, u1, v1;
        double[] q;

        for (i = 0; i < C.lpc.p; i++) {
            PKP.pr[P[C.lpc.kkt_idx[i]]] = -C.lpc.v[i] - GlblOpts.DELTASTAT;
        }

        for (i = 0; i < C.nsoc; i++) {
            SOCone soc = C.soc[i];
            conesize = soc.p;
            conesize_m1 = conesize - 1;
            eta_square = soc.eta_square;
            d1 = soc.d1;
            u0 = soc.u0;
            u1 = soc.u1;
            v1 = soc.v1;
            q = soc.q;

            PKP.pr[P[soc.Didx[0]]] = -eta_square * d1 - GlblOpts.DELTASTAT;
            for (k = 1; k < conesize; k++) {
                PKP.pr[P[soc.Didx[k]]] = -eta_square - GlblOpts.DELTASTAT;
            }

            j = 1;
            for (k = 0; k < conesize_m1; k++) {
                PKP.pr[P[soc.Didx[conesize_m1] + j++]] = -eta_square * v1 * q[k];
            }
            PKP.pr[P[soc.Didx[conesize_m1] + j++]] = -eta_square;

            PKP.pr[P[soc.Didx[conesize_m1] + j++]] = -eta_square * u0;
            for (k = 0; k < conesize_m1; k++) {
                PKP.pr[P[soc.Didx[conesize_m1] + j++]] = -eta_square * u1 * q[k];
            }
            PKP.pr[P[soc.Didx[conesize_m1] + j++]] = +eta_square + GlblOpts.DELTASTAT;
        }

        for (i = 0; i < C.nexc; i++) {
            PKP.pr[P[C.expc[i].colstart[0]]] = -C.expc[i].v[0] - GlblOpts.DELTASTAT;
            PKP.pr[P[C.expc[i].colstart[1]]] = -C.expc[i].v[1];
            PKP.pr[P[C.expc[i].colstart[1] + 1]] = -C.expc[i].v[2] - GlblOpts.DELTASTAT;
            PKP.pr[P[C.expc[i].colstart[2]]] = -C.expc[i].v[3];
            PKP.pr[P[C.expc[i].colstart[2] + 1]] = -C.expc[i].v[4];
            PKP.pr[P[C.expc[i].colstart[2] + 2]] = -C.expc[i].v[5] - GlblOpts.DELTASTAT;
        }
    }

    public static void kkt_init(SpMat PKP, int[] P, Cone C) {
        int i, k, conesize, conesize_m1;

        for (i = 0; i < C.lpc.p; i++) {
            PKP.pr[P[C.lpc.kkt_idx[i]]] = -1.0;
        }

        for (i = 0; i < C.nsoc; i++) {
            SOCone soc = C.soc[i];
            conesize = soc.p;
            conesize_m1 = conesize - 1;

            PKP.pr[P[soc.Didx[0]]] = -1.0;
            for (k = 1; k < conesize; k++) {
                PKP.pr[P[soc.Didx[k]]] = -1.0;
            }

            int j = 1;
            for (k = 0; k < conesize_m1; k++) {
                PKP.pr[P[soc.Didx[conesize_m1] + j++]] = 0.0;
            }
            PKP.pr[P[soc.Didx[conesize_m1] + j++]] = -1.0;

            PKP.pr[P[soc.Didx[conesize_m1] + j++]] = 0.0;
            for (k = 0; k < conesize_m1; k++) {
                PKP.pr[P[soc.Didx[conesize_m1] + j++]] = 0.0;
            }
            PKP.pr[P[soc.Didx[conesize_m1] + j++]] = +1.0;
        }
    }
}
