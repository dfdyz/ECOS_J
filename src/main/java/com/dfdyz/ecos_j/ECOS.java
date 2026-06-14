package com.dfdyz.ecos_j;

public class ECOS {

    private static final double[] t1_arr = new double[1];
    private static final double[] t2_arr = new double[1];

    static int compareStatistics(Stats infoA, Stats infoB) {
        if (infoA.pinfres != GlblOpts.ECOS_NAN && infoA.kapovert > 1) {
            if (infoB.pinfres != GlblOpts.ECOS_NAN) {
                if ((infoA.gap > 0 && infoB.gap > 0 && infoA.gap < infoB.gap) &&
                    (infoA.pinfres > 0 && infoA.pinfres < infoB.pres) &&
                    (infoA.mu > 0 && infoA.mu < infoB.mu)) {
                    return 1;
                } else {
                    return 0;
                }
            } else {
                if ((infoA.gap > 0 && infoB.gap > 0 && infoA.gap < infoB.gap) &&
                    (infoA.mu > 0 && infoA.mu < infoB.mu)) {
                    return 1;
                } else {
                    return 0;
                }
            }
        } else {
            if ((infoA.gap > 0 && infoB.gap > 0 && infoA.gap < infoB.gap) &&
                (infoA.pres > 0 && infoA.pres < infoB.pres) &&
                (infoA.dres > 0 && infoA.dres < infoB.dres) &&
                (infoA.kapovert > 0 && infoA.kapovert < infoB.kapovert) &&
                (infoA.mu > 0 && infoA.mu < infoB.mu)) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    static void saveIterateAsBest(PWork w) {
        int i;
        for (i = 0; i < w.n; i++) { w.best_x[i] = w.x[i]; }
        for (i = 0; i < w.p; i++) { w.best_y[i] = w.y[i]; }
        for (i = 0; i < w.m; i++) { w.best_z[i] = w.z[i]; }
        for (i = 0; i < w.m; i++) { w.best_s[i] = w.s[i]; }
        w.best_kap = w.kap;
        w.best_tau = w.tau;
        w.best_cx = w.cx;
        w.best_by = w.by;
        w.best_hz = w.hz;
        w.best_info.pcost = w.info.pcost;
        w.best_info.dcost = w.info.dcost;
        w.best_info.pres = w.info.pres;
        w.best_info.dres = w.info.dres;
        w.best_info.pinfres = w.info.pinfres;
        w.best_info.dinfres = w.info.dinfres;
        w.best_info.gap = w.info.gap;
        w.best_info.relgap = w.info.relgap;
        w.best_info.mu = w.info.mu;
        w.best_info.kapovert = w.info.kapovert;
        w.best_info.iter = w.info.iter;
    }

    static void restoreBestIterate(PWork w) {
        int i;
        for (i = 0; i < w.n; i++) { w.x[i] = w.best_x[i]; }
        for (i = 0; i < w.p; i++) { w.y[i] = w.best_y[i]; }
        for (i = 0; i < w.m; i++) { w.z[i] = w.best_z[i]; }
        for (i = 0; i < w.m; i++) { w.s[i] = w.best_s[i]; }
        w.kap = w.best_kap;
        w.tau = w.best_tau;
        w.cx = w.best_cx;
        w.by = w.best_by;
        w.hz = w.best_hz;
        w.info.pcost = w.best_info.pcost;
        w.info.dcost = w.best_info.dcost;
        w.info.pres = w.best_info.pres;
        w.info.dres = w.best_info.dres;
        w.info.pinfres = w.best_info.pinfres;
        w.info.dinfres = w.best_info.dinfres;
        w.info.gap = w.best_info.gap;
        w.info.relgap = w.best_info.relgap;
        w.info.mu = w.best_info.mu;
        w.info.kapovert = w.best_info.kapovert;
    }

    static int checkExitConditions(PWork w, int mode) {
        double feastol, abstol, reltol;

        if (mode == 0) {
            feastol = w.stgs.feastol;
            abstol = w.stgs.abstol;
            reltol = w.stgs.reltol;
        } else {
            feastol = w.stgs.feastol_inacc;
            abstol = w.stgs.abstol_inacc;
            reltol = w.stgs.reltol_inacc;
        }

        if ((-w.cx > 0 || -w.by - w.hz >= -abstol) &&
            (w.info.pres < feastol && w.info.dres < feastol) &&
            (w.info.gap < abstol || w.info.relgap < reltol)) {
            if (GlblOpts.PRINTLEVEL > 0) {
                if (w.stgs.verbose != 0) {
                    if (mode == 0) {
                        GlblOpts.printText("\nOPTIMAL (within feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).",
                            Math.max(w.info.dres, w.info.pres), w.info.relgap, w.info.gap);
                    } else {
                        GlblOpts.printText("\nClose to OPTIMAL (within feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).",
                            Math.max(w.info.dres, w.info.pres), w.info.relgap, w.info.gap);
                    }
                }
            }
            w.info.pinf = 0;
            w.info.dinf = 0;
            return GlblOpts.ECOS_OPTIMAL + mode;
        } else if ((w.info.dinfres != GlblOpts.ECOS_NAN) && (w.info.dinfres < feastol) && (w.tau < w.kap)) {
            if (GlblOpts.PRINTLEVEL > 0) {
                if (w.stgs.verbose != 0) {
                    if (mode == 0) {
                        GlblOpts.printText("\nUNBOUNDED (within feastol=%3.1e).", w.info.dinfres);
                    } else {
                        GlblOpts.printText("\nClose to UNBOUNDED (within feastol=%3.1e).", w.info.dinfres);
                    }
                }
            }
            w.info.pinf = 0;
            w.info.dinf = 1;
            return GlblOpts.ECOS_DINF + mode;
        } else if (((w.info.pinfres != GlblOpts.ECOS_NAN && w.info.pinfres < feastol) && (w.tau < w.kap)) ||
                   (w.tau < w.stgs.feastol && w.kap < w.stgs.feastol && w.info.pinfres < w.stgs.feastol)) {
            if (GlblOpts.PRINTLEVEL > 0) {
                if (w.stgs.verbose != 0) {
                    if (mode == 0) {
                        GlblOpts.printText("\nPRIMAL INFEASIBLE (within feastol=%3.1e).", w.info.pinfres);
                    } else {
                        GlblOpts.printText("\nClose to PRIMAL INFEASIBLE (within feastol=%3.1e).", w.info.pinfres);
                    }
                }
            }
            w.info.pinf = 1;
            w.info.dinf = 0;
            return GlblOpts.ECOS_PINF + mode;
        } else {
            return GlblOpts.ECOS_NOT_CONVERGED_YET;
        }
    }

    static int init(PWork w) {
        int i, j, k, l, KKT_FACTOR_RETURN_CODE;
        int[] Pinv = w.KKT.Pinv;
        double rx, ry, rz;

        Timer tfactor = new Timer();
        Timer tkktsolve = new Timer();

        w.KKT.delta = w.stgs.delta;

        KKTOps.kkt_init(w.KKT.PKPt, w.KKT.PK, w.C);

        k = 0;
        j = 0;
        for (i = 0; i < w.n; i++) { w.KKT.RHS1[w.KKT.Pinv[k++]] = 0; }
        for (i = 0; i < w.p; i++) { w.KKT.RHS1[w.KKT.Pinv[k++]] = w.b[i]; }
        for (i = 0; i < w.C.lpc.p; i++) { w.KKT.RHS1[w.KKT.Pinv[k++]] = w.h[i]; j++; }
        for (l = 0; l < w.C.nsoc; l++) {
            for (i = 0; i < w.C.soc[l].p; i++) { w.KKT.RHS1[w.KKT.Pinv[k++]] = w.h[j++]; }
            w.KKT.RHS1[w.KKT.Pinv[k++]] = 0;
            w.KKT.RHS1[w.KKT.Pinv[k++]] = 0;
        }
        for (l = 0; l < w.C.nexc; l++) {
            for (i = 0; i < 3; i++) { w.KKT.RHS1[w.KKT.Pinv[k++]] = w.h[j++]; }
        }

        for (i = 0; i < w.n; i++) { w.KKT.RHS2[w.KKT.Pinv[i]] = -w.c[i]; }
        for (i = w.n; i < w.KKT.PKPt.n; i++) { w.KKT.RHS2[w.KKT.Pinv[i]] = 0; }

        rx = SPLA.norm2(w.c, w.n);
        w.resx0 = Math.max(1, rx);
        ry = SPLA.norm2(w.b, w.p);
        w.resy0 = Math.max(1, ry);
        rz = SPLA.norm2(w.h, w.m);
        w.resz0 = Math.max(1, rz);

        if (w.C.nexc == 0) {
            if (GlblOpts.PROFILING > 1) {
                TimerOps.tic(tfactor);
                KKT_FACTOR_RETURN_CODE = KKTOps.kkt_factor(w.KKT, w.stgs.eps, w.stgs.delta, t1_arr, t2_arr);
                w.info.tfactor += TimerOps.toc(tfactor);
                w.info.tfactor_t1 = t1_arr[0];
                w.info.tfactor_t2 = t2_arr[0];
            } else {
                KKT_FACTOR_RETURN_CODE = KKTOps.kkt_factor(w.KKT, w.stgs.eps, w.stgs.delta);
            }

            if (KKT_FACTOR_RETURN_CODE != GlblOpts.KKT_OK) {
                if (GlblOpts.PRINTLEVEL > 0) {
                    if (w.stgs.verbose != 0) GlblOpts.printText("\nProblem in factoring KKT system, aborting.");
                }
                return GlblOpts.ECOS_FATAL;
            }

            if (GlblOpts.PROFILING > 1) TimerOps.tic(tkktsolve);
            w.info.nitref1 = KKTOps.kkt_solve(w.KKT, w.A, w.G, w.KKT.RHS1, w.KKT.dx1, w.KKT.dy1, w.KKT.dz1, w.n, w.p, w.m, w.C, 1, w.stgs.nitref, w.stgs.staticreg);
            if (GlblOpts.PROFILING > 1) w.info.tkktsolve += TimerOps.toc(tkktsolve);

            for (i = 0; i < w.n; i++) { w.x[i] = w.KKT.dx1[i]; }
            for (i = 0; i < w.m; i++) { w.KKT.work1[i] = -w.KKT.dz1[i]; }
            ConeOps.bring2cone(w.C, w.KKT.work1, w.s);

            if (GlblOpts.PROFILING > 1) TimerOps.tic(tkktsolve);
            w.info.nitref2 = KKTOps.kkt_solve(w.KKT, w.A, w.G, w.KKT.RHS2, w.KKT.dx2, w.KKT.dy2, w.KKT.dz2, w.n, w.p, w.m, w.C, 1, w.stgs.nitref, w.stgs.staticreg);
            if (GlblOpts.PROFILING > 1) w.info.tkktsolve += TimerOps.toc(tkktsolve);

            for (i = 0; i < w.p; i++) { w.y[i] = w.KKT.dy2[i]; }
            ConeOps.bring2cone(w.C, w.KKT.dz2, w.z);
        } else {
            ConeOps.unitInitialization(w.C, w.s, w.z, 1.0);
            for (i = 0; i < w.p; i++) { w.y[i] = 0.0; }
            for (i = 0; i < w.n; i++) { w.x[i] = 0.0; }
            w.info.nitref1 = 0;
            w.info.nitref2 = 0;
        }

        for (i = 0; i < w.n; i++) { w.KKT.RHS1[Pinv[i]] = -w.c[i]; }

        w.kap = 1.0;
        w.tau = 1.0;
        w.info.step = 0;
        w.info.step_aff = 0;
        w.info.dinf = 0;
        w.info.pinf = 0;

        return 0;
    }

    static void computeResiduals(PWork w) {
        if (w.p > 0) {
            SPLA.sparseMtVm(w.A, w.y, w.rx, 1, 0);
            SPLA.sparseMtVm(w.G, w.z, w.rx, 0, 0);
        } else {
            SPLA.sparseMtVm(w.G, w.z, w.rx, 1, 0);
        }
        w.hresx = SPLA.norm2(w.rx, w.n);
        SPLA.vsubscale(w.n, w.tau, w.c, w.rx);

        if (w.p > 0) {
            SPLA.sparseMV(w.A, w.x, w.ry, 1, 1);
            w.hresy = SPLA.norm2(w.ry, w.p);
            SPLA.vsubscale(w.p, w.tau, w.b, w.ry);
        } else {
            w.hresy = 0;
            w.ry = null;
        }

        SPLA.sparseMV(w.G, w.x, w.rz, 1, 1);
        SPLA.vadd(w.m, w.s, w.rz);
        w.hresz = SPLA.norm2(w.rz, w.m);
        SPLA.vsubscale(w.m, w.tau, w.h, w.rz);

        w.cx = SPLA.eddot(w.n, w.c, w.x);
        w.by = w.p > 0 ? SPLA.eddot(w.p, w.b, w.y) : 0.0;
        w.hz = SPLA.eddot(w.m, w.h, w.z);
        w.rt = w.kap + w.cx + w.by + w.hz;

        w.nx = SPLA.norm2(w.x, w.n);
        w.ny = SPLA.norm2(w.y, w.p);
        w.ns = SPLA.norm2(w.s, w.m);
        w.nz = SPLA.norm2(w.z, w.m);
    }

    static void updateStatistics(PWork w) {
        double nry, nrz;
        Stats info = w.info;

        info.gap = SPLA.eddot(w.m, w.s, w.z);
        info.mu = (info.gap + w.kap * w.tau) / (w.D + 1);

        info.kapovert = w.kap / w.tau;
        info.pcost = w.cx / w.tau;
        info.dcost = -(w.hz + w.by) / w.tau;

        if (info.pcost < 0) { info.relgap = info.gap / (-info.pcost); }
        else if (info.dcost > 0) { info.relgap = info.gap / info.dcost; }
        else info.relgap = GlblOpts.ECOS_NAN;

        nry = w.p > 0 ? SPLA.norm2(w.ry, w.p) / Math.max(w.resy0 + w.nx, 1) : 0.0;
        nrz = SPLA.norm2(w.rz, w.m) / Math.max(w.resz0 + w.nx + w.ns, 1);
        info.pres = Math.max(nry, nrz) / w.tau;
        info.dres = SPLA.norm2(w.rx, w.n) / Math.max(w.resx0 + w.ny + w.nz, 1) / w.tau;

        info.pinfres = (w.hz + w.by) / Math.max(w.ny + w.nz, 1) < -w.stgs.reltol ? w.hresx / Math.max(w.ny + w.nz, 1) : GlblOpts.ECOS_NAN;
        info.dinfres = w.cx / Math.max(w.nx, 1) < -w.stgs.reltol ? Math.max(w.hresy / Math.max(w.nx, 1), w.hresz / Math.max(w.nx + w.ns, 1)) : GlblOpts.ECOS_NAN;
    }

    static void printProgress(Stats info) {
        if (info.iter == 0) {
            if (GlblOpts.PRINTLEVEL == 2) {
                GlblOpts.printText("\nECOS 2.0.10 - (C) embotech GmbH, Zurich Switzerland, 2012-15. Web: www.embotech.com/ECOS\n");
                GlblOpts.printText("\n");
            }
            GlblOpts.printText("It     pcost       dcost      gap   pres   dres    k/t    mu     step   sigma     IR    |   BT\n");
            GlblOpts.printText("%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e    ---    ---   %2d %2d  - |  -  - \n",
                (int) info.iter, info.pcost, info.dcost, info.gap, info.pres, info.dres, info.kapovert, info.mu,
                (int) info.nitref1, (int) info.nitref2);
        } else {
            GlblOpts.printText("%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e  %6.4f  %2.0e  %2d %2d %2d | %2d %2d\n",
                (int) info.iter, info.pcost, info.dcost, info.gap, info.pres, info.dres, info.kapovert, info.mu,
                info.step, info.sigma,
                (int) info.nitref1, (int) info.nitref2, (int) info.nitref3,
                (int) info.affBack, (int) info.cmbBack);
        }
    }

    static void RHS_affine(PWork w) {
        double[] RHS = w.KKT.RHS2;
        int n = w.n;
        int p = w.p;
        int i, j, k, l;
        int[] Pinv = w.KKT.Pinv;

        j = 0;
        for (i = 0; i < n; i++) { RHS[Pinv[j++]] = w.rx[i]; }
        for (i = 0; i < p; i++) { RHS[Pinv[j++]] = -w.ry[i]; }
        for (i = 0; i < w.C.lpc.p; i++) { RHS[Pinv[j++]] = w.s[i] - w.rz[i]; }
        k = w.C.lpc.p;
        for (l = 0; l < w.C.nsoc; l++) {
            for (i = 0; i < w.C.soc[l].p; i++) {
                RHS[Pinv[j++]] = w.s[k] - w.rz[k]; k++;
            }
            RHS[Pinv[j++]] = 0;
            RHS[Pinv[j++]] = 0;
        }
        for (l = 0; l < w.C.nexc; l++) {
            for (i = 0; i < 3; i++) {
                RHS[Pinv[j++]] = w.s[k] - w.rz[k]; k++;
            }
        }
    }

    static void RHS_combined(PWork w) {
        double[] ds1 = w.KKT.work1;
        double[] ds2 = w.KKT.work2;
        int i, j, k, l;
        double sigmamu = w.info.sigma * w.info.mu;
        double one_minus_sigma = 1.0 - w.info.sigma;
        int[] Pinv = w.KKT.Pinv;

        ConeOps.conicProduct(w.lambda, w.lambda, w.C, ds1);
        ConeOps.conicProduct(w.dsaff_by_W, w.W_times_dzaff, w.C, ds2);
        for (i = 0; i < w.C.lpc.p; i++) { ds1[i] += ds2[i] - sigmamu; }
        k = w.C.lpc.p;
        for (i = 0; i < w.C.nsoc; i++) {
            ds1[k] += ds2[k] - sigmamu; k++;
            for (j = 1; j < w.C.soc[i].p; j++) { ds1[k] += ds2[k]; k++; }
        }

        ConeOps.conicDivision(w.lambda, ds1, w.C, w.dsaff_by_W);
        ConeOps.scale(w.dsaff_by_W, w.C, ds1);

        j = 0;
        for (i = 0; i < w.n; i++) { w.KKT.RHS2[Pinv[j++]] *= one_minus_sigma; }
        for (i = 0; i < w.p; i++) { w.KKT.RHS2[Pinv[j++]] *= one_minus_sigma; }
        for (i = 0; i < w.C.lpc.p; i++) { w.KKT.RHS2[Pinv[j++]] = -one_minus_sigma * w.rz[i] + ds1[i]; }
        k = w.C.lpc.p;
        for (l = 0; l < w.C.nsoc; l++) {
            for (i = 0; i < w.C.soc[l].p; i++) {
                w.KKT.RHS2[Pinv[j++]] = -one_minus_sigma * w.rz[k] + ds1[k];
                k++;
            }
            w.KKT.RHS2[Pinv[j++]] = 0;
            w.KKT.RHS2[Pinv[j++]] = 0;
        }

        k = w.C.fexv;
        for (l = 0; l < w.C.nexc; l++) {
            w.C.expc[l].g[0] = w.s[k] + sigmamu * w.C.expc[l].g[0];
            w.KKT.RHS2[Pinv[j++]] = -one_minus_sigma * w.rz[k] + w.C.expc[l].g[0];
            k++;
            w.C.expc[l].g[1] = w.s[k] + sigmamu * w.C.expc[l].g[1];
            w.KKT.RHS2[Pinv[j++]] = -one_minus_sigma * w.rz[k] + w.C.expc[l].g[1];
            k++;
            w.C.expc[l].g[2] = w.s[k] + sigmamu * w.C.expc[l].g[2];
            w.KKT.RHS2[Pinv[j++]] = -one_minus_sigma * w.rz[k] + w.C.expc[l].g[2];
            k++;
        }
    }

    static double expConeLineSearch(PWork w, double dtau, double dkappa, int affine) {
        double[] ws = w.KKT.work1;
        double[] wz = w.KKT.work2;
        double[] s = w.s;
        double[] z = w.z;
        double[] ds = w.dsaff;
        double[] dz = w.KKT.dz2;
        double tau = w.tau;
        double kap = w.kap;
        double gamma = w.stgs.gamma;

        int bk_iter, j;
        double mui, mu;
        double alpha;
        double[] centrality = new double[1];
        double barrier;
        double cent_constant = w.D + 1;
        int fc = w.C.fexv;

        if (affine == 1) { alpha = w.info.step_aff; }
        else { alpha = w.info.step; }

        w.info.pob = 0;
        w.info.cb = 0;
        w.info.cob = 0;
        w.info.pb = 0;
        w.info.db = 0;

        for (bk_iter = 0; bk_iter < w.stgs.max_bk_iter; bk_iter++) {
            mu = 0.0;
            for (j = 0; j < w.m; j++) {
                ws[j] = s[j] + alpha * ds[j];
                wz[j] = z[j] + alpha * dz[j];
                mu += ws[j] * wz[j];
            }
            mu = mu + (tau + alpha * dtau) * (kap + alpha * dkappa);
            mu = mu / (w.D + 1);

            if (ExpConeOps.evalExpDualFeas(wz, fc, w.C.nexc) == 1) {
                if (ExpConeOps.evalExpPrimalFeas(ws, fc, w.C.nexc) == 1) {
                    j = w.C.fexv;
                    mui = ws[j] * wz[j] + ws[j + 1] * wz[j + 1] + ws[j + 2] * wz[j + 2];
                    mui /= 3.0;
                    while (mui > GlblOpts.MIN_DISTANCE * mu && j < w.m - 2) {
                        j += 3;
                        if (j < w.m) {
                            mui = ws[j] * wz[j] + ws[j + 1] * wz[j + 1] + ws[j + 2] * wz[j + 2];
                            mui /= 3.0;
                        }
                    }
                    if (j == w.m) {
                        barrier = ExpConeOps.evalBarrierValue(ws, wz, w.C.fexv, w.C.nexc);
                        barrier += ConeOps.evalSymmetricBarrierValue(ws, wz, tau + alpha * dtau, kap + alpha * dkappa, w.C, w.D);
                        centrality[0] = barrier + cent_constant * Math.log(mu) + cent_constant;

                        if (centrality[0] < w.stgs.centrality) {
                            w.info.centrality = centrality[0];
                            return alpha * gamma;
                        } else {
                            w.info.cb++;
                        }
                    } else {
                        w.info.cob++;
                    }
                } else {
                    w.info.pb++;
                }
            } else {
                w.info.db++;
            }
            alpha = alpha * w.stgs.bk_scale;
        }
        return -1;
    }

    static double lineSearch(double[] lambda, double[] ds, double[] dz, double tau, double dtau, double kap, double dkap, Cone C, KKT KKT) {
        int i, j, cone_start, conesize;
        double rhomin, sigmamin, alpha, lknorm2, lknorm, lknorminv, rhonorm, sigmanorm, conic_step, temp;
        double lkbar_times_dsk, lkbar_times_dzk, factor;
        double[] lk;
        double[] dsk;
        double[] dzk;
        double[] lkbar = KKT.work1;
        double[] rho = KKT.work2;
        double minus_tau_by_dtau = -tau / dtau;
        double minus_kap_by_dkap = -kap / dkap;

        if (C.lpc.p > 0) {
            rhomin = ds[0] / lambda[0];
            sigmamin = dz[0] / lambda[0];
            for (i = 1; i < C.lpc.p; i++) {
                rho[0] = ds[i] / lambda[i]; if (rho[0] < rhomin) { rhomin = rho[0]; }
                rho[0] = dz[i] / lambda[i]; if (rho[0] < sigmamin) { sigmamin = rho[0]; }
            }
            if (-sigmamin > -rhomin) {
                alpha = sigmamin < 0 ? 1.0 / (-sigmamin) : 1.0 / GlblOpts.EPS;
            } else {
                alpha = rhomin < 0 ? 1.0 / (-rhomin) : 1.0 / GlblOpts.EPS;
            }
        } else {
            alpha = 10;
        }

        if (minus_tau_by_dtau > 0 && minus_tau_by_dtau < alpha) { alpha = minus_tau_by_dtau; }
        if (minus_kap_by_dkap > 0 && minus_kap_by_dkap < alpha) { alpha = minus_kap_by_dkap; }

        cone_start = C.lpc.p;
        for (i = 0; i < C.nsoc; i++) {
            conesize = C.soc[i].p;
            lk = lambda;
            dsk = ds;
            dzk = dz;
            int lkOff = cone_start;
            int dskOff = cone_start;
            int dzkOff = cone_start;

            lknorm2 = lk[lkOff] * lk[lkOff] - SPLA.eddot(conesize - 1, lk, lk, lkOff + 1, lkOff + 1);
            if (lknorm2 <= 0.0) { cone_start += C.soc[i].p; continue; }

            lknorm = Math.sqrt(lknorm2);
            for (j = 0; j < conesize; j++) { lkbar[j] = lk[lkOff + j] / lknorm; }
            lknorminv = 1.0 / lknorm;

            lkbar_times_dsk = lkbar[0] * dsk[dskOff];
            for (j = 1; j < conesize; j++) { lkbar_times_dsk -= lkbar[j] * dsk[dskOff + j]; }
            lkbar_times_dzk = lkbar[0] * dzk[dzkOff];
            for (j = 1; j < conesize; j++) { lkbar_times_dzk -= lkbar[j] * dzk[dzkOff + j]; }

            rho[0] = lknorminv * lkbar_times_dsk;
            factor = (lkbar_times_dsk + dsk[dskOff]) / (lkbar[0] + 1);
            for (j = 1; j < conesize; j++) { rho[j] = lknorminv * (dsk[dskOff + j] - factor * lkbar[j]); }
            rhonorm = SPLA.norm2(rho, 1, conesize - 1) - rho[0];

            rho[0] = lknorminv * lkbar_times_dzk;
            factor = (lkbar_times_dzk + dzk[dzkOff]) / (lkbar[0] + 1);
            for (j = 1; j < conesize; j++) { rho[j] = lknorminv * (dzk[dzkOff + j] - factor * lkbar[j]); }
            sigmanorm = SPLA.norm2(rho, 1, conesize - 1) - rho[0];

            conic_step = 0;
            if (rhonorm > conic_step) { conic_step = rhonorm; }
            if (sigmanorm > conic_step) { conic_step = sigmanorm; }
            if (conic_step != 0) {
                temp = 1.0 / conic_step;
                if (temp < alpha) { alpha = temp; }
            }

            cone_start += C.soc[i].p;
        }

        if (alpha > GlblOpts.STEPMAX) alpha = GlblOpts.STEPMAX;
        if (alpha < GlblOpts.STEPMIN) alpha = GlblOpts.STEPMIN;

        return alpha;
    }

    static void backscale(PWork w) {
        int i;
        if (w.stgs.equilibrate > 0) {
            for (i = 0; i < w.n; i++) { w.x[i] /= (w.xequil[i] * w.tau); }
            for (i = 0; i < w.p; i++) { w.y[i] /= (w.Aequil[i] * w.tau); }
            for (i = 0; i < w.m; i++) { w.z[i] /= (w.Gequil[i] * w.tau); }
            for (i = 0; i < w.m; i++) { w.s[i] *= (w.Gequil[i] / w.tau); }
            for (i = 0; i < w.n; i++) { w.c[i] *= w.xequil[i]; }
        } else {
            for (i = 0; i < w.n; i++) { w.x[i] /= w.tau; }
            for (i = 0; i < w.p; i++) { w.y[i] /= w.tau; }
            for (i = 0; i < w.m; i++) { w.z[i] /= w.tau; }
            for (i = 0; i < w.m; i++) { w.s[i] /= w.tau; }
        }
    }

    public static int ECOS_solve(PWork w) {
        int i, initcode, KKT_FACTOR_RETURN_CODE;
        double dtau_denom, dtauaff, dkapaff, sigma, dtau, dkap, bkap;
        int exitcode = GlblOpts.ECOS_FATAL, interrupted = 0;
        double pres_prev = GlblOpts.ECOS_NAN;
        int fc = w.C.fexv;
        int k_outer;

        Timer tsolve = new Timer();
        Timer tfactor = new Timer();
        Timer tkktsolve = new Timer();

        if (w.stgs.equilibrate > 0) {
            for (i = 0; i < w.n; i++) { w.c[i] /= w.xequil[i]; }
        }

        if (GlblOpts.PROFILING > 0) TimerOps.tic(tsolve);

        initcode = init(w);
        if (initcode == GlblOpts.ECOS_FATAL) {
            if (GlblOpts.PRINTLEVEL > 0) {
                if (w.stgs.verbose != 0) GlblOpts.printText("\nFatal error during initialization, aborting.");
            }
            return GlblOpts.ECOS_FATAL;
        }

        for (w.info.iter = 0; w.info.iter <= w.stgs.maxit; w.info.iter++) {
            computeResiduals(w);
            updateStatistics(w);

            if (GlblOpts.PRINTLEVEL > 1) {
                if (w.stgs.verbose != 0) printProgress(w.info);
            }

            if (w.info.iter > 0 && (w.info.pres > GlblOpts.SAFEGUARD * pres_prev || w.info.gap < 0)) {
                if (GlblOpts.PRINTLEVEL > 1) {
                    if (w.stgs.verbose != 0)
                        GlblOpts.printText("Unreliable search direction detected, recovering best iterate (%d) and stopping.\n", (int) w.best_info.iter);
                }
                restoreBestIterate(w);
                exitcode = checkExitConditions(w, GlblOpts.ECOS_INACC_OFFSET);
                if (exitcode == GlblOpts.ECOS_NOT_CONVERGED_YET) {
                    exitcode = GlblOpts.ECOS_NUMERICS;
                    if (GlblOpts.PRINTLEVEL > 0) {
                        if (w.stgs.verbose != 0)
                            GlblOpts.printText("\nNUMERICAL PROBLEMS (reached feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).",
                                Math.max(w.info.dres, w.info.pres), w.info.relgap, w.info.gap);
                    }
                    break;
                } else {
                    break;
                }
            }
            pres_prev = w.info.pres;

            exitcode = checkExitConditions(w, 0);
            if (exitcode == GlblOpts.ECOS_NOT_CONVERGED_YET) {
                if (w.info.iter > 0 && w.info.step == GlblOpts.STEPMIN * GlblOpts.GAMMA) {
                    if (GlblOpts.PRINTLEVEL > 0) {
                        if (w.stgs.verbose != 0)
                            GlblOpts.printText("No further progress possible, recovering best iterate (%d) and stopping.", (int) w.best_info.iter);
                    }
                    restoreBestIterate(w);
                    exitcode = checkExitConditions(w, GlblOpts.ECOS_INACC_OFFSET);
                    if (exitcode == GlblOpts.ECOS_NOT_CONVERGED_YET) {
                        exitcode = GlblOpts.ECOS_NUMERICS;
                        if (GlblOpts.PRINTLEVEL > 0) {
                            if (w.stgs.verbose != 0)
                                GlblOpts.printText("\nNUMERICAL PROBLEMS (reached feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).",
                                    Math.max(w.info.dres, w.info.pres), w.info.relgap, w.info.gap);
                        }
                    }
                    break;
                } else if (interrupted != 0 || w.info.iter == w.stgs.maxit) {
                    if (compareStatistics(w.info, w.best_info) != 0) {
                    } else {
                        restoreBestIterate(w);
                    }
                    exitcode = checkExitConditions(w, GlblOpts.ECOS_INACC_OFFSET);
                    if (exitcode == GlblOpts.ECOS_NOT_CONVERGED_YET) {
                        exitcode = interrupted != 0 ? GlblOpts.ECOS_SIGINT : GlblOpts.ECOS_MAXIT;
                    }
                    break;
                } else if (GlblOpts.isNan(w.info.pcost)) {
                    if (compareStatistics(w.info, w.best_info) != 0) {
                    } else {
                        restoreBestIterate(w);
                    }
                    exitcode = checkExitConditions(w, GlblOpts.ECOS_INACC_OFFSET);
                    if (exitcode == GlblOpts.ECOS_NOT_CONVERGED_YET) {
                        exitcode = GlblOpts.ECOS_NUMERICS;
                    }
                    break;
                }
            } else {
                break;
            }

            if (w.info.iter == 0) {
                saveIterateAsBest(w);
            } else if (compareStatistics(w.info, w.best_info) != 0) {
                saveIterateAsBest(w);
            }

            if (ConeOps.updateScalings(w.C, w.s, w.z, w.lambda, w.info.mu) == GlblOpts.OUTSIDE_CONE) {
                if (GlblOpts.PRINTLEVEL > 0) {
                    if (w.stgs.verbose != 0)
                        GlblOpts.printText("Slacks/multipliers leaving the cone, recovering best iterate (%d) and stopping.\n", (int) w.best_info.iter);
                }
                restoreBestIterate(w);
                exitcode = checkExitConditions(w, GlblOpts.ECOS_INACC_OFFSET);
                if (exitcode == GlblOpts.ECOS_NOT_CONVERGED_YET) {
                    if (GlblOpts.PRINTLEVEL > 0) {
                        if (w.stgs.verbose != 0)
                            GlblOpts.printText("\nNUMERICAL PROBLEMS (reached feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).",
                                Math.max(w.info.dres, w.info.pres), w.info.relgap, w.info.gap);
                    }
                    return GlblOpts.ECOS_OUTCONE;
                } else {
                    break;
                }
            }

            KKTOps.kkt_update(w.KKT.PKPt, w.KKT.PK, w.C);

            if (GlblOpts.PROFILING > 1) {
                TimerOps.tic(tfactor);
                KKT_FACTOR_RETURN_CODE = KKTOps.kkt_factor(w.KKT, w.stgs.eps, w.stgs.delta, t1_arr, t2_arr);
                w.info.tfactor += TimerOps.toc(tfactor);
                w.info.tfactor_t1 = t1_arr[0];
                w.info.tfactor_t2 = t2_arr[0];
            } else {
                KKT_FACTOR_RETURN_CODE = KKTOps.kkt_factor(w.KKT, w.stgs.eps, w.stgs.delta);
            }

            if (KKT_FACTOR_RETURN_CODE != GlblOpts.KKT_OK) {
                if (GlblOpts.PRINTLEVEL > 0) {
                    if (w.stgs.verbose != 0) GlblOpts.printText("\nProblem in factoring KKT system, aborting.");
                }
                return GlblOpts.ECOS_FATAL;
            }

            if (GlblOpts.PROFILING > 1) TimerOps.tic(tkktsolve);
            w.info.nitref1 = KKTOps.kkt_solve(w.KKT, w.A, w.G, w.KKT.RHS1, w.KKT.dx1, w.KKT.dy1, w.KKT.dz1, w.n, w.p, w.m, w.C, 0, w.stgs.nitref, w.stgs.staticreg);
            if (GlblOpts.PROFILING > 1) w.info.tkktsolve += TimerOps.toc(tkktsolve);

            RHS_affine(w);
            if (GlblOpts.PROFILING > 1) TimerOps.tic(tkktsolve);
            w.info.nitref2 = KKTOps.kkt_solve(w.KKT, w.A, w.G, w.KKT.RHS2, w.KKT.dx2, w.KKT.dy2, w.KKT.dz2, w.n, w.p, w.m, w.C, 0, w.stgs.nitref, w.stgs.staticreg);
            if (GlblOpts.PROFILING > 1) w.info.tkktsolve += TimerOps.toc(tkktsolve);

            dtau_denom = w.kap / w.tau - SPLA.eddot(w.n, w.c, w.KKT.dx1) - SPLA.eddot(w.p, w.b, w.KKT.dy1) - SPLA.eddot(w.m, w.h, w.KKT.dz1);

            dtauaff = (w.rt - w.kap + SPLA.eddot(w.n, w.c, w.KKT.dx2) + SPLA.eddot(w.p, w.b, w.KKT.dy2) + SPLA.eddot(w.m, w.h, w.KKT.dz2)) / dtau_denom;

            for (i = 0; i < w.m; i++) { w.KKT.dz2[i] = w.KKT.dz2[i] + dtauaff * w.KKT.dz1[i]; }
            ConeOps.scale(w.KKT.dz2, w.C, w.W_times_dzaff);

            for (i = 0; i < w.m; i++) { w.dsaff_by_W[i] = -w.W_times_dzaff[i] - w.lambda[i]; }

            dkapaff = -w.kap - w.kap / w.tau * dtauaff;

            w.info.step_aff = lineSearch(w.lambda, w.dsaff_by_W, w.W_times_dzaff, w.tau, dtauaff, w.kap, dkapaff, w.C, w.KKT);

            ConeOps.scale(w.dsaff_by_W, w.C, w.dsaff);

            for (i = fc; i < w.m; i++) { w.dsaff[i] = 0.0; }
            ExpConeOps.scaleToAddExpcone(w.dsaff, w.KKT.dz2, w.C.expc, w.C.nexc, fc);
            for (i = fc; i < w.m; i++) { w.dsaff[i] = -w.dsaff[i] - w.s[i]; }

            w.info.affBack = 0;
            if (w.C.nexc > 0) {
                w.info.step_aff = expConeLineSearch(w, dtauaff, dkapaff, 1);
                w.info.affBack += w.info.pob;
                w.info.affBack += w.info.cb;
                w.info.affBack += w.info.cob;
                w.info.affBack += w.info.pb;
                w.info.affBack += w.info.db;
                if (w.info.step_aff == -1) { w.info.step_aff = 0.0; }
            }

            sigma = 1.0 - w.info.step_aff;
            sigma = sigma * sigma * sigma;
            if (sigma > GlblOpts.SIGMAMAX) sigma = GlblOpts.SIGMAMAX;
            if (sigma < GlblOpts.SIGMAMIN) sigma = GlblOpts.SIGMAMIN;
            w.info.sigma = sigma;

            RHS_combined(w);
            if (GlblOpts.PROFILING > 1) TimerOps.tic(tkktsolve);
            w.info.nitref3 = KKTOps.kkt_solve(w.KKT, w.A, w.G, w.KKT.RHS2, w.KKT.dx2, w.KKT.dy2, w.KKT.dz2, w.n, w.p, w.m, w.C, 0, w.stgs.nitref, w.stgs.staticreg);
            if (GlblOpts.PROFILING > 1) w.info.tkktsolve += TimerOps.toc(tkktsolve);

            bkap = w.kap * w.tau + dkapaff * dtauaff - sigma * w.info.mu;

            dtau = ((1 - sigma) * w.rt - bkap / w.tau + SPLA.eddot(w.n, w.c, w.KKT.dx2) + SPLA.eddot(w.p, w.b, w.KKT.dy2) + SPLA.eddot(w.m, w.h, w.KKT.dz2)) / dtau_denom;

            for (i = 0; i < w.n; i++) { w.KKT.dx2[i] += dtau * w.KKT.dx1[i]; }
            for (i = 0; i < w.p; i++) { w.KKT.dy2[i] += dtau * w.KKT.dy1[i]; }
            for (i = 0; i < w.m; i++) { w.KKT.dz2[i] += dtau * w.KKT.dz1[i]; }

            ConeOps.scale(w.KKT.dz2, w.C, w.W_times_dzaff);
            for (i = 0; i < w.m; i++) { w.dsaff_by_W[i] = -(w.dsaff_by_W[i] + w.W_times_dzaff[i]); }

            dkap = -(bkap + w.kap * dtau) / w.tau;

            w.info.step = lineSearch(w.lambda, w.dsaff_by_W, w.W_times_dzaff, w.tau, dtau, w.kap, dkap, w.C, w.KKT) * w.stgs.gamma;

            ConeOps.scale(w.dsaff_by_W, w.C, w.dsaff);

            for (i = fc; i < w.m; i++) { w.dsaff[i] = 0.0; }
            ExpConeOps.scaleToAddExpcone(w.dsaff, w.KKT.dz2, w.C.expc, w.C.nexc, fc);

            k_outer = fc;
            for (i = 0; i < w.C.nexc; i++) {
                w.dsaff[k_outer] = -w.dsaff[k_outer] - w.C.expc[i].g[0]; k_outer++;
                w.dsaff[k_outer] = -w.dsaff[k_outer] - w.C.expc[i].g[1]; k_outer++;
                w.dsaff[k_outer] = -w.dsaff[k_outer] - w.C.expc[i].g[2]; k_outer++;
            }

            w.info.cmbBack = 0;
            if (w.C.nexc > 0) {
                w.info.step = expConeLineSearch(w, dtau, dkap, 0);
                w.info.cmbBack += w.info.cb;
                w.info.cmbBack += w.info.cob;
                w.info.cmbBack += w.info.pb;
                w.info.cmbBack += w.info.db;

                if (w.info.step == -1) {
                    restoreBestIterate(w);
                    exitcode = checkExitConditions(w, GlblOpts.ECOS_INACC_OFFSET);
                    if (exitcode == GlblOpts.ECOS_NOT_CONVERGED_YET) {
                        exitcode = GlblOpts.ECOS_NUMERICS;
                        if (GlblOpts.PRINTLEVEL > 0) {
                            if (w.stgs.verbose != 0)
                                GlblOpts.printText("\nNUMERICAL PROBLEMS (reached feastol=%3.1e, reltol=%3.1e, abstol=%3.1e).",
                                    Math.max(w.info.dres, w.info.pres), w.info.relgap, w.info.gap);
                        }
                        break;
                    } else {
                        break;
                    }
                }
            }

            for (i = 0; i < w.n; i++) { w.x[i] += w.info.step * w.KKT.dx2[i]; }
            for (i = 0; i < w.p; i++) { w.y[i] += w.info.step * w.KKT.dy2[i]; }
            for (i = 0; i < w.m; i++) { w.z[i] += w.info.step * w.KKT.dz2[i]; }
            for (i = 0; i < w.m; i++) { w.s[i] += w.info.step * w.dsaff[i]; }
            w.kap += w.info.step * dkap;
            w.tau += w.info.step * dtau;
        }

        backscale(w);

        if (GlblOpts.PROFILING > 0) w.info.tsolve = TimerOps.toc(tsolve);

        if (GlblOpts.PRINTLEVEL > 0) {
            if (GlblOpts.PROFILING > 0) {
                if (w.stgs.verbose != 0) GlblOpts.printText("\nRuntime: %f seconds.", w.info.tsetup + w.info.tsolve);
            }
            if (w.stgs.verbose != 0) GlblOpts.printText("\n\n");
        }

        return exitcode;
    }

    public static PWork ECOS_setup(int n, int m, int p, int l, int ncones, int[] q, int nexc,
                                    double[] Gpr, int[] Gjc, int[] Gir,
                                    double[] Apr, int[] Ajc, int[] Air,
                                    double[] c, double[] h, double[] b) {
        int i, cidx, conesize, lnz, amd_result, nK;
        int[] P, Pinv, Sign;
        PWork mywork;
        double[] Control = new double[AMD.AMD_CONTROL];
        double[] Info = new double[AMD.AMD_INFO];
        double[] Lpr;
        SpMat At, Gt, KU;
        int[] AtoAt, GtoGt, AttoK, GttoK;

        Timer tsetup = new Timer();
        Timer tcreatekkt = new Timer();
        Timer tmattranspose = new Timer();
        Timer tordering = new Timer();

        if (GlblOpts.PROFILING > 0) TimerOps.tic(tsetup);

        mywork = new PWork();
        mywork.n = n;
        mywork.m = m;
        mywork.p = p;
        mywork.D = l + ncones + 3 * nexc;

        mywork.x = new double[n];
        mywork.y = new double[p];
        mywork.z = new double[m];
        mywork.s = new double[m];
        mywork.lambda = new double[m];
        mywork.dsaff_by_W = new double[m];
        mywork.dsaff = new double[m];
        mywork.dzaff = new double[m];
        mywork.saff = new double[m];
        mywork.zaff = new double[m];
        mywork.W_times_dzaff = new double[m];

        mywork.best_x = new double[n];
        mywork.best_y = new double[p];
        mywork.best_z = new double[m];
        mywork.best_s = new double[m];
        mywork.best_info = new Stats();

        mywork.C = new Cone();
        mywork.C.lpc = new LPCone();
        mywork.C.lpc.p = l;
        if (l > 0) {
            mywork.C.lpc.w = new double[l];
            mywork.C.lpc.v = new double[l];
            mywork.C.lpc.kkt_idx = new int[l];
        } else {
            mywork.C.lpc.w = null;
            mywork.C.lpc.v = null;
            mywork.C.lpc.kkt_idx = null;
        }

        mywork.C.soc = (ncones == 0) ? null : new SOCone[ncones];
        mywork.C.nsoc = ncones;
        cidx = 0;
        for (i = 0; i < ncones; i++) {
            mywork.C.soc[i] = new SOCone();
            conesize = q[i];
            mywork.C.soc[i].p = conesize;
            mywork.C.soc[i].a = 0;
            mywork.C.soc[i].eta = 0;
            mywork.C.soc[i].q = new double[conesize - 1];
            mywork.C.soc[i].skbar = new double[conesize];
            mywork.C.soc[i].zkbar = new double[conesize];
            mywork.C.soc[i].Didx = new int[conesize];
            cidx += conesize;
        }

        mywork.C.nexc = nexc;
        mywork.C.expc = (nexc == 0) ? null : new ExpCone[nexc];
        for (i = 0; i < nexc; i++) { mywork.C.expc[i] = new ExpCone(); }
        mywork.C.fexv = cidx + l;

        if (cidx + l + 3 * nexc != m) {
            return null;
        }

        mywork.info = new Stats();
        if (GlblOpts.PROFILING > 1) {
            mywork.info.tfactor = 0;
            mywork.info.tkktsolve = 0;
            mywork.info.tfactor_t1 = 0;
            mywork.info.tfactor_t2 = 0;
        }

        mywork.stgs = new Settings();
        mywork.stgs.maxit = GlblOpts.MAXIT;
        mywork.stgs.gamma = GlblOpts.GAMMA;
        mywork.stgs.delta = GlblOpts.DELTA;
        mywork.stgs.eps = GlblOpts.EPS;
        mywork.stgs.nitref = GlblOpts.NITREF;
        mywork.stgs.abstol = GlblOpts.ABSTOL;
        mywork.stgs.feastol = GlblOpts.FEASTOL;
        mywork.stgs.reltol = GlblOpts.RELTOL;
        mywork.stgs.abstol_inacc = GlblOpts.ATOL_INACC;
        mywork.stgs.feastol_inacc = GlblOpts.FTOL_INACC;
        mywork.stgs.reltol_inacc = GlblOpts.RTOL_INACC;
        mywork.stgs.verbose = GlblOpts.VERBOSE;
        mywork.stgs.max_bk_iter = GlblOpts.MAX_BK;
        mywork.stgs.bk_scale = GlblOpts.BK_SCALE;
        mywork.stgs.centrality = GlblOpts.CENTRALITY;

        if (mywork.stgs.equilibrate > 0) {
            mywork.xequil = new double[n];
            mywork.Aequil = new double[p];
            mywork.Gequil = new double[m];
        }

        mywork.c = c;
        mywork.h = h;
        mywork.b = b;

        if (Apr != null && Ajc != null && Air != null) {
            mywork.A = SPLAMM.ecosCreateSparseMatrix(p, n, Ajc[n], Ajc, Air, Apr);
        } else {
            mywork.A = null;
        }
        if (Gpr != null && Gjc != null && Gir != null) {
            mywork.G = SPLAMM.ecosCreateSparseMatrix(m, n, Gjc[n], Gjc, Gir, Gpr);
        } else {
            mywork.G = SPLAMM.ecosCreateSparseMatrix(m, n, 0, Gjc, Gir, Gpr);
        }

        if (mywork.stgs.equilibrate > 0) { Equilibration.set_equilibration(mywork); }

        if (GlblOpts.PROFILING > 1) {
            mywork.info.ttranspose = 0;
            TimerOps.tic(tmattranspose);
        }
        if (mywork.A != null) {
            AtoAt = new int[mywork.A.nnz];
            At = SPLAMM.transposeSparseMatrix(mywork.A, AtoAt);
        } else {
            At = null;
            AtoAt = null;
        }
        if (GlblOpts.PROFILING > 1) mywork.info.ttranspose += TimerOps.toc(tmattranspose);

        if (GlblOpts.PROFILING > 1) TimerOps.tic(tmattranspose);
        GtoGt = new int[mywork.G.nnz];
        Gt = SPLAMM.transposeSparseMatrix(mywork.G, GtoGt);
        if (GlblOpts.PROFILING > 1) mywork.info.ttranspose += TimerOps.toc(tmattranspose);

        if (GlblOpts.PROFILING > 1) TimerOps.tic(tcreatekkt);
        if (mywork.A != null) AttoK = new int[mywork.A.nnz]; else AttoK = null;
        GttoK = new int[mywork.G.nnz];

        int[][] SignOut = new int[1][];
        SpMat[] KUOut = new SpMat[1];
        createKKT_U(Gt, At, mywork.C, SignOut, KUOut, AttoK, GttoK, mywork.stgs.staticreg);
        Sign = SignOut[0];
        KU = KUOut[0];
        if (GlblOpts.PROFILING > 1) mywork.info.tkktcreate = TimerOps.toc(tcreatekkt);

        if (mywork.A != null) {
            mywork.AtoK = new int[mywork.A.nnz];
            for (i = 0; i < mywork.A.nnz; i++) { mywork.AtoK[i] = AttoK[AtoAt[i]]; }
        } else {
            mywork.AtoK = null;
        }
        mywork.GtoK = new int[mywork.G.nnz];
        for (i = 0; i < mywork.G.nnz; i++) { mywork.GtoK[i] = GttoK[GtoGt[i]]; }

        nK = KU.n;

        mywork.KKT = new KKT();
        mywork.KKT.D = new double[nK];
        mywork.KKT.Parent = new int[nK];
        mywork.KKT.Pinv = new int[nK];
        mywork.KKT.work1 = new double[nK];
        mywork.KKT.work2 = new double[nK];
        mywork.KKT.work3 = new double[nK];
        mywork.KKT.work4 = new double[nK];
        mywork.KKT.work5 = new double[nK];
        mywork.KKT.work6 = new double[nK];
        mywork.KKT.Flag = new int[nK];
        mywork.KKT.Pattern = new int[nK];
        mywork.KKT.Lnz = new int[nK];
        mywork.KKT.RHS1 = new double[nK];
        mywork.KKT.RHS2 = new double[nK];
        mywork.KKT.dx1 = new double[mywork.n];
        mywork.KKT.dx2 = new double[mywork.n];
        mywork.KKT.dy1 = new double[mywork.p];
        mywork.KKT.dy2 = new double[mywork.p];
        mywork.KKT.dz1 = new double[mywork.m];
        mywork.KKT.dz2 = new double[mywork.m];
        mywork.KKT.Sign = new int[nK];
        mywork.KKT.PKPt = SPLAMM.newSparseMatrix(nK, nK, KU.nnz);
        mywork.KKT.PK = new int[KU.nnz];

        P = new int[nK];
        if (GlblOpts.PROFILING > 1) TimerOps.tic(tordering);
        AMD.amd_defaults(Control);
        amd_result = AMD.amd_order(nK, KU.jc, KU.ir, P, Control, Info);
        if (GlblOpts.PROFILING > 1) mywork.info.torder = TimerOps.toc(tordering);

        if (amd_result == AMD.AMD_OK) {
            if (GlblOpts.PRINTLEVEL > 2) {
                GlblOpts.printText("AMD ordering successfully computed.\n");
                AMD.amd_info(Info);
            }
        } else {
            if (GlblOpts.PRINTLEVEL > 2) {
                GlblOpts.printText("Problem in AMD ordering, exiting.\n");
                AMD.amd_info(Info);
            }
            return null;
        }

        SPLAMM.pinv(nK, P, mywork.KKT.Pinv);
        Pinv = mywork.KKT.Pinv;
        SPLAMM.permuteSparseSymmetricMatrix(KU, mywork.KKT.Pinv, mywork.KKT.PKPt, mywork.KKT.PK);

        for (i = 0; i < nK; i++) { mywork.KKT.Sign[Pinv[i]] = Sign[i]; }

        int[] Ljc = new int[nK + 1];
        LDL.ldl_symbolic2(
            mywork.KKT.PKPt.n,
            mywork.KKT.PKPt.jc,
            mywork.KKT.PKPt.ir,
            Ljc,
            mywork.KKT.Parent,
            mywork.KKT.Lnz,
            mywork.KKT.Flag
        );

        lnz = Ljc[nK];
        int[] Lir = new int[lnz];
        Lpr = new double[lnz];
        mywork.KKT.L = SPLAMM.ecosCreateSparseMatrix(nK, nK, lnz, Ljc, Lir, Lpr);

        SPLAMM.permuteSparseSymmetricMatrix(KU, mywork.KKT.Pinv, mywork.KKT.PKPt, null);

        mywork.rx = (n == 0) ? null : new double[n];
        mywork.ry = (p == 0) ? null : new double[p];
        mywork.rz = (m == 0) ? null : new double[m];

        mywork.KKT.P = P;
        Sign = null;
        if (At != null) {
            SPLAMM.freeSparseMatrix(At);
            AtoAt = null;
            AttoK = null;
        }
        SPLAMM.freeSparseMatrix(Gt);
        SPLAMM.freeSparseMatrix(KU);
        GtoGt = null;
        GttoK = null;

        if (GlblOpts.PROFILING > 0) mywork.info.tsetup = TimerOps.toc(tsetup);

        return mywork;
    }

    static void createKKT_U(SpMat Gt, SpMat At, Cone C, int[][] SOut, SpMat[] KOut,
                              int[] AttoK, int[] GttoK, int staticreg) {
        int i, j, k, l, r, row_stop, row, cone_strt, ks, conesize, exp_cone_strt;
        int n = Gt.m;
        int m = Gt.n;
        int p = (At != null) ? At.n : 0;
        int nK, nnzK;
        double[] Kpr;
        int[] Kjc, Kir, Sign;

        nK = n + p + m + 2 * C.nsoc;

        nnzK = ((At != null) ? At.nnz : 0) + Gt.nnz + C.lpc.p;
        if (staticreg > 0) {
            nnzK += n + p;
        }
        for (i = 0; i < C.nsoc; i++) {
            nnzK += 3 * C.soc[i].p + 1;
        }
        nnzK += 6 * C.nexc;

        Kpr = new double[nnzK];
        Kir = new int[nnzK];
        Kjc = new int[nK + 1];
        Sign = new int[nK];

        for (ks = 0; ks < n; ks++) { Sign[ks] = +1; }
        for (ks = n; ks < n + p; ks++) { Sign[ks] = -1; }
        for (ks = n + p; ks < n + p + C.lpc.p; ks++) { Sign[ks] = -1; }
        ks = n + p + C.lpc.p;
        for (l = 0; l < C.nsoc; l++) {
            for (i = 0; i < C.soc[l].p; i++) { Sign[ks++] = -1; }
            Sign[ks++] = -1;
            Sign[ks++] = +1;
        }
        for (ks = nK - 3 * C.nexc; ks < nK; ) { Sign[ks++] = -1; }

        k = 0;
        if (staticreg > 0) {
            for (j = 0; j < n; j++) {
                Kjc[j] = j;
                Kir[j] = j;
                Kpr[k++] = GlblOpts.DELTASTAT;
            }
        } else {
            for (j = 0; j < n; j++) {
                Kjc[j] = 0;
            }
        }

        i = 0;
        if (At != null) {
            for (j = 0; j < p; j++) {
                row = At.jc[j];
                row_stop = At.jc[j + 1];
                if (row <= row_stop) {
                    Kjc[n + j] = k;
                    while (row++ < row_stop) {
                        Kir[k] = At.ir[i];
                        Kpr[k] = At.pr[i];
                        AttoK[i++] = k++;
                    }
                }
                if (staticreg > 0) {
                    Kir[k] = n + j;
                    Kpr[k++] = -GlblOpts.DELTASTAT;
                }
            }
        }

        i = 0;
        for (j = 0; j < C.lpc.p; j++) {
            row = Gt.jc[j];
            row_stop = Gt.jc[j + 1];
            if (row <= row_stop) {
                Kjc[n + p + j] = k;
                while (row++ < row_stop) {
                    Kir[k] = Gt.ir[i];
                    Kpr[k] = Gt.pr[i];
                    GttoK[i++] = k++;
                }
            }
            C.lpc.kkt_idx[j] = k;
            Kir[k] = n + p + j;
            Kpr[k] = -1.0;
            k++;
        }

        cone_strt = C.lpc.p;
        for (l = 0; l < C.nsoc; l++) {
            conesize = C.soc[l].p;
            for (j = 0; j < conesize; j++) {
                row = Gt.jc[cone_strt + j];
                row_stop = Gt.jc[cone_strt + j + 1];
                if (row <= row_stop) {
                    Kjc[n + p + cone_strt + 2 * l + j] = k;
                    while (row++ < row_stop) {
                        Kir[k] = Gt.ir[i];
                        Kpr[k] = Gt.pr[i];
                        GttoK[i++] = k++;
                    }
                }
                Kir[k] = n + p + cone_strt + 2 * l + j;
                Kpr[k] = -1.0;
                C.soc[l].Didx[j] = k;
                k++;
            }

            Kjc[n + p + cone_strt + 2 * l + conesize] = k;
            for (r = 1; r < conesize; r++) {
                Kir[k] = n + p + cone_strt + 2 * l + r;
                Kpr[k] = 0;
                k++;
            }
            Kir[k] = n + p + cone_strt + 2 * l + conesize;
            Kpr[k] = -1;
            k++;

            Kjc[n + p + cone_strt + 2 * l + conesize + 1] = k;
            for (r = 0; r < conesize; r++) {
                Kir[k] = n + p + cone_strt + 2 * l + r;
                Kpr[k] = 0;
                k++;
            }
            Kir[k] = n + p + cone_strt + 2 * l + conesize + 1;
            Kpr[k] = +1;
            k++;

            cone_strt += C.soc[l].p;
        }

        exp_cone_strt = cone_strt + 2 * C.nsoc;
        for (l = 0; l < C.nexc; l++) {
            for (j = 0; j < 3; j++) {
                row = Gt.jc[cone_strt + j];
                row_stop = Gt.jc[cone_strt + j + 1];
                if (row <= row_stop) {
                    Kjc[n + p + exp_cone_strt + j] = k;
                    while (row++ < row_stop) {
                        Kir[k] = Gt.ir[i];
                        Kpr[k] = Gt.pr[i];
                        GttoK[i++] = k++;
                    }
                }
                C.expc[l].colstart[j] = k;
                for (r = 0; r < j; r++) {
                    Kir[k] = n + p + exp_cone_strt + r;
                    Kpr[k] = 0.0;
                    k++;
                }
                Kir[k] = n + p + exp_cone_strt + j;
                Kpr[k++] = -1.0;
            }
            cone_strt += 3;
            exp_cone_strt += 3;
        }

        for (int ii = k; ii < nnzK; ii++) { Kjc[nK] = nnzK; }

        SOut[0] = Sign;
        KOut[0] = SPLAMM.ecosCreateSparseMatrix(nK, nK, nnzK, Kjc, Kir, Kpr);
    }

    public static void ECOS_cleanup(PWork w, int keepvars) {
        int i, eqflag;
        eqflag = (w.stgs != null) ? w.stgs.equilibrate : 0;
        if (eqflag > 0) { Equilibration.unset_equilibration(w); }

        w.KKT.D = null;
        w.KKT.dx1 = null;
        w.KKT.dx2 = null;
        w.KKT.dy1 = null;
        w.KKT.dy2 = null;
        w.KKT.dz1 = null;
        w.KKT.dz2 = null;
        w.KKT.Flag = null;
        SPLAMM.freeSparseMatrix(w.KKT.L);
        w.KKT.Lnz = null;
        w.KKT.Parent = null;
        w.KKT.Pattern = null;
        w.KKT.Sign = null;
        w.KKT.Pinv = null;
        w.KKT.P = null;
        w.KKT.PK = null;
        SPLAMM.freeSparseMatrix(w.KKT.PKPt);
        w.KKT.RHS1 = null;
        w.KKT.RHS2 = null;
        w.KKT.work1 = null;
        w.KKT.work2 = null;
        w.KKT.work3 = null;
        w.KKT.work4 = null;
        w.KKT.work5 = null;
        w.KKT.work6 = null;
        w.KKT = null;
        if (w.A != null) { w.AtoK = null; }
        w.GtoK = null;

        if (w.C.lpc.p > 0) {
            w.C.lpc.kkt_idx = null;
            w.C.lpc.v = null;
            w.C.lpc.w = null;
        }
        w.C.lpc = null;

        for (i = 0; i < w.C.nsoc; i++) {
            w.C.soc[i].q = null;
            w.C.soc[i].skbar = null;
            w.C.soc[i].zkbar = null;
            w.C.soc[i].Didx = null;
        }
        if (w.C.nsoc > 0) { w.C.soc = null; }
        if (w.C.nexc > 0) { w.C.expc = null; }
        w.C = null;

        w.W_times_dzaff = null;
        w.dsaff_by_W = null;
        w.dzaff = null;
        w.dsaff = null;
        w.zaff = null;
        w.saff = null;
        w.info = null;
        w.best_info = null;
        w.lambda = null;
        w.rx = null;
        w.ry = null;
        w.rz = null;
        w.stgs = null;
        w.G = null;
        if (w.A != null) w.A = null;
        w.best_z = null;
        w.best_s = null;
        w.best_y = null;
        w.best_x = null;
        if (keepvars < 4) { w.z = null; }
        if (keepvars < 3) { w.s = null; }
        if (keepvars < 2) { w.y = null; }
        if (keepvars < 1) { w.x = null; }
        if (eqflag > 0) {
            w.xequil = null;
            w.Aequil = null;
            w.Gequil = null;
        }
    }

    public static String ECOS_ver() {
        return GlblOpts.ECOS_VERSION;
    }

    public static void ecos_updateDataEntry_h(PWork w, int idx, double value) {
        if (w.stgs.equilibrate > 0) {
            w.h[idx] = value / w.Gequil[idx];
        } else {
            w.h[idx] = value;
        }
    }

    public static void ecos_updateDataEntry_c(PWork w, int idx, double value) {
        w.c[idx] = value;
    }

    public static void ECOS_updateData(PWork w, double[] Gpr, double[] Apr,
                                        double[] c, double[] h, double[] b) {
        int k;

        if (w.stgs.equilibrate > 0) { Equilibration.unset_equilibration(w); }

        if (w.G != null) {
            w.G.pr = Gpr;
            w.h = h;
        }
        if (w.A != null) {
            w.A.pr = Apr;
            w.b = b;
        }
        w.c = c;

        if (w.stgs.equilibrate > 0) { Equilibration.set_equilibration(w); }

        if (w.A != null) {
            for (k = 0; k < w.A.nnz; k++) {
                w.KKT.PKPt.pr[w.KKT.PK[w.AtoK[k]]] = Apr[k];
            }
        }
        if (w.G != null) {
            for (k = 0; k < w.G.nnz; k++) {
                w.KKT.PKPt.pr[w.KKT.PK[w.GtoK[k]]] = Gpr[k];
            }
        }
    }
}
