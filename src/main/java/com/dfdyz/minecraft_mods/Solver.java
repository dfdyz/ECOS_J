package com.dfdyz.minecraft_mods;

import com.dfdyz.ecos_j.ECOS;
import com.dfdyz.ecos_j.PWork;

public class Solver implements AutoCloseable {
    private final PWork work;
    private final double[] gpr, apr, c, h, b;
    private final int[] gprBuf, aprBuf;

    public Solver(int n, int m, int p, int l, int ncones, int[] q, int nexc,
                  double[] Gpr, int[] Gjc, int[] Gir,
                  double[] Apr, int[] Ajc, int[] Air,
                  double[] c, double[] h, double[] b) {
        this.gpr = Gpr;
        this.apr = Apr;
        this.c = c;
        this.h = h;
        this.b = b;
        this.work = ECOS.ECOS_setup(n, m, p, l, ncones, q, nexc,
                                    Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
        this.gprBuf = null;
        this.aprBuf = null;
    }

    public double[] getGprBuffer() { return gpr; }
    public double[] getAprBuffer() { return apr; }
    public double[] getHBuffer() { return h; }
    public double[] getBBuffer() { return b; }
    public double[] getCBuffer() { return c; }
    public double[] getX() { return work.x; }

    public int solve() {
        return ECOS.ECOS_solve(work);
    }

    public void updateData(double[] gprNew, double[] aprNew, double[] cNew, double[] hNew, double[] bNew) {
        ECOS.ECOS_updateData(work, gprNew, aprNew, cNew, hNew, bNew);
    }

    public void updateDataInline() {
        if (work.A != null) {
            for (int k = 0; k < work.A.nnz; k++) {
                work.KKT.PKPt.pr[work.KKT.PK[work.AtoK[k]]] = apr[k];
            }
        }
        if (work.G != null) {
            for (int k = 0; k < work.G.nnz; k++) {
                work.KKT.PKPt.pr[work.KKT.PK[work.GtoK[k]]] = gpr[k];
            }
        }
    }

    @Override
    public void close() {
        ECOS.ECOS_cleanup(work, 0);
    }
}
