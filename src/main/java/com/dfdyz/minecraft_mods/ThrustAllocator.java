package com.dfdyz.minecraft_mods;

import org.joml.Vector3d;
import java.util.List;

public class ThrustAllocator implements AutoCloseable {
    private final Ship ship;
    private final int N, R, F, nVars, mCones, lDim;
    private final Solver solver;
    private final int[] Acols;
    private final boolean[] isRot;
    private final int[] rotOf, fixOf;
    private final Vector3d[] dirs;
    private final Vector3d[] r;

    public ThrustAllocator(Ship ship) {
        this.ship = ship;
        List<Thrust> thrusters = ship.getThrusts();
        N = thrusters.size();
        int rCnt = 0, fCnt = 0;
        isRot = new boolean[N];
        rotOf = new int[N];
        fixOf = new int[N];
        dirs = new Vector3d[N];
        for (int i = 0; i < N; i++) {
            Thrust t = thrusters.get(i);
            dirs[i] = t.getDir();
            isRot[i] = t.getLimitType() == LimitType.HALF_BALL;
            if (isRot[i]) { rotOf[i] = rCnt; fixOf[i] = -1; rCnt++; }
            else { fixOf[i] = fCnt; rotOf[i] = -1; fCnt++; }
        }
        R = rCnt; F = fCnt;
        nVars = 4 * R + F;
        mCones = 6 * R + 2 * F;
        lDim = 2 * R + 2 * F;
        int nSoc = R;

        r = new Vector3d[N];
        for (int i = 0; i < N; i++) r[i] = new Vector3d();
        Vector3d[] initR = leverArms();

        int[] Gcols = new int[nVars + 1];
        for (int i = 0; i < N; i++) {
            if (isRot[i]) {
                int ri = rotOf[i];
                Gcols[3*ri]=2; Gcols[3*ri+1]=2; Gcols[3*ri+2]=2;
                Gcols[3*R+ri]=2;
            } else {
                Gcols[4*R+fixOf[i]]=2;
            }
        }
        int tg = 0;
        for (int j = 0; j < nVars; j++) { int c = Gcols[j]; Gcols[j] = tg; tg += c; }
        Gcols[nVars] = tg;

        int[] Gir = new int[tg];
        double[] Gpr = new double[tg];
        int[] Gp = new int[nVars];
        for (int j = 0; j < nVars; j++) Gp[j] = Gcols[j];

        for (int i = 0; i < N; i++) {
            if (!isRot[i]) continue;
            int rr = rotOf[i];
            int vx = 3*rr, vy = 3*rr+1, vz = 3*rr+2, tt = 3*R+rr;
            int lp = rr, sb = 2*R+2*F+4*rr;
            int p;
            p = Gp[vx]; Gir[p] = lp; Gpr[p] = -dirs[i].x; Gir[p+1] = sb+1; Gpr[p+1] = -1.0; Gp[vx] = p+2;
            p = Gp[vy]; Gir[p] = lp; Gpr[p] = -dirs[i].y; Gir[p+1] = sb+2; Gpr[p+1] = -1.0; Gp[vy] = p+2;
            p = Gp[vz]; Gir[p] = lp; Gpr[p] = -dirs[i].z; Gir[p+1] = sb+3; Gpr[p+1] = -1.0; Gp[vz] = p+2;
            p = Gp[tt]; Gir[p] = R+rr; Gpr[p] = 1.0; Gir[p+1] = sb; Gpr[p+1] = -1.0; Gp[tt] = p+2;
        }
        for (int i = 0; i < N; i++) {
            if (isRot[i]) continue;
            int fj = fixOf[i];
            int ca = 4*R+fj;
            int p = Gp[ca];
            Gir[p] = 2*R+fj; Gpr[p] = -1.0;
            Gir[p+1] = 2*R+F+fj; Gpr[p+1] = 1.0;
            Gp[ca] = p+2;
        }

        Acols = new int[nVars+1];
        for (int i = 0; i < N; i++) {
            if (isRot[i]) {
                int ri = rotOf[i];
                Acols[3*ri]=3; Acols[3*ri+1]=3; Acols[3*ri+2]=3;
                Acols[3*R+ri]=0;
            } else {
                Acols[4*R+fixOf[i]]=6;
            }
        }
        int ta = 0;
        for (int j = 0; j < nVars; j++) { int c = Acols[j]; Acols[j] = ta; ta += c; }
        Acols[nVars] = ta;

        int[] Air = new int[ta];
        double[] Apr = new double[ta];
        int[] Ap = new int[nVars];
        System.arraycopy(Acols, 0, Ap, 0, nVars);
        fillApr(Apr, Air, Ap, initR);

        double[] c = new double[nVars];
        for (int ri = 0; ri < R; ri++) c[3*R+ri] = 1.0;
        for (int fj = 0; fj < F; fj++) c[4*R+fj] = 1.0;

        double[] h = new double[mCones];
        for (int i = 0; i < N; i++) {
            double fm = ship.getThrusts().get(i).getFmax();
            if (isRot[i]) h[R+rotOf[i]] = fm;
            else h[2*R+F+fixOf[i]] = fm;
        }

        int[] q = new int[nSoc];
        for (int ri = 0; ri < R; ri++) q[ri] = 4;

        double[] bInit = new double[]{0,0,0,0,0,0};

        solver = new Solver(nVars, mCones, 6, lDim, nSoc, q, 0,
            Gpr, Gcols, Gir, Apr, Acols, Air, c, h, bInit);
    }

    public int lastExitFlag;

    public void solve() { solve(ship.getTargetF(), ship.getTargetT()); }

    public void solve(Vector3d F_target, Vector3d T_target) {
        double[] gpr = solver.getGprBuffer();
        for (int i = 0; i < N; i++) {
            if (isRot[i]) {
                int ri = rotOf[i];
                gpr[6*ri] = -dirs[i].x;
                gpr[6*ri+2] = -dirs[i].y;
                gpr[6*ri+4] = -dirs[i].z;
            }
        }
        leverArms();
        updateAprRaw();

        double[] hBuf = solver.getHBuffer();
        for (int i = 0; i < N; i++) {
            double fm = ship.getThrusts().get(i).getFmax();
            if (isRot[i]) hBuf[R+rotOf[i]] = fm;
            else hBuf[2*R+F+fixOf[i]] = fm;
        }

        double[] b = solver.getBBuffer();
        b[0] = F_target.x; b[1] = F_target.y; b[2] = F_target.z;
        b[3] = T_target.x; b[4] = T_target.y; b[5] = T_target.z;

        /*solver.updateData(solver.getGprBuffer(), solver.getAprBuffer(),
            solver.getCBuffer(), solver.getHBuffer(), solver.getBBuffer());*/
        solver.updateDataInline();

        lastExitFlag = solver.solve();

        double[] x = solver.getX();
        for (int i = 0; i < N; i++) {
            Vector3d outF = ship.getThrusts().get(i).getForce();
            if (isRot[i]) {
                int ri = rotOf[i];
                outF.set(x[3*ri], x[3*ri+1], x[3*ri+2]);
            } else {
                double a = x[4*R+fixOf[i]];
                outF.set(a*dirs[i].x, a*dirs[i].y, a*dirs[i].z);
            }
        }
    }

    @Override public void close() { }

    private Vector3d[] leverArms() {
        List<Thrust> ts = ship.getThrusts();
        var c = ship.getCenter();
        for (int i = 0; i < N; i++) r[i].set(ts.get(i).getPos()).sub(c);
        return r;
    }

    private void updateAprRaw() {
        double[] b = solver.getAprBuffer();
        if (b == null) return;
        for (int i = 0; i < N; i++) {
            double rx = r[i].x, ry = r[i].y, rz = r[i].z;
            if (isRot[i]) {
                int ri = rotOf[i];
                int p;
                p = Acols[3*ri]; b[p]=1.0; b[p+1]=rz; b[p+2]=-ry;
                p = Acols[3*ri+1]; b[p]=1.0; b[p+1]=-rz; b[p+2]=rx;
                p = Acols[3*ri+2]; b[p]=1.0; b[p+1]=ry; b[p+2]=-rx;
            } else {
                int fj = fixOf[i];
                int pa = Acols[4*R+fj];
                double dx=dirs[i].x, dy=dirs[i].y, dz=dirs[i].z;
                b[pa]=dx; b[pa+1]=dy; b[pa+2]=dz;
                b[pa+3]=ry*dz-rz*dy; b[pa+4]=rz*dx-rx*dz; b[pa+5]=rx*dy-ry*dx;
            }
        }
    }

    private void fillApr(double[] Apr, int[] Air, int[] Ap, Vector3d[] r) {
        for (int i = 0; i < N; i++) {
            double rx=r[i].x, ry=r[i].y, rz=r[i].z;
            if (isRot[i]) {
                int ri = rotOf[i];
                int p;
                p = Ap[3*ri]; Air[p]=0; Apr[p]=1.0; Air[p+1]=4; Apr[p+1]=rz; Air[p+2]=5; Apr[p+2]=-ry; Ap[3*ri]=p+3;
                p = Ap[3*ri+1]; Air[p]=1; Apr[p]=1.0; Air[p+1]=3; Apr[p+1]=-rz; Air[p+2]=5; Apr[p+2]=rx; Ap[3*ri+1]=p+3;
                p = Ap[3*ri+2]; Air[p]=2; Apr[p]=1.0; Air[p+1]=3; Apr[p+1]=ry; Air[p+2]=4; Apr[p+2]=-rx; Ap[3*ri+2]=p+3;
            } else {
                int fj = fixOf[i];
                int p = Ap[4*R+fj];
                double dx=dirs[i].x, dy=dirs[i].y, dz=dirs[i].z;
                Air[p]=0; Apr[p]=dx; Air[p+1]=1; Apr[p+1]=dy; Air[p+2]=2; Apr[p+2]=dz;
                Air[p+3]=3; Apr[p+3]=ry*dz-rz*dy;
                Air[p+4]=4; Apr[p+4]=rz*dx-rx*dz;
                Air[p+5]=5; Apr[p+5]=rx*dy-ry*dx;
                Ap[4*R+fj]=p+6;
            }
        }
    }
}
