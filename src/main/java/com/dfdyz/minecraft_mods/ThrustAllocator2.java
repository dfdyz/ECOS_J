package com.dfdyz.minecraft_mods;


import org.joml.Vector3d;
import java.util.List;

import static com.dfdyz.minecraft_mods.LimitType.ANGLE_LIMIT;
import static com.dfdyz.minecraft_mods.LimitType.*;

public class ThrustAllocator2 implements AutoCloseable {
    private final int S, R, A, F, nVars, mCones, lDim;
    private final Ship ship;
    private final Solver solver;
    private final int[] Acols;
    private final LimitType[] ltype;
    private final int[] socIdx, fixIdx;
    private final Vector3d[] dirs;
    private final int[] tColOffset;
    private final Vector3d[] r;

    public ThrustAllocator2(Ship ship) {
        this.ship = ship;
        List<Thrust> ts = ship.getThrusts();
        int N = ts.size();
        ltype = new LimitType[N];
        socIdx = new int[N];
        fixIdx = new int[N];
        dirs = new Vector3d[N];
        int rCnt=0, aCnt=0, fCnt=0;
        for (int i = 0; i < N; i++) {
            Thrust t = ts.get(i);
            ltype[i] = t.getLimitType();
            dirs[i] = t.getDir();
            switch (ltype[i]) {
                case HALF_BALL: socIdx[i] = rCnt++; fixIdx[i] = -1; break;
                case ANGLE_LIMIT: socIdx[i] = rCnt + aCnt; aCnt++; fixIdx[i] = -1; break;
                case STATIC: socIdx[i] = -1; fixIdx[i] = fCnt++; break;
            }
        }
        R = rCnt; A = aCnt; F = fCnt; S = R + A;
        nVars = 4*S + F;
        mCones = 6*S + 2*F;
        lDim = 2*S + 2*F;
        int nSoc = S;

        r = new Vector3d[N];
        for (int i = 0; i < N; i++) r[i] = new Vector3d();
        Vector3d[] initR = leverArms();

        tColOffset = new int[S];
        int tcp = 6*S;
        for (int gi = 0; gi < S; gi++) { tColOffset[gi] = tcp; tcp += 2; }
        for (int i = 0; i < N; i++) {
            if (ltype[i] == ANGLE_LIMIT) {
                int gi = socIdx[i];
                for (int g = gi+1; g < S; g++) tColOffset[g]++;
            }
        }

        int[] Gcols = new int[nVars+1];
        for (int i = 0; i < N; i++) {
            if (ltype[i] == STATIC) {
                Gcols[4*S + fixIdx[i]] = 2;
            } else {
                int gi = socIdx[i];
                Gcols[3*gi] = 2; Gcols[3*gi+1] = 2; Gcols[3*gi+2] = 2;
                Gcols[3*S + gi] = (ltype[i] == ANGLE_LIMIT) ? 3 : 2;
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
            if (ltype[i] != HALF_BALL) continue;
            int gi = socIdx[i];
            int vx=3*gi, vy=3*gi+1, vz=3*gi+2, tt=3*S+gi;
            int hemiRow=gi, maxRow=S+gi, socBase=2*S+2*F+4*gi;
            int p;
            p=Gp[vx]; Gir[p]=hemiRow; Gpr[p]=-dirs[i].x; Gir[p+1]=socBase+1; Gpr[p+1]=-1.0; Gp[vx]=p+2;
            p=Gp[vy]; Gir[p]=hemiRow; Gpr[p]=-dirs[i].y; Gir[p+1]=socBase+2; Gpr[p+1]=-1.0; Gp[vy]=p+2;
            p=Gp[vz]; Gir[p]=hemiRow; Gpr[p]=-dirs[i].z; Gir[p+1]=socBase+3; Gpr[p+1]=-1.0; Gp[vz]=p+2;
            p=Gp[tt]; Gir[p]=maxRow; Gpr[p]=1.0; Gir[p+1]=socBase; Gpr[p+1]=-1.0; Gp[tt]=p+2;
        }
        for (int i = 0; i < N; i++) {
            if (ltype[i] != ANGLE_LIMIT) continue;
            int gi = socIdx[i];
            int vx=3*gi, vy=3*gi+1, vz=3*gi+2, tt=3*S+gi;
            int angleRow=gi, maxRow=S+gi, socBase=2*S+2*F+4*gi;
            double ct = ts.get(i).getCosTheta();
            int p;
            p=Gp[vx]; Gir[p]=angleRow; Gpr[p]=-dirs[i].x; Gir[p+1]=socBase+1; Gpr[p+1]=-1.0; Gp[vx]=p+2;
            p=Gp[vy]; Gir[p]=angleRow; Gpr[p]=-dirs[i].y; Gir[p+1]=socBase+2; Gpr[p+1]=-1.0; Gp[vy]=p+2;
            p=Gp[vz]; Gir[p]=angleRow; Gpr[p]=-dirs[i].z; Gir[p+1]=socBase+3; Gpr[p+1]=-1.0; Gp[vz]=p+2;
            p=Gp[tt]; Gir[p]=maxRow; Gpr[p]=1.0; Gir[p+1]=socBase; Gpr[p+1]=-1.0;
                        Gir[p+2]=angleRow; Gpr[p+2]=ct; Gp[tt]=p+3;
        }
        for (int i = 0; i < N; i++) {
            if (ltype[i] != STATIC) continue;
            int fj = fixIdx[i];
            int ca = 4*S + fj;
            int p = Gp[ca];
            Gir[p] = 2*S+fj; Gpr[p] = -1.0;
            Gir[p+1] = 2*S+F+fj; Gpr[p+1] = 1.0;
            Gp[ca] = p+2;
        }

        Acols = new int[nVars+1];
        for (int i = 0; i < N; i++) {
            if (ltype[i] == STATIC) {
                Acols[4*S + fixIdx[i]] = 6;
            } else {
                int gi = socIdx[i];
                Acols[3*gi]=3; Acols[3*gi+1]=3; Acols[3*gi+2]=3;
                Acols[3*S+gi]=0;
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
        for (int gi = 0; gi < S; gi++) c[3*S+gi] = 1.0;
        for (int fj = 0; fj < F; fj++) c[4*S+fj] = 1.0;

        double[] h = new double[mCones];
        for (int i = 0; i < N; i++) {
            double fm = ship.getThrusts().get(i).getFmax();
            if (ltype[i] != STATIC) h[S + socIdx[i]] = fm;
            else h[2*S + F + fixIdx[i]] = fm;
        }

        int[] q = new int[nSoc];
        for (int gi = 0; gi < S; gi++) q[gi] = 4;

        double[] bInit = new double[]{0,0,0,0,0,0};

        solver = new Solver(nVars, mCones, 6, lDim, nSoc, q, 0,
            Gpr, Gcols, Gir, Apr, Acols, Air, c, h, bInit);
    }

    public int lastExitFlag;

    public void solve() { solve(ship.getTargetF(), ship.getTargetT()); }

    public void solve(Vector3d F_target, Vector3d T_target) {
        double[] gpr = solver.getGprBuffer();
        for (int i = 0; i < ship.getThrusts().size(); i++) {
            if (ltype[i] == STATIC) continue;
            int gi = socIdx[i];
            gpr[6*gi] = -dirs[i].x;
            gpr[6*gi+2] = -dirs[i].y;
            gpr[6*gi+4] = -dirs[i].z;
            if (ltype[i] == ANGLE_LIMIT) {
                double ct = ship.getThrusts().get(i).getCosTheta();
                gpr[tColOffset[gi] + 2] = ct;
            }
        }
        leverArms();
        updateAprRaw();

        double[] hBuf = solver.getHBuffer();
        List<Thrust> ts = ship.getThrusts();
        for (int i = 0; i < ts.size(); i++) {
            double fm = ts.get(i).getFmax();
            if (ltype[i] != STATIC) hBuf[S + socIdx[i]] = fm;
            else hBuf[2*S + F + fixIdx[i]] = fm;
        }

        double[] b = solver.getBBuffer();
        b[0]=F_target.x; b[1]=F_target.y; b[2]=F_target.z;
        b[3]=T_target.x; b[4]=T_target.y; b[5]=T_target.z;

        solver.updateData(solver.getGprBuffer(), solver.getAprBuffer(),
            solver.getCBuffer(), solver.getHBuffer(), solver.getBBuffer());

        lastExitFlag = solver.solve();

        double[] x = solver.getX();
        for (int i = 0; i < ts.size(); i++) {
            Vector3d f = ts.get(i).getForce();
            if (ltype[i] != STATIC) {
                int gi = socIdx[i];
                f.set(x[3*gi], x[3*gi+1], x[3*gi+2]);
            } else {
                double a = x[4*S + fixIdx[i]];
                f.set(a*dirs[i].x, a*dirs[i].y, a*dirs[i].z);
            }
        }
    }

    @Override public void close() {  }

    private Vector3d[] leverArms() {
        List<Thrust> ts = ship.getThrusts();
        var c = ship.getCenter();
        for (int i = 0; i < ts.size(); i++) r[i].set(ts.get(i).getPos()).sub(c);
        return r;
    }

    private void updateAprRaw() {
        double[] b = solver.getAprBuffer();
        if (b == null) return;
        List<Thrust> ts = ship.getThrusts();
        for (int i = 0; i < ts.size(); i++) {
            double rx=r[i].x, ry=r[i].y, rz=r[i].z;
            if (ltype[i] != STATIC) {
                int gi = socIdx[i];
                int p;
                p = Acols[3*gi]; b[p]=1.0; b[p+1]=rz; b[p+2]=-ry;
                p = Acols[3*gi+1]; b[p]=1.0; b[p+1]=-rz; b[p+2]=rx;
                p = Acols[3*gi+2]; b[p]=1.0; b[p+1]=ry; b[p+2]=-rx;
            } else {
                int fj = fixIdx[i];
                int pa = Acols[4*S+fj];
                double dx=dirs[i].x, dy=dirs[i].y, dz=dirs[i].z;
                b[pa]=dx; b[pa+1]=dy; b[pa+2]=dz;
                b[pa+3]=ry*dz-rz*dy; b[pa+4]=rz*dx-rx*dz; b[pa+5]=rx*dy-ry*dx;
            }
        }
    }

    private void fillApr(double[] Apr, int[] Air, int[] Ap, Vector3d[] r) {
        for (int i = 0; i < ship.getThrusts().size(); i++) {
            double rx=r[i].x, ry=r[i].y, rz=r[i].z;
            if (ltype[i] != STATIC) {
                int gi = socIdx[i];
                int p;
                p = Ap[3*gi]; Air[p]=0; Apr[p]=1.0; Air[p+1]=4; Apr[p+1]=rz; Air[p+2]=5; Apr[p+2]=-ry; Ap[3*gi]=p+3;
                p = Ap[3*gi+1]; Air[p]=1; Apr[p]=1.0; Air[p+1]=3; Apr[p+1]=-rz; Air[p+2]=5; Apr[p+2]=rx; Ap[3*gi+1]=p+3;
                p = Ap[3*gi+2]; Air[p]=2; Apr[p]=1.0; Air[p+1]=3; Apr[p+1]=ry; Air[p+2]=4; Apr[p+2]=-rx; Ap[3*gi+2]=p+3;
            } else {
                int fj = fixIdx[i];
                int p = Ap[4*S+fj];
                double dx=dirs[i].x, dy=dirs[i].y, dz=dirs[i].z;
                Air[p]=0; Apr[p]=dx; Air[p+1]=1; Apr[p+1]=dy; Air[p+2]=2; Apr[p+2]=dz;
                Air[p+3]=3; Apr[p+3]=ry*dz-rz*dy;
                Air[p+4]=4; Apr[p+4]=rz*dx-rx*dz;
                Air[p+5]=5; Apr[p+5]=rx*dy-ry*dx;
                Ap[4*S+fj] = p+6;
            }
        }
    }
}
