package com.dfdyz.ecos_j;

public class PWork {
    public int n;
    public int m;
    public int p;
    public int D;

    public double[] x;
    public double[] y;
    public double[] z;
    public double[] s;
    public double[] lambda;
    public double kap;
    public double tau;

    public double[] best_x;
    public double[] best_y;
    public double[] best_z;
    public double[] best_s;
    public double best_kap;
    public double best_tau;
    public double best_cx;
    public double best_by;
    public double best_hz;
    public Stats best_info;

    public double[] dsaff;
    public double[] dzaff;
    public double[] W_times_dzaff;
    public double[] dsaff_by_W;
    public double[] saff;
    public double[] zaff;

    public Cone C;

    public SpMat A;
    public SpMat G;
    public double[] c;
    public double[] b;
    public double[] h;

    public int[] AtoK;
    public int[] GtoK;

    public double[] xequil;
    public double[] Aequil;
    public double[] Gequil;

    public double resx0;
    public double resy0;
    public double resz0;

    public double[] rx;
    public double[] ry;
    public double[] rz;
    public double rt;
    public double hresx;
    public double hresy;
    public double hresz;

    public double nx, ny, nz, ns;

    public double cx;
    public double by;
    public double hz;
    public double sz;

    public KKT KKT;

    public Stats info;
    public Settings stgs;
}
