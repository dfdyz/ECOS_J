package com.dfdyz.ecos_j;

public final class GlblOpts {
    private GlblOpts() {}

    public static final int PRINTLEVEL = 0;
    public static final int PROFILING = 1;
    public static final int DEBUG = 0;

    public static final double ECOS_INFINITY = Double.MAX_VALUE + Double.MAX_VALUE;
    public static final double ECOS_NAN = ECOS_INFINITY - ECOS_INFINITY;

    public static final String ECOS_VERSION = "2.0.10";

    public static final int MAXIT = 150;
    public static final double FEASTOL = 1E-8;
    public static final double ABSTOL = 1E-8;
    public static final double RELTOL = 1E-8;
    public static final double FTOL_INACC = 1E-4;
    public static final double ATOL_INACC = 5E-5;
    public static final double RTOL_INACC = 5E-5;
    public static final double GAMMA = 0.99;
    public static final int STATICREG = 1;
    public static final double DELTASTAT = 7E-8;
    public static final double DELTA = 2E-7;
    public static final double EPS = 1E-13;
    public static final int VERBOSE = 1;
    public static final int NITREF = 9;
    public static final int IRERRFACT = 6;
    public static final double LINSYSACC = 1E-14;
    public static final double SIGMAMIN = 1E-4;
    public static final double SIGMAMAX = 1.0;
    public static final double STEPMIN = 1E-6;
    public static final double STEPMAX = 0.999;
    public static final double SAFEGUARD = 500;

    public static final int MAX_BK = 90;
    public static final double BK_SCALE = 0.8;
    public static final double MIN_DISTANCE = 0.1;
    public static final double CENTRALITY = 1;

    public static final int EQUILIBRATE = 0;
    public static final int EQUIL_ITERS = 3;

    public static final int CONEMODE = 0;

    public static final int ECOS_OPTIMAL = 0;
    public static final int ECOS_PINF = 1;
    public static final int ECOS_DINF = 2;
    public static final int ECOS_INACC_OFFSET = 10;
    public static final int ECOS_MAXIT = -1;
    public static final int ECOS_NUMERICS = -2;
    public static final int ECOS_OUTCONE = -3;
    public static final int ECOS_SIGINT = -4;
    public static final int ECOS_FATAL = -7;
    public static final int ECOS_NOT_CONVERGED_YET = -87;

    public static final int INSIDE_CONE = 0;
    public static final int OUTSIDE_CONE = 1;

    public static final int KKT_PROBLEM = 0;
    public static final int KKT_OK = 1;

    public static double max(double a, double b) { return a < b ? b : a; }
    public static int max(int a, int b) { return a < b ? b : a; }

    public static double safeDivPos(double x, double y) {
        return (y < EPS ? (x / EPS) : (x / y));
    }

    public static double fabs(double x) { return x < 0 ? -x : x; }

    public static void printText(String format, Object... args) {
        System.out.printf(format, args);
    }

    public static boolean isInf(double x) { return Double.isInfinite(x); }
    public static boolean isNan(double x) { return Double.isNaN(x); }
}
