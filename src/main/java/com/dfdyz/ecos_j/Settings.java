package com.dfdyz.ecos_j;

public class Settings {
    public double gamma;
    public double delta;
    public double eps;
    public double feastol;
    public double abstol;
    public double reltol;
    public double feastol_inacc;
    public double abstol_inacc;
    public double reltol_inacc;
    public int nitref;
    public int maxit;
    public int verbose;
    public int max_bk_iter;
    public double bk_scale;
    public double centrality;
    public int equilibrate;
    public int staticreg;

    public Settings() {
        this.gamma = GlblOpts.GAMMA;
        this.delta = GlblOpts.DELTA;
        this.eps = GlblOpts.EPS;
        this.feastol = GlblOpts.FEASTOL;
        this.abstol = GlblOpts.ABSTOL;
        this.reltol = GlblOpts.RELTOL;
        this.feastol_inacc = GlblOpts.FTOL_INACC;
        this.abstol_inacc = GlblOpts.ATOL_INACC;
        this.reltol_inacc = GlblOpts.RTOL_INACC;
        this.nitref = GlblOpts.NITREF;
        this.maxit = GlblOpts.MAXIT;
        this.verbose = GlblOpts.VERBOSE;
        this.max_bk_iter = GlblOpts.MAX_BK;
        this.bk_scale = GlblOpts.BK_SCALE;
        this.centrality = GlblOpts.CENTRALITY;
        this.equilibrate = GlblOpts.EQUILIBRATE;
        this.staticreg = GlblOpts.STATICREG;
    }
}
