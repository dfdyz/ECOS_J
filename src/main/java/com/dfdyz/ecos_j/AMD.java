package com.dfdyz.ecos_j;

public class AMD {

    public static final int EMPTY = -1;
    public static final int AMD_CONTROL = 5;
    public static final int AMD_INFO = 20;
    public static final int AMD_DENSE = 0;
    public static final int AMD_AGGRESSIVE = 1;
    public static final double AMD_DEFAULT_DENSE = 10.0;
    public static final int AMD_DEFAULT_AGGRESSIVE = 1;
    public static final int AMD_STATUS = 0;
    public static final int AMD_N = 1;
    public static final int AMD_NZ = 2;
    public static final int AMD_SYMMETRY = 3;
    public static final int AMD_NZDIAG = 4;
    public static final int AMD_NZ_A_PLUS_AT = 5;
    public static final int AMD_NDENSE = 6;
    public static final int AMD_MEMORY = 7;
    public static final int AMD_NCMPA = 8;
    public static final int AMD_LNZ = 9;
    public static final int AMD_NDIV = 10;
    public static final int AMD_NMULTSUBS_LDL = 11;
    public static final int AMD_NMULTSUBS_LU = 12;
    public static final int AMD_DMAX = 13;
    public static final int AMD_OK = 0;
    public static final int AMD_OUT_OF_MEMORY = -1;
    public static final int AMD_INVALID = -2;
    public static final int AMD_OK_BUT_JUMBLED = 1;

    static int FLIP(int i) { return (-(i) - 2); }
    static int UNFLIP(int i) { return (i < EMPTY) ? FLIP(i) : (i); }

    public static void amd_defaults(double[] Control) {
        if (Control != null) {
            for (int i = 0; i < AMD_CONTROL; i++) { Control[i] = 0; }
            Control[AMD_DENSE] = AMD_DEFAULT_DENSE;
            Control[AMD_AGGRESSIVE] = AMD_DEFAULT_AGGRESSIVE;
        }
    }

    public static void amd_info(double[] Info) {
        if (Info == null) return;
        double n = Info[AMD_N];
        double ndiv = Info[AMD_NDIV];
        double nmultsubs_ldl = Info[AMD_NMULTSUBS_LDL];
        double nmultsubs_lu = Info[AMD_NMULTSUBS_LU];
        double lnz = Info[AMD_LNZ];
        double lnzd = (n >= 0 && lnz >= 0) ? (n + lnz) : (-1);
        if (GlblOpts.PRINTLEVEL > 2) {
            GlblOpts.printText("\nAMD version 2.3.1, Jun 20, 2012, results:\n");
            GlblOpts.printText("    status: OK\n");
            if (n >= 0) GlblOpts.printText("    n, dimension of A:                                  %.20g\n", n);
            GlblOpts.printText("    nonzeros in pattern of A+A' (excl. diagonal):       %.20g\n", Info[AMD_NZ_A_PLUS_AT]);
            GlblOpts.printText("    memory used, in bytes:                              %.20g\n", Info[AMD_MEMORY]);
            GlblOpts.printText("    # of memory compactions:                            %.20g\n", Info[AMD_NCMPA]);
            if (lnz >= 0) GlblOpts.printText("    nonzeros in L (excluding diagonal):                 %.20g\n", lnz);
            if (lnzd >= 0) GlblOpts.printText("    nonzeros in L (including diagonal):                 %.20g\n", lnzd);
            if (ndiv >= 0) GlblOpts.printText("    # divide operations for LDL' or LU:                 %.20g\n", ndiv);
            if (nmultsubs_ldl >= 0) GlblOpts.printText("    # multiply-subtract operations for LDL':            %.20g\n", nmultsubs_ldl);
            GlblOpts.printText("\n");
        }
    }

    public static int amd_order(int n, int[] Ap, int[] Ai, int[] P, double[] Control, double[] Info) {
        int[] Len, S = null;
        int nz, i, status, infoFlag;
        long nzaat, slen;
        double mem = 0;

        infoFlag = (Info != null) ? 1 : 0;
        if (infoFlag != 0) {
            for (i = 0; i < AMD_INFO; i++) { Info[i] = EMPTY; }
            Info[AMD_N] = n;
            Info[AMD_STATUS] = AMD_OK;
        }

        if (Ai == null || Ap == null || P == null || n < 0) {
            if (infoFlag != 0) Info[AMD_STATUS] = AMD_INVALID;
            return AMD_INVALID;
        }
        if (n == 0) { return AMD_OK; }

        nz = Ap[n];
        if (infoFlag != 0) { Info[AMD_NZ] = nz; }
        if (nz < 0) {
            if (infoFlag != 0) Info[AMD_STATUS] = AMD_INVALID;
            return AMD_INVALID;
        }

        status = amd_valid(n, n, Ap, Ai);
        if (status == AMD_INVALID) {
            if (infoFlag != 0) Info[AMD_STATUS] = AMD_INVALID;
            return AMD_INVALID;
        }

        Len = new int[n];
        int[] Pinv = new int[n];
        mem += n;
        mem += n;

        int[] Cp, Ci;
        int[] Rp = null, Ri = null;
        if (status == AMD_OK_BUT_JUMBLED) {
            Rp = new int[n + 1];
            Ri = new int[Math.max(nz, 1)];
            mem += (n + 1);
            mem += Math.max(nz, 1);
            amd_preprocess(n, Ap, Ai, Rp, Ri, Len, Pinv);
            Cp = Rp;
            Ci = Ri;
        } else {
            Cp = Ap;
            Ci = Ai;
        }

        nzaat = amd_aat(n, Cp, Ci, Len, P, Info);

        slen = nzaat;
        slen += nzaat / 5;
        for (i = 0; i < 7; i++) { slen += n; }
        mem += slen;

        if (slen < Integer.MAX_VALUE) {
            S = new int[(int) slen];
        }

        if (S == null) {
            Rp = null; Ri = null;
            Len = null; Pinv = null;
            if (infoFlag != 0) Info[AMD_STATUS] = AMD_OUT_OF_MEMORY;
            return AMD_OUT_OF_MEMORY;
        }

        if (infoFlag != 0) {
            Info[AMD_MEMORY] = mem * 4; // sizeof(int) = 4
        }

        amd_1(n, Cp, Ci, P, Pinv, Len, (int) slen, S, Control, Info);

        Rp = null; Ri = null;
        Len = null; Pinv = null;
        S = null;
        if (infoFlag != 0) Info[AMD_STATUS] = status;
        return status;
    }

    static long amd_aat(int n, int[] Ap, int[] Ai, int[] Len, int[] Tp, double[] Info) {
        int p1, p2, p, i, j, pj, pj2, k, nzdiag, nzboth, nz;
        double sym;
        long nzaat;

        if (Info != null) {
            for (i = 0; i < AMD_INFO; i++) { Info[i] = EMPTY; }
            Info[AMD_STATUS] = AMD_OK;
        }

        for (k = 0; k < n; k++) { Len[k] = 0; }

        nzdiag = 0;
        nzboth = 0;
        nz = Ap[n];

        for (k = 0; k < n; k++) {
            p1 = Ap[k];
            p2 = Ap[k + 1];
            for (p = p1; p < p2; ) {
                j = Ai[p];
                if (j < k) {
                    Len[j]++;
                    Len[k]++;
                    p++;
                } else if (j == k) {
                    p++;
                    nzdiag++;
                    break;
                } else {
                    break;
                }
                pj2 = Ap[j + 1];
                for (pj = Tp[j]; pj < pj2; ) {
                    i = Ai[pj];
                    if (i < k) {
                        Len[i]++;
                        Len[j]++;
                        pj++;
                    } else if (i == k) {
                        pj++;
                        nzboth++;
                        break;
                    } else {
                        break;
                    }
                }
                Tp[j] = pj;
            }
            Tp[k] = p;
        }

        for (j = 0; j < n; j++) {
            for (pj = Tp[j]; pj < Ap[j + 1]; pj++) {
                i = Ai[pj];
                Len[i]++;
                Len[j]++;
            }
        }

        if (nz == nzdiag) { sym = 1; }
        else { sym = (2 * (double) nzboth) / ((double) (nz - nzdiag)); }

        nzaat = 0;
        for (k = 0; k < n; k++) { nzaat += Len[k]; }

        if (Info != null) {
            Info[AMD_STATUS] = AMD_OK;
            Info[AMD_N] = n;
            Info[AMD_NZ] = nz;
            Info[AMD_SYMMETRY] = sym;
            Info[AMD_NZDIAG] = nzdiag;
            Info[AMD_NZ_A_PLUS_AT] = nzaat;
        }

        return nzaat;
    }

    static int amd_valid(int n_row, int n_col, int[] Ap, int[] Ai) {
        int nz, j, p1, p2, ilast, i, p, result = AMD_OK;
        if (n_row < 0 || n_col < 0 || Ap == null || Ai == null) { return AMD_INVALID; }
        nz = Ap[n_col];
        if (Ap[0] != 0 || nz < 0) { return AMD_INVALID; }
        for (j = 0; j < n_col; j++) {
            p1 = Ap[j];
            p2 = Ap[j + 1];
            if (p1 > p2) { return AMD_INVALID; }
            ilast = EMPTY;
            for (p = p1; p < p2; p++) {
                i = Ai[p];
                if (i < 0 || i >= n_row) { return AMD_INVALID; }
                if (i <= ilast) { result = AMD_OK_BUT_JUMBLED; }
                ilast = i;
            }
        }
        return result;
    }

    static void amd_preprocess(int n, int[] Ap, int[] Ai, int[] Rp, int[] Ri, int[] W, int[] Flag) {
        int i, j, p, p2;
        for (i = 0; i < n; i++) { W[i] = 0; Flag[i] = EMPTY; }
        for (j = 0; j < n; j++) {
            p2 = Ap[j + 1];
            for (p = Ap[j]; p < p2; p++) {
                i = Ai[p];
                if (Flag[i] != j) { W[i]++; Flag[i] = j; }
            }
        }
        Rp[0] = 0;
        for (i = 0; i < n; i++) { Rp[i + 1] = Rp[i] + W[i]; }
        for (i = 0; i < n; i++) { W[i] = Rp[i]; Flag[i] = EMPTY; }
        for (j = 0; j < n; j++) {
            p2 = Ap[j + 1];
            for (p = Ap[j]; p < p2; p++) {
                i = Ai[p];
                if (Flag[i] != j) { Ri[W[i]++] = j; Flag[i] = j; }
            }
        }
    }

    static void amd_1(int n, int[] Ap, int[] Ai, int[] P, int[] Pinv,
                      int[] Len, int slen, int[] S, double[] Control, double[] Info) {
        int i, j, k, p, pfree, iwlen, pj, p1, p2, pj2;
        int[] Iw, Pe, Nv, Head, Elen, Degree, W, Sp, Tp;

        iwlen = slen - 6 * n;
        int sIdx = 0;
        Pe = new int[n]; System.arraycopy(S, sIdx, Pe, 0, n); sIdx += n;
        Nv = new int[n]; // used as Sp
        Head = new int[n];
        Elen = new int[n];
        Degree = new int[n];
        W = new int[n];
        Iw = new int[iwlen];

        Sp = Nv;
        Tp = W;
        pfree = 0;
        for (j = 0; j < n; j++) {
            Pe[j] = pfree;
            Sp[j] = pfree;
            pfree += Len[j];
        }

        for (k = 0; k < n; k++) {
            p1 = Ap[k];
            p2 = Ap[k + 1];
            for (p = p1; p < p2; ) {
                j = Ai[p];
                if (j < k) {
                    Iw[Sp[j]++] = k;
                    Iw[Sp[k]++] = j;
                    p++;
                } else if (j == k) { p++; break; }
                else { break; }
                pj2 = Ap[j + 1];
                for (pj = Tp[j]; pj < pj2; ) {
                    i = Ai[pj];
                    if (i < k) {
                        Iw[Sp[i]++] = j;
                        Iw[Sp[j]++] = i;
                        pj++;
                    } else if (i == k) { pj++; break; }
                    else { break; }
                }
                Tp[j] = pj;
            }
            Tp[k] = p;
        }

        for (j = 0; j < n; j++) {
            for (pj = Tp[j]; pj < Ap[j + 1]; pj++) {
                i = Ai[pj];
                Iw[Sp[i]++] = j;
                Iw[Sp[j]++] = i;
            }
        }

        amd_2(n, Pe, Iw, Len, iwlen, pfree, Nv, Pinv, P, Head, Elen, Degree, W, Control, Info);
    }

    static int clear_flag(int wflg, int wbig, int[] W, int n) {
        if (wflg < 2 || wflg >= wbig) {
            for (int x = 0; x < n; x++) { if (W[x] != 0) W[x] = 1; }
            wflg = 2;
        }
        return wflg;
    }

    static void amd_2(int n, int[] Pe, int[] Iw, int[] Len, int iwlen, int pfree,
                      int[] Nv, int[] Next, int[] Last, int[] Head, int[] Elen,
                      int[] Degree, int[] W, double[] Control, double[] Info) {
        int deg, degme, dext, lemax, e, elenme, eln, i, ilast, inext, j,
            jlast, jnext, k, knt1, knt2, knt3, lenj, ln, me, mindeg, nel, nleft,
            nvi, nvj, nvpiv, slenme, wbig, we, wflg, wnvi, ok, ndense, ncmpa,
            dense, aggressive;
        int hash;

        double f, r, ndiv, s, nms_lu, nms_ldl, dmax, alpha, lnz, lnzme;
        int p, p1, p2, p3, p4, pdst, pend, pj, pme, pme1, pme2, pn, psrc;

        lnz = 0;
        ndiv = 0;
        nms_lu = 0;
        nms_ldl = 0;
        dmax = 1;
        me = EMPTY;
        mindeg = 0;
        ncmpa = 0;
        nel = 0;
        lemax = 0;

        if (Control != null) {
            alpha = Control[AMD_DENSE];
            aggressive = (Control[AMD_AGGRESSIVE] != 0) ? 1 : 0;
        } else {
            alpha = AMD_DEFAULT_DENSE;
            aggressive = AMD_DEFAULT_AGGRESSIVE;
        }
        if (alpha < 0) { dense = n - 2; }
        else { dense = (int) (alpha * Math.sqrt(n)); }
        dense = Math.max(16, dense);
        dense = Math.min(n, dense);

        for (i = 0; i < n; i++) {
            Last[i] = EMPTY; Head[i] = EMPTY; Next[i] = EMPTY;
            Nv[i] = 1; W[i] = 1; Elen[i] = 0; Degree[i] = Len[i];
        }

        wbig = Integer.MAX_VALUE - n;
        wflg = clear_flag(0, wbig, W, n);

        ndense = 0;
        for (i = 0; i < n; i++) {
            deg = Degree[i];
            if (deg == 0) {
                Elen[i] = FLIP(1); nel++; Pe[i] = EMPTY; W[i] = 0;
            } else if (deg > dense) {
                ndense++; Nv[i] = 0; Elen[i] = EMPTY; nel++; Pe[i] = EMPTY;
            } else {
                inext = Head[deg];
                if (inext != EMPTY) Last[inext] = i;
                Next[i] = inext; Head[deg] = i;
            }
        }

        while (nel < n) {
            for (deg = mindeg; deg < n; deg++) { me = Head[deg]; if (me != EMPTY) break; }
            mindeg = deg;
            inext = Next[me];
            if (inext != EMPTY) Last[inext] = EMPTY;
            Head[deg] = inext;
            elenme = Elen[me];
            nvpiv = Nv[me];
            nel += nvpiv;

            Nv[me] = -nvpiv;
            degme = 0;

            if (elenme == 0) {
                pme1 = Pe[me];
                pme2 = pme1 - 1;
                for (p = pme1; p <= pme1 + Len[me] - 1; p++) {
                    i = Iw[p];
                    nvi = Nv[i];
                    if (nvi > 0) {
                        degme += nvi; Nv[i] = -nvi; Iw[++pme2] = i;
                        ilast = Last[i]; inext = Next[i];
                        if (inext != EMPTY) Last[inext] = ilast;
                        if (ilast != EMPTY) { Next[ilast] = inext; }
                        else { Head[Degree[i]] = inext; }
                    }
                }
            } else {
                p = Pe[me];
                pme1 = pfree;
                slenme = Len[me] - elenme;
                for (knt1 = 1; knt1 <= elenme + 1; knt1++) {
                    if (knt1 > elenme) { e = me; pj = p; ln = slenme; }
                    else { e = Iw[p++]; pj = Pe[e]; ln = Len[e]; }
                    for (knt2 = 1; knt2 <= ln; knt2++) {
                        i = Iw[pj++];
                        nvi = Nv[i];
                        if (nvi > 0) {
                            if (pfree >= iwlen) {
                                Pe[me] = p; Len[me] -= knt1;
                                if (Len[me] == 0) Pe[me] = EMPTY;
                                Pe[e] = pj; Len[e] = ln - knt2;
                                if (Len[e] == 0) Pe[e] = EMPTY;
                                ncmpa++;
                                for (j = 0; j < n; j++) {
                                    pn = Pe[j];
                                    if (pn >= 0) { Pe[j] = Iw[pn]; Iw[pn] = FLIP(j); }
                                }
                                psrc = 0; pdst = 0; pend = pme1 - 1;
                                while (psrc <= pend) {
                                    j = FLIP(Iw[psrc++]);
                                    if (j >= 0) {
                                        Iw[pdst] = Pe[j]; Pe[j] = pdst++;
                                        lenj = Len[j];
                                        for (knt3 = 0; knt3 <= lenj - 2; knt3++) { Iw[pdst++] = Iw[psrc++]; }
                                    }
                                }
                                p1 = pdst;
                                for (psrc = pme1; psrc <= pfree - 1; psrc++) { Iw[pdst++] = Iw[psrc]; }
                                pme1 = p1; pfree = pdst; pj = Pe[e]; p = Pe[me];
                            }
                            degme += nvi; Nv[i] = -nvi; Iw[pfree++] = i;
                            ilast = Last[i]; inext = Next[i];
                            if (inext != EMPTY) Last[inext] = ilast;
                            if (ilast != EMPTY) { Next[ilast] = inext; }
                            else { Head[Degree[i]] = inext; }
                        }
                    }
                    if (e != me) { Pe[e] = FLIP(me); W[e] = 0; }
                }
                pme2 = pfree - 1;
            }

            Degree[me] = degme; Pe[me] = pme1; Len[me] = pme2 - pme1 + 1;
            Elen[me] = FLIP(nvpiv + degme);
            wflg = clear_flag(wflg, wbig, W, n);

            for (pme = pme1; pme <= pme2; pme++) {
                i = Iw[pme]; eln = Elen[i];
                if (eln > 0) {
                    nvi = -Nv[i]; wnvi = wflg - nvi;
                    for (p = Pe[i]; p <= Pe[i] + eln - 1; p++) {
                        e = Iw[p]; we = W[e];
                        if (we >= wflg) { we -= nvi; }
                        else if (we != 0) { we = Degree[e] + wnvi; }
                        W[e] = we;
                    }
                }
            }

            for (pme = pme1; pme <= pme2; pme++) {
                i = Iw[pme]; p1 = Pe[i]; p2 = p1 + Elen[i] - 1; pn = p1;
                hash = 0; deg = 0;
                if (aggressive != 0) {
                    for (p = p1; p <= p2; p++) {
                        e = Iw[p]; we = W[e];
                        if (we != 0) {
                            dext = we - wflg;
                            if (dext > 0) { deg += dext; Iw[pn++] = e; hash += e; }
                            else { Pe[e] = FLIP(me); W[e] = 0; }
                        }
                    }
                } else {
                    for (p = p1; p <= p2; p++) {
                        e = Iw[p]; we = W[e];
                        if (we != 0) { dext = we - wflg; deg += dext; Iw[pn++] = e; hash += e; }
                    }
                }
                Elen[i] = pn - p1 + 1;
                p3 = pn; p4 = p1 + Len[i];
                for (p = p2 + 1; p < p4; p++) {
                    j = Iw[p]; nvj = Nv[j];
                    if (nvj > 0) { deg += nvj; Iw[pn++] = j; hash += j; }
                }
                if (Elen[i] == 1 && p3 == pn) {
                    Pe[i] = FLIP(me); nvi = -Nv[i];
                    degme -= nvi; nvpiv += nvi; nel += nvi;
                    Nv[i] = 0; Elen[i] = EMPTY;
                } else {
                    Degree[i] = Math.min(Degree[i], deg);
                    Iw[pn] = Iw[p3]; Iw[p3] = Iw[p1]; Iw[p1] = me;
                    Len[i] = pn - p1 + 1;
                    hash = (hash & 0x7FFFFFFF) % n;
                    j = Head[hash];
                    if (j <= EMPTY) { Next[i] = FLIP(j); Head[hash] = FLIP(i); }
                    else { Next[i] = Last[j]; Last[j] = i; }
                    Last[i] = hash;
                }
            }

            Degree[me] = degme;
            lemax = Math.max(lemax, degme);
            wflg += lemax;
            wflg = clear_flag(wflg, wbig, W, n);

            for (pme = pme1; pme <= pme2; pme++) {
                i = Iw[pme];
                if (Nv[i] < 0) {
                    hash = Last[i];
                    j = Head[hash];
                    if (j == EMPTY) { i = EMPTY; }
                    else if (j < EMPTY) { i = FLIP(j); Head[hash] = EMPTY; }
                    else { i = Last[j]; Last[j] = EMPTY; }

                    while (i != EMPTY && Next[i] != EMPTY) {
                        ln = Len[i]; eln = Elen[i];
                        for (p = Pe[i] + 1; p <= Pe[i] + ln - 1; p++) { W[Iw[p]] = wflg; }
                        jlast = i; j = Next[i];
                        while (j != EMPTY) {
                            ok = (Len[j] == ln) && (Elen[j] == eln) ? 1 : 0;
                            for (p = Pe[j] + 1; ok != 0 && p <= Pe[j] + ln - 1; p++) {
                                if (W[Iw[p]] != wflg) ok = 0;
                            }
                            if (ok != 0) {
                                Pe[j] = FLIP(i); Nv[i] += Nv[j]; Nv[j] = 0;
                                Elen[j] = EMPTY; j = Next[j]; Next[jlast] = j;
                            } else { jlast = j; j = Next[j]; }
                        }
                        wflg++; i = Next[i];
                    }
                }
            }

            p = pme1; nleft = n - nel;
            for (pme = pme1; pme <= pme2; pme++) {
                i = Iw[pme]; nvi = -Nv[i];
                if (nvi > 0) {
                    Nv[i] = nvi;
                    deg = Degree[i] + degme - nvi;
                    deg = Math.min(deg, nleft - nvi);
                    inext = Head[deg];
                    if (inext != EMPTY) Last[inext] = i;
                    Next[i] = inext; Last[i] = EMPTY; Head[deg] = i;
                    mindeg = Math.min(mindeg, deg); Degree[i] = deg;
                    Iw[p++] = i;
                }
            }

            Nv[me] = nvpiv;
            Len[me] = p - pme1;
            if (Len[me] == 0) { Pe[me] = EMPTY; W[me] = 0; }
            if (elenme != 0) { pfree = p; }

            if (Info != null) {
                f = nvpiv; r = degme + ndense;
                dmax = Math.max(dmax, f + r);
                lnzme = f * r + (f - 1) * f / 2; lnz += lnzme;
                ndiv += lnzme;
                s = f * r * r + r * (f - 1) * f + (f - 1) * f * (2 * f - 1) / 6;
                nms_lu += s; nms_ldl += (s + lnzme) / 2;
            }
        }

        if (Info != null) {
            f = ndense; dmax = Math.max(dmax, ndense);
            lnzme = (f - 1) * f / 2; lnz += lnzme;
            ndiv += lnzme;
            s = (f - 1) * f * (2 * f - 1) / 6; nms_lu += s; nms_ldl += (s + lnzme) / 2;
            Info[AMD_LNZ] = lnz; Info[AMD_NDIV] = ndiv;
            Info[AMD_NMULTSUBS_LDL] = nms_ldl; Info[AMD_NMULTSUBS_LU] = nms_lu;
            Info[AMD_NDENSE] = ndense; Info[AMD_DMAX] = dmax;
            Info[AMD_NCMPA] = ncmpa; Info[AMD_STATUS] = AMD_OK;
        }

        for (i = 0; i < n; i++) { Pe[i] = FLIP(Pe[i]); }
        for (i = 0; i < n; i++) { Elen[i] = FLIP(Elen[i]); }

        for (i = 0; i < n; i++) {
            if (Nv[i] == 0) {
                j = Pe[i];
                if (j == EMPTY) continue;
                while (Nv[j] == 0) { j = Pe[j]; }
                e = j;
                j = i;
                while (Nv[j] == 0) { int jnext2 = Pe[j]; Pe[j] = e; j = jnext2; }
            }
        }

        amd_postorder(n, Pe, Nv, Elen, W, Head, Next, Last);

        for (k = 0; k < n; k++) { Head[k] = EMPTY; Next[k] = EMPTY; }
        for (e = 0; e < n; e++) { k = W[e]; if (k != EMPTY) { Head[k] = e; } }
        nel = 0;
        for (k = 0; k < n; k++) {
            e = Head[k]; if (e == EMPTY) break;
            Next[e] = nel; nel += Nv[e];
        }
        for (i = 0; i < n; i++) {
            if (Nv[i] == 0) {
                e = Pe[i];
                if (e != EMPTY) { Next[i] = Next[e]; Next[e]++; }
                else { Next[i] = nel++; }
            }
        }
        for (i = 0; i < n; i++) { k = Next[i]; Last[k] = i; }
    }

    static void amd_postorder(int nn, int[] Parent, int[] Nv, int[] Fsize,
                              int[] Order, int[] Child, int[] Sibling, int[] Stack) {
        int i, j, k, parent, frsize, f, fprev, maxfrsize, bigfprev, bigf, fnext;

        for (j = 0; j < nn; j++) { Child[j] = EMPTY; Sibling[j] = EMPTY; }
        for (j = nn - 1; j >= 0; j--) {
            if (Nv[j] > 0) {
                parent = Parent[j];
                if (parent != EMPTY) { Sibling[j] = Child[parent]; Child[parent] = j; }
            }
        }

        for (i = 0; i < nn; i++) {
            if (Nv[i] > 0 && Child[i] != EMPTY) {
                fprev = EMPTY; maxfrsize = EMPTY; bigfprev = EMPTY; bigf = EMPTY;
                for (f = Child[i]; f != EMPTY; f = Sibling[f]) {
                    frsize = Fsize[f];
                    if (frsize >= maxfrsize) { maxfrsize = frsize; bigfprev = fprev; bigf = f; }
                    fprev = f;
                }
                fnext = Sibling[bigf];
                if (fnext != EMPTY) {
                    if (bigfprev == EMPTY) { Child[i] = fnext; }
                    else { Sibling[bigfprev] = fnext; }
                    Sibling[bigf] = EMPTY; Sibling[fprev] = bigf;
                }
            }
        }

        for (i = 0; i < nn; i++) { Order[i] = EMPTY; }
        k = 0;
        for (i = 0; i < nn; i++) {
            if (Parent[i] == EMPTY && Nv[i] > 0) {
                k = amd_post_tree(i, k, Child, Sibling, Order, Stack, nn);
            }
        }
    }

    static int amd_post_tree(int root, int k, int[] Child, int[] Sibling,
                             int[] Order, int[] Stack, int nn) {
        int f, head, h, i;
        head = 0; Stack[0] = root;
        while (head >= 0) {
            i = Stack[head];
            if (Child[i] != EMPTY) {
                for (f = Child[i]; f != EMPTY; f = Sibling[f]) { head++; }
                h = head;
                for (f = Child[i]; f != EMPTY; f = Sibling[f]) { Stack[h--] = f; }
                Child[i] = EMPTY;
            } else {
                head--; Order[i] = k++;
            }
        }
        return k;
    }
}
