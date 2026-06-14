package com.dfdyz.ecos_j;

public class TimerOps {
    public static void tic(Timer t) {
        t.tic = System.nanoTime();
        t.freq = 1000000000L;
    }

    public static double toc(Timer t) {
        t.toc = System.nanoTime();
        return ((double)(t.toc - t.tic) / (double)t.freq);
    }
}
