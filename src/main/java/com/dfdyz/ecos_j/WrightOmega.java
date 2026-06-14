package com.dfdyz.ecos_j;

public class WrightOmega {

    public static double wrightOmega(double z) {
        double w = 0.0;
        double r = 0.0;
        double q = 0.0;
        double zi = 0.0;

        if (z < 0.0) return -1;

        if (z < 1.0 + Math.PI) {
            q = z - 1;
            r = q;
            w = 1 + 0.5 * r;
            r *= q;
            w += 1.0 / 16.0 * r;
            r *= q;
            w -= 1.0 / 192.0 * r;
            r *= q;
            w -= 1.0 / 3072.0 * q;
            r *= q;
            w += 13.0 / 61440.0 * q;
        } else {
            r = Math.log(z);
            q = r;
            zi = 1.0 / z;
            w = z - r;
            q = r * zi;
            w += q;
            q = q * zi;
            w += q * (0.5 * r - 1);
            q = q * zi;
            w += q * (1.0 / 3.0 * r * r - 3.0 / 2.0 * r + 1);
        }
        r = z - w - Math.log(w);
        z = (1 + w);
        q = z + 2.0 / 3.0 * r;
        w *= 1 + r / z * (z * q - 0.5 * r) / (z * q - r);
        r = (2 * w * w - 8 * w - 1) / (72.0 * (z * z * z * z * z * z)) * r * r * r * r;
        z = (1 + w);
        q = z + 2.0 / 3.0 * r;
        w *= 1 + r / z * (z * q - 0.5 * r) / (z * q - r);
        r = (2 * w * w - 8 * w - 1) / (72.0 * (z * z * z * z * z * z)) * r * r * r * r;

        return w;
    }
}
