package com.dfdyz.minecraft_mods;

import org.joml.Vector3d;

public class Utils {


    public static String getStr(double[] arr) {
        StringBuilder stringBuilder = new StringBuilder().append("[");
        for (int i = 0; i < arr.length; i++) {
            stringBuilder.append(String.format("%3f", arr[i]));
            if (i != arr.length - 1) {
                stringBuilder.append(", ");
            }
        }
        return stringBuilder.append("]").toString();
    }

    public static String getStr(double[][] arr) {
        StringBuilder stringBuilder = new StringBuilder().append("[");
        for (int i = 0; i < arr.length; i++) {
            stringBuilder.append(getStr(arr[i]));
            if (i != arr.length - 1) {
                stringBuilder.append(",\n");
            }
        }
        return stringBuilder.append("]").toString();
    }

    public static String toStr(Vector3d vector) {
        return String.format("( %.5f, %.5f, %.5f )", vector.x, vector.y, vector.z);
    }

}
