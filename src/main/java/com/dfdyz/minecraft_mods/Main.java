package com.dfdyz.minecraft_mods;


public class Main {
    public static void main(String[] args) {

// 初始化（一次）
        double[][] r = {{0,0,1}, {0,0,-1}};  // 喷口位置
        double[][] d = {{0,0,1}, {0,0,-1}};  // 原始朝向
        double[]   Fmax = {10, 10};           // 最大推力

// 高频循环（零分配）
        while (true) {
            //double[] Ft = getNewForceTarget();
            //double[] Tt = getNewTorqueTarget();
           // double[][] v = alloc.solve(Ft, Tt);  // → v[i][3] 为各喷口推力矢量
        }

    }
}