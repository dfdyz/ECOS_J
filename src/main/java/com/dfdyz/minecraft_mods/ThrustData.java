package com.dfdyz.minecraft_mods;

import org.joml.Vector3d;

public class ThrustData implements Thrust {
    public final Vector3d pos;
    public final Vector3d dir;
    public final LimitType limitType;

    public void setTheta(double theta) {
        this.theta = theta;
        this.cosTheta = Math.cos(theta);
    }

    protected double theta;        // angle limit (radians), only used for ANGLE_LIMIT
    protected double cosTheta;
    public double Fmax;
    public final Vector3d force = new Vector3d();

    protected ThrustData(Vector3d pos, Vector3d dir, LimitType limitType, double Fmax) {
        this.pos = pos;
        this.dir = dir;
        this.limitType = limitType;
        this.Fmax = Fmax;
        this.theta = Math.toRadians(30);
    }

    protected ThrustData(Vector3d pos, Vector3d dir, double Fmax, double theta) {
        this.pos = pos;
        this.dir = dir;
        this.limitType = LimitType.ANGLE_LIMIT;
        this.Fmax = Fmax;
        this.theta = theta;
    }

    public static ThrustData Static(Vector3d pos, Vector3d dir, double Fmax){
        return new ThrustData(pos, dir, LimitType.STATIC, Fmax);
    }

    public static ThrustData HalfBall(Vector3d pos, Vector3d dir, double Fmax){
        return new ThrustData(pos, dir, LimitType.HALF_BALL, Fmax);
    }

    public static ThrustData AngleLimit(Vector3d pos, Vector3d dir, double Fmax, double theta){
        return new ThrustData(pos, dir, Fmax, theta);
    }

    @Override public Vector3d  getPos()       { return pos; }
    @Override public Vector3d  getDir()       { return dir; }
    @Override public LimitType getLimitType() { return limitType; }

    @Override
    public double getCosTheta() {return cosTheta;}
    @Override public double    getFmax()      { return Fmax; }
    @Override public void      setFmax(double f) { Fmax = f; }
    @Override public Vector3d  getForce()     { return force; }
}
