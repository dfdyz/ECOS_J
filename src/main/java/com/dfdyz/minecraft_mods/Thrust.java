package com.dfdyz.minecraft_mods;

import org.joml.Vector3d;

public interface Thrust {
    Vector3d getPos();
    Vector3d getDir();
    LimitType getLimitType();
    double getCosTheta();         // max angle (radians) for ANGLE_LIMIT; ignored otherwise
    double getFmax();
    void setFmax(double fmax);
    Vector3d getForce();
}
