package com.dfdyz.minecraft_mods;

import org.joml.Vector3d;

import java.util.List;

public interface Ship {
    Vector3d getCenter();
    Vector3d getTargetF();
    Vector3d getTargetT();
    List<Thrust> getThrusts();
}
