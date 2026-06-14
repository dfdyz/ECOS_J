package com.dfdyz.minecraft_mods;

import org.joml.Vector3d;
import java.util.ArrayList;
import java.util.List;

public class ShipData implements Ship {
    public final Vector3d center = new Vector3d();
    public final Vector3d targetF = new Vector3d();
    public final Vector3d targetT = new Vector3d();
    public final List<Thrust> thrusts = new ArrayList<>();

    @Override public Vector3d getCenter()  { return center; }
    @Override public Vector3d getTargetF() { return targetF; }
    @Override public Vector3d getTargetT() { return targetT; }
    @Override public List<Thrust> getThrusts() { return thrusts; }
}
