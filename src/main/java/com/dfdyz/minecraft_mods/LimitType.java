package com.dfdyz.minecraft_mods;

public enum LimitType {
    STATIC,       // fixed direction, output = α · d
    HALF_BALL,    // v · d ≥ 0 (hemisphere of the original direction)
    ANGLE_LIMIT   // v · d ≥ ||v|| · cos(θ),  max angle = θ
}
