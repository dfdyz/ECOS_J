import com.dfdyz.minecraft_mods.*;
import org.joml.Vector3d;
import java.util.Random;

public class MixedThrustTest {
    public static void main(String[] args) {
        ShipData ship = new ShipData();
        ship.thrusts.add(ThrustData.HalfBall(new Vector3d(1,0,0), new Vector3d(0,0,1), 10));
        ship.thrusts.add(ThrustData.HalfBall(new Vector3d(-1,0,0), new Vector3d(0,0,1), 10));
        ship.thrusts.add(ThrustData.Static(new Vector3d(0,1,0), new Vector3d(1,0,0), 8));
        ship.thrusts.add(ThrustData.Static(new Vector3d(0,-1,0), new Vector3d(0,1,0), 8));

        try (ThrustAllocator alloc = new ThrustAllocator(ship)) {
            ship.targetF.set(0,0,4); ship.targetT.set(0,0,0);
            alloc.solve();
            System.out.println("--- Test 1: F=(0,0,4), T=(0,0,0) exit=" + alloc.lastExitFlag + " ---");
            check("  ", ship);

            ship.targetF.set(3,2,0); ship.targetT.set(0,0,1);
            alloc.solve();
            System.out.println("--- Test 2: F=(3,2,0), T=(0,0,1) exit=" + alloc.lastExitFlag + " ---");
            check("  ", ship);

            ship.center.x = 1.0;
            ship.targetF.set(0,0,4); ship.targetT.set(0,0,2);
            alloc.solve();
            System.out.println("--- Test 3: CG=(1,0,0), F=(0,0,4), T=(0,0,2) exit=" + alloc.lastExitFlag + " ---");
            check("  ", ship);
            ship.center.x = 0;

            ship.center.set(0,0,0);
            ship.targetF.set(0,0,4); ship.targetT.set(0,0,1);
            ship.thrusts.get(0).getDir().set(0,1,1).normalize();
            alloc.solve();
            System.out.println("--- Test 4: dir changed, F=(0,0,4), T=(0,0,1) exit=" + alloc.lastExitFlag + " ---");
            check("  ", ship);
            ship.thrusts.get(0).getDir().set(0,0,1);

            System.out.println("=== Test 5: random mixed (50) ===");
            Random rng = new Random(42);
            ShipData ship2 = new ShipData();
            final double fmax = 1e6;

            ship2.thrusts.add(ThrustData.HalfBall(new Vector3d(3,1,2), new Vector3d(0,0,-1), fmax));
            ship2.thrusts.add(ThrustData.HalfBall(new Vector3d(-3,1,2), new Vector3d(0,0,-1), fmax));
            ship2.thrusts.add(ThrustData.HalfBall(new Vector3d(3,1,-2), new Vector3d(0,0,1), fmax));
            ship2.thrusts.add(ThrustData.HalfBall(new Vector3d(-3,1,-2), new Vector3d(0,0,1), fmax));

            ship2.thrusts.add(ThrustData.HalfBall(new Vector3d(3,-1,2), new Vector3d(0,0,-1), fmax));
            ship2.thrusts.add(ThrustData.HalfBall(new Vector3d(-3,-1,2), new Vector3d(0,0,-1), fmax));
            ship2.thrusts.add(ThrustData.HalfBall(new Vector3d(3,1,-2), new Vector3d(0,0,1), fmax));
            ship2.thrusts.add(ThrustData.HalfBall(new Vector3d(-3,1,-2), new Vector3d(0,0,1), fmax));

            /*for (int[] p : new int[][]{{1,1},{ -1,1},{1,-1},{ -1,-1}}) {
                ship2.thrusts.add(ThrustData.Static(new Vector3d(p[0],p[1],0), new Vector3d(0,p[1]>0?1:-1,0), 800));
                ship2.thrusts.add(ThrustData.Static(new Vector3d(p[0],0,p[1]), new Vector3d(0,0,p[1]), 800));
            }
            ship2.thrusts.add(ThrustData.Static(new Vector3d(1,0,0), new Vector3d(1,0,0), 800));
            ship2.thrusts.add(ThrustData.Static(new Vector3d(-1,0,0), new Vector3d(-1,0,0), 800));*/

            try(var alloc2 = new ThrustAllocator2(ship2)) {
                ship2.center.set(0,1,0);
                ship2.targetF.set(0,20,0);
                ship2.targetT.set(0,0,0);
                alloc2.solve();
                System.out.println("--- static test: flag=" + alloc2.lastExitFlag + " ---");
                check("  ", ship2);

                int opt = 0;
                long start = System.nanoTime();
                for (int iter = 0; iter < 50; iter++) {
                    ship2.center.set((rng.nextDouble()-0.5)*1, (rng.nextDouble()-0.5)*1, (rng.nextDouble()-0.5)*1);
                    ship2.targetF.set((rng.nextDouble()-0.5)*2, (rng.nextDouble()-0.5)*2, (rng.nextDouble()-0.5)*2);
                    ship2.targetT.set(0, (rng.nextDouble()-0.5)*2, (rng.nextDouble()-0.5)*2);
                    alloc2.solve();
                    if (alloc2.lastExitFlag == 0 || alloc2.lastExitFlag == 10) opt++;
                    else {
                        System.out.printf("> %d/50 flag=%d\n", iter, alloc2.lastExitFlag);
                        check("  ", ship2);
                    }
                }
                long elapsed = System.nanoTime() - start;
                System.out.printf(">>> Optimal: %d/50 | avg %.1f us/solve | total %.3fms%n", opt, elapsed/1000.0/50, elapsed/1000.0/1000.0);
            }
        }
        System.out.println("\n=== Done ===");
    }

    static void check(String p, Ship ship) {
        int N = ship.getThrusts().size();
        Vector3d sf = new Vector3d(), st = new Vector3d();
        double totalF = 0;
        for (int i = 0; i < N; i++) {
            Vector3d fv = ship.getThrusts().get(i).getForce();
            sf.add(fv);
            totalF += fv.length();
            Vector3d r = new Vector3d(ship.getThrusts().get(i).getPos()).sub(ship.getCenter());
            Vector3d t = new Vector3d(r).cross(fv);
            st.add(t);
        }
        var fErr = sf.sub(ship.getTargetF(), new Vector3d());
        var tErr = st.sub(ship.getTargetT(), new Vector3d());

        System.out.println(p + "sum |f_i|: " + totalF);
        System.out.println(p + "f  : " + Utils.toStr(ship.getTargetF()));
        System.out.println(p + "sf  : " + Utils.toStr(sf));
        System.out.println(p + "fErr: " + Utils.toStr(fErr));
        System.out.println(p + "t  : " + Utils.toStr(ship.getTargetT()));
        System.out.println(p + "st  : " + Utils.toStr(st));
        System.out.println(p + "tErr: " + Utils.toStr(tErr));
        boolean ok = fErr.length() < 1e-3 && tErr.length() < 1e-3;
        for (int i = 0; i < N; i++) {
            Thrust t = ship.getThrusts().get(i);
            Vector3d fv = t.getForce();
            if (t.getLimitType() != LimitType.STATIC) {
                if (fv.dot(t.getDir()) < -1e-3) { System.out.println(p + "hemisphere fail " + i); ok = false; }
            } else {
                double len = fv.length();
                if (len > 1e-10) {
                    double cos = fv.dot(t.getDir()) / len;
                    if (Math.abs(cos) < 1-1e-3) { System.out.println(p + "fixed dir fail " + i + " cos="+cos); ok = false; }
                }
            }
        }
        System.out.println(p + (ok ? "PASS" : "FAIL"));
    }
}
