package micycle.pts;

import org.locationtech.jts.geom.Geometry;
import org.twak.camp.Corner;
import org.twak.camp.Edge;
import org.twak.camp.Machine;
import org.twak.camp.Output.Face;
import org.twak.camp.Skeleton;
import org.twak.utils.collections.Loop;
import org.twak.utils.collections.LoopL;

import processing.core.PShape;

/**
 * https://github.com/twak/campskeleton/blob/wiki/headless.md
 * 
 * @author MCarleton
 *
 */
public class SSkeleton {

	public static void test(PShape shape) {

		Corner c1 = new Corner(0, 0), c2 = new Corner(100, -100), c3 = new Corner(100, 0);

		LoopL<Corner> corners = new LoopL<Corner>();

		Machine speed1 = new Machine(Math.PI / 4);
		Machine speed2 = new Machine(Math.PI / 3);

		Edge e1 = new Edge(c1, c2), e2 = new Edge(c2, c3), e3 = new Edge(c3, c1);

		corners.add(new Loop<Corner>(c1, c2, c3));

		Loop<Edge> loop1 = new Loop<Edge>();
		loop1.append(e1);
		loop1.append(e2);
		loop1.append(e3);

		e1.machine = speed1;
		e2.machine = speed1;
		e3.machine = speed2;

		Skeleton skel = new Skeleton(loop1.singleton(), true);
//		Skeleton skel = new Skeleton(corners);
		skel.skeleton();

//		JTS: create a GeometryCollection of the polygons that you wish to join and then call .union() on it.

		for (Face face : skel.output.faces.values()) {
			System.out.println("face:");
			for (Loop<javax.vecmath.Point3d> lp3 : face.points)
				for (javax.vecmath.Point3d pt : lp3)
					System.out.println(pt);
		}
	}

//	private static Loop<Edge> edgesFromGeometry(Geometry geometry) {
//		
//	}

}
