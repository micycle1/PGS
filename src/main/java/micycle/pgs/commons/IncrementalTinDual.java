package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;

import micycle.pgs.color.RGB;
import processing.core.PConstants;
import processing.core.PShape;

/**
 * Produces a (barycentric) dual graph from a TinFour Delaunay Triangulation.
 * Triangle centroids form the vertices of the dual graph, and each
 * triangulation vertex corresponds to a face in the dual graph; the dual graph
 * has similarity to the voronoi diagram, but it isn't the same.
 * <p>
 * The dual of (fairly) regular triangulation is a regular polygonal (hex/oct)
 * mesh.
 * 
 * @author Michael Carleton
 *
 */
public class IncrementalTinDual {

	// TODO option to maintain original perimeter

	private final IIncrementalTin tin;
	private final boolean hasConstraints;
	private PShape mesh;

	/**
	 * List of all edges in the dual graph.
	 */
	public final List<DualEdge> edges;
	/**
	 * Maps triangulation edges to their dual edge.
	 */
	public final Map<IQuadEdge, DualEdge> edgeDuals;
	/**
	 * Maps triangulation vertices to their corresponding/surrounding face in the
	 * dual graph.
	 */
	public final Map<Vertex, PShape> vertexDuals;

	public IncrementalTinDual(IIncrementalTin tin) {
		this.tin = tin;
		hasConstraints = !tin.getConstraints().isEmpty();
		edges = new ArrayList<>();
		edgeDuals = new HashMap<>();
		vertexDuals = new HashMap<>();
		createDualEdges();
	}

	/**
	 * Compute list of dual edges.
	 */
	private void createDualEdges() {
		tin.getEdgeIterator().forEachRemaining(e -> { // iterates base edges only
			if (hasConstraints && !e.isConstrainedRegionInterior()) {
				// TODO keep unconstrained edges, but mark as unconstrained?
				return;
			}

			final SimpleTriangle a = new SimpleTriangle(tin, e);
			final SimpleTriangle b = new SimpleTriangle(tin, e.getDual());

			if (a.isGhost() || b.isGhost()) {
				// ignore dual edges whose triangle(s) lie outside the bounds of the
				// triangulation
				return;
			}

			/*
			 * Dualedges are created from the centroids of the two triangles that share a
			 * triangulation edge.
			 */
			final Vertex vA = centroid(a);
			final Vertex vB = centroid(b);
			final DualEdge dualEdge = new DualEdge(vA, vB, e);
			edges.add(dualEdge);
			edgeDuals.put(e, dualEdge);
		});
	}

	/**
	 * Generate mesh of polygonal dual faces.
	 */
	public PShape getMesh() {
		if (mesh == null) {
			mesh = new PShape(PConstants.GROUP);

			final HashSet<Vertex> seenHubs = new HashSet<>();

			edgeDuals.keySet().forEach(edge -> {
				if (hasConstraints && !edge.isConstrainedRegionInterior()) {
					return;
				}

				/*
				 * Compute dual face for both the edge and its dual. If using base edge only,
				 * then faces are sometimes skipped (see
				 * github.com/gwlucastrig/Tinfour/discussions/80).
				 */
				if (seenHubs.add(edge.getA())) {
					PShape face = dualFace(edge);
					if (face != null) {
						mesh.addChild(face);
					}
				}
				if (seenHubs.add(edge.getB())) {
					PShape face = dualFace(edge.getDual());
					if (face != null) {
						mesh.addChild(face);
					}
				}
			});
		}
		return mesh;
	}

	/**
	 * Creates a face from the dual graph that is the dual of a vertex in the
	 * triangulation (the A vertex of the given edge is used).
	 * 
	 * @param edge
	 * @return face shape, or null if no face could be made (such as if the edge
	 *         lies on exterior, etc.)
	 */
	private PShape dualFace(IQuadEdge edge) {
		// Use TreeSet to order face vertices by angle compared to central vertex.
		final Set<Vertex> faceVertices = new TreeSet<>((o1, o2) -> {
			final double a1 = atan2Quick(o1.y - edge.getA().y, o1.x - edge.getA().x);
			final double a2 = atan2Quick(o2.y - edge.getA().y, o2.x - edge.getA().x);
			return a1 > a2 ? 1 : -1;
		});

		/*
		 * Use pinwheel on each edge (around each hub vertex) to find the dualedges of
		 * the face that encloses the hub (a faster alternative to JTS Polygonizer).
		 */
		for (IQuadEdge p : edge.pinwheel()) {
			final DualEdge dual = edgeDuals.get(p.getBaseReference());
			if (dual == null) {
				return null; // no face to be made
			}
			faceVertices.add(dual.a);
			faceVertices.add(dual.b);
		}

		final PShape dualFace = new PShape(PShape.GEOMETRY);
		dualFace.setFill(true);
		dualFace.setFill(RGB.composeColor(0, 150, 200, 128));
		dualFace.setStroke(true);
		dualFace.setStroke(RGB.composeColor(255, 255, 255));
		dualFace.setStrokeWeight(2);
		dualFace.beginShape();
		faceVertices.forEach(v -> dualFace.vertex((float) v.getX(), (float) v.getY()));
		dualFace.endShape(PConstants.CLOSE);
		vertexDuals.put(edge.getA(), dualFace);
		return dualFace;
	}

	/**
	 * Computes the centroid/barycentre of a triangle.
	 */
	private static Vertex centroid(final SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new Vertex(x, y, 0);
	}

	private static double atan2Quick(final double y, final double x) {
		final double THREE_QRTR_PI = Math.PI * 0.75;
		final double QRTR_PI = Math.PI * 0.25;

		double r, angle;
		final double abs_y = Math.abs(y) + 1e-10f; // kludge to prevent 0/0 condition

		if (x < 0.0f) {
			r = (x + abs_y) / (abs_y - x); // (3)
			angle = THREE_QRTR_PI; // (4)
		} else {
			r = (x - abs_y) / (x + abs_y); // (1)
			angle = QRTR_PI; // (2)
		}
		angle += (0.1963f * r * r - 0.9817f) * r; // (2 | 4)
		if (y < 0.0f)
			return (-angle); // negate if in quad III or IV
		else
			return (angle);
	}

	public class DualEdge {

		public final Vertex a, b;
		/** The DualEdge instance is the dual of this quadedge */
		final IQuadEdge e;

		private DualEdge(Vertex a, Vertex b, IQuadEdge e) {
			this.a = a;
			this.b = b;
			this.e = e;
		}

		@Override
		public int hashCode() {
			return Float.floatToIntBits((float) (b.y + a.y)) ^ Float.floatToIntBits((float) (b.x + a.x - 1));
		}

		@Override
		public boolean equals(Object obj) {
			if (obj instanceof DualEdge) {
				DualEdge other = (DualEdge) obj;
				return (other.a.equals(a) && other.b.equals(b)) || (other.a.equals(b) && other.b.equals(a));
			}
			return false;
		}

	}

}
