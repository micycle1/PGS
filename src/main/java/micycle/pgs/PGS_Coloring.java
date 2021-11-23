package micycle.pgs;

import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.jgrapht.alg.color.ColorRefinementAlgorithm;
import org.jgrapht.alg.color.LargestDegreeFirstColoring;
import org.jgrapht.alg.color.RandomGreedyColoring;
import org.jgrapht.alg.color.SaturationDegreeColoring;
import org.jgrapht.alg.color.SmallestDegreeLastColoring;
import org.jgrapht.alg.interfaces.VertexColoringAlgorithm.Coloring;
import org.jgrapht.graph.AbstractBaseGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.noding.BasicSegmentString;
import org.locationtech.jts.noding.SegmentString;
import micycle.pgs.utility.PEdge;
import micycle.pgs.color.RGB;
import micycle.pgs.utility.RLFColoring;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Coloring mesh-like shapes / shapes that have neighbouring parts via
 * graph-coloring technique.
 * 
 * Assigns colors to the faces of a mesh such that no two adjacent faces are of
 * the same color.
 * 
 * Intelligently colors meshes (mesh-like shapes), such that no two adjacent
 * faces are of the same color.
 * 
 * PGS_MeshColoring ? PGS_Colorizer ?
 * 
 * <p>
 * The smallest number of colors needed to color a graph G is called the
 * chromatic number
 * 
 * <p>
 * The methods in this class distinguish between mesh-like shapes and
 * non-mesh-like shapes. In mesh-like shapes, and adjacent cells not only share
 * edges, but the edges are identical!; that is, an edge from each has identical
 * coordinates (such as a triangulation). This distinction is necessary because
 * non-mesh-like shapes require a single step of pre-processing ("noding") to
 * find neighbouring cells whose neighbouring edges are not identical (i.e. they
 * overlap but are not equal). In none-mesh-like shape, lines map overlap, but
 * they may not have identical start and points.
 * 
 * In non-mesh shapes, edges adjacent cells may be shared to some extent (but
 * not the whole length).
 * 
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Coloring {

	private PGS_Coloring() {
	}

	public enum ColoringAlgorithm {
		/**
		 * The greedy coloring algorithm with a random vertex ordering.
		 */
		RANDOM,
		/**
		 * The largest degree first greedy coloring algorithm.
		 * 
		 * <p>
		 * This is the greedy coloring algorithm which orders the vertices by
		 * non-increasing degree.
		 */
		LARGEST_DEGREE_FIRST, // aka Largest Degree Ordering ?
		/**
		 * The smallest degree last greedy coloring algorithm.
		 * 
		 * This is the greedy coloring algorithm with the smallest-last ordering of the
		 * vertices.
		 * 
		 */
		SMALLEST_DEGREE_LAST,
		/**
		 * This is the greedy coloring algorithm using saturation degree ordering. The
		 * saturation degree of a vertex is defined as the number of different colors to
		 * which it is adjacent. The algorithm selects always the vertex with the
		 * largest saturation degree. If multiple vertices have the same maximum
		 * saturation degree, a vertex of maximum degree in the uncolored subgraph is
		 * selected.
		 */
		DSATUR,
		/**
		 * Finds the coarsest coloring of a graph (via color refinement algorithm).
		 */
		COARSE,
		/**
		 * Recursive largest-first (recommended).
		 */
		RLF
	}

	// TODO color point set by embedding as triangulatation then compute coloring of
	// triangulation
	private static void colorPoints(List<PVector> points) {
		// color a gabriel graph?
	}

	/**
	 * 
	 * @param meshShape GROUP PShape, whose children make up a mesh (more formally,
	 *                  a planar straight-line graph) (share line segments)
	 * @return
	 */
	public static Map<PShape, Integer> colorMesh(PShape meshShape, ColoringAlgorithm coloringAlgorithm) {
//		if (meshShape.getFamily() != PConstants.GROUP) {
////			System.err.println("colorMesh(): Input shape is not a group.");
//			return new HashMap<>();
//		}
		return colorMesh(PGS_Conversion.getChildren(meshShape), coloringAlgorithm);
	}

	/**
	 * 
	 * @param shapes a collection of shapes representing a mesh (more formally, a
	 *               planar straight-line graph).
	 * @return Get the color map; a mapping from each face to its color class
	 * @see #colorMesh(PShape, int[])
	 */
	public static Map<PShape, Integer> colorMesh(Collection<PShape> shapes, ColoringAlgorithm coloringAlgorithm) {
		final AbstractBaseGraph<PShape, DefaultEdge> graph = prepareGraph(shapes);
		final Coloring<PShape> coloring;

		switch (coloringAlgorithm) {
			case RANDOM : // randomly ordered sequential
				coloring = new RandomGreedyColoring<>(graph, new Random()).getColoring();
				break;
			case SMALLEST_DEGREE_LAST :
				/*
				 * LDO orders the vertices in decreasing order of degree, the idea being that
				 * the large degree vertices can be colored more easily.
				 */
				coloring = new SmallestDegreeLastColoring<>(graph).getColoring();
				break;
			case LARGEST_DEGREE_FIRST :
				coloring = new LargestDegreeFirstColoring<>(graph).getColoring(); // Greedy Welsh-Powell
				break;
			case DSATUR : // Degree of Saturation Algorithm (DSATUR)
				/*
				 * SDO (saturation degree ordering) is a variant on LDO where the vertices are
				 * ordered in decreasing order by "saturation degree", defined as the number of
				 * distinct colors in the vertex neighborhood. NOTE IS BETTER THAN RANDOM
				 */
				coloring = new SaturationDegreeColoring<>(graph).getColoring();
				break;
			case COARSE :
				coloring = new ColorRefinementAlgorithm<>(graph).getColoring();
				break;
			case RLF :
				coloring = new RLFColoring<>(graph).getColoring();
				break;
			default :
				return new HashMap<>();
		}
		System.out.println(coloring.getNumberColors());
		return coloring.getColors();
	}

	/**
	 * Applies the coloring (given by an array of colors) to the mesh-like shape.
	 * This method is void; it mutates the fill colour of the input shape.
	 * 
	 * @param shape
	 * @param colorPalette
	 */
	public static void colorMesh(PShape shape, ColoringAlgorithm coloringAlgorithm, int[] colorPalette) {
		PGS_Coloring.colorMesh(shape, coloringAlgorithm).forEach((face, color) -> {
			int c = colorPalette[color % colorPalette.length]; // NOTE remove for now to help debugging
			face.setFill(c);
		});
	}

	/**
	 * 
	 * @param shape
	 * @param coloringAlgorithm
	 * @param colorPalette      ["#FFFFFF", "#...", etc.]
	 */
	public static void colorMesh(PShape shape, ColoringAlgorithm coloringAlgorithm, String[] colorPalette) {
		colorMesh(shape, coloringAlgorithm, hexToColor(colorPalette));
	}

	public static Map<PShape, Integer> colorNonMesh(PShape shape, ColoringAlgorithm coloringAlgorithm) {
		final PShape mesh = toPShape(nodeNonMesh(shape));
		return colorMesh(mesh, coloringAlgorithm);
	}

	/**
	 * NOTE unlike {@link #colorMesh(PShape, ColoringAlgorithm, int[])} this method
	 * returns a PShape, since it doesn't apply coloring to the non-mesh input shape
	 * (rather a processed mesh version of it)
	 * 
	 * @param shape
	 * @param coloringAlgorithm
	 * @param colors
	 * @return
	 */
	public static PShape colorNonMesh(PShape shape, ColoringAlgorithm coloringAlgorithm, int[] colors) {
		final PShape mesh = toPShape(nodeNonMesh(shape));
		colorMesh(mesh, coloringAlgorithm, colors);
		PGS_Conversion.setAllStrokeColor(mesh, RGB.WHITE, 2);
		return mesh;
	}

	public static PShape colorNonMesh(PShape shape, ColoringAlgorithm coloringAlgorithm, String[] colors) {
		final PShape mesh = toPShape(nodeNonMesh(shape));
		colorMesh(mesh, coloringAlgorithm, colors);
		PGS_Conversion.setAllStrokeColor(mesh, RGB.WHITE, 2);
		return mesh;
	}

	/**
	 * Generates a (dual) graph representing the mesh (where graph vertices
	 * represent mesh faces and graph edges represent a shared edge between faces).
	 * 
	 * @param meshFaces
	 * @return dual-graph of the mesh
	 */
	private static AbstractBaseGraph<PShape, DefaultEdge> prepareGraph(Collection<PShape> meshFaces) {
		final AbstractBaseGraph<PShape, DefaultEdge> graph = new SimpleGraph<>(DefaultEdge.class);
		// map of which edge belong to each face; used to detect half-edges
		final HashMap<PEdge, PShape> edgesMap = new HashMap<>();

		for (PShape face : meshFaces) {
			graph.addVertex(face); // always add child so disconnected shapes are colored
			for (int i = 0; i < face.getVertexCount(); i++) {
				final PVector a = face.getVertex(i);
				final PVector b = face.getVertex((i + 1) % face.getVertexCount());
				if (a.equals(b)) {
					continue;
				}
				final PEdge e = new PEdge(a, b);
				final PShape neighbour = edgesMap.get(e);

				if (neighbour != null) {
					// edge seen before, so faces must be adjacent; create edge between faces
					if (neighbour.equals(face)) { // probably bad input (3 edges the same)
						System.err.println("PGS_Coloring: Bad input — saw the same edge 3 times.");
						continue; // continue to prevent self-loop in graph
					}
					graph.addEdge(neighbour, face);
				} else {
					edgesMap.put(e, face); // edge is new
				}
			}
		}
		return graph;
	}

	/**
	 * Prepare a non-mesh shape into a graphable shape by "noding" it. This
	 * essentially means splitting edges into two at points where they intersect
	 * (touch) another edge.
	 * 
	 * @param shape a GROUP PShape
	 * @return
	 */
	private static Collection<Geometry> nodeNonMesh(PShape shape) {
		final List<SegmentString> segmentStrings = new ArrayList<>(shape.getChildCount() * 3);

		for (PShape face : shape.getChildren()) {
			for (int i = 0; i < face.getVertexCount(); i++) {
				final PVector a = face.getVertex(i);
				final PVector b = face.getVertex((i + 1) % face.getVertexCount());
				if (!a.equals(b)) {
					segmentStrings.add(new BasicSegmentString(new Coordinate[] { PGS.coordFromPVector(a), PGS.coordFromPVector(b) }, null));
				}
			}
		}
		return PGS.polygonizeSegments(segmentStrings);
	}

	public static PShape nodeShape(PShape shape) {
		return toPShape(nodeNonMesh(shape));
	}

	/**
	 * Converts hex color strings to Processing integer colors (RRGGBB).
	 * 
	 * @param colors list of hex e.g "021d34" or "#FF6F91"
	 * @return
	 */
	static int[] hexToColor(String[] colors) {
		int[] out = new int[colors.length];
		for (int i = 0; i < colors.length; i++) {
			if (colors[i].charAt(0) == '#') { // handle leading '#'
				colors[i] = colors[i].substring(1);
			}
			out[i] = -16777216 + (int) (Long.parseLong(colors[i], 16));
		}
		return out;
	}

}
