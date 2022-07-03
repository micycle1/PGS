package micycle.pgs;

import java.util.ArrayList;
import java.util.Collection;
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
import org.locationtech.jts.noding.SegmentString;

import micycle.pgs.color.RGB;
import micycle.pgs.commons.GeneticColoring;
import micycle.pgs.commons.RLFColoring;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Intelligently color meshes (or mesh-like shapes) such that no two adjacent
 * faces have the same color, while minimising the number of colors used.
 * <p>
 * The methods in this class distinguish between mesh-like shapes (<i>conforming
 * meshes</i>) and non-mesh-like shapes (<i>non-conforming meshes</i>). This
 * distinction is necessary because shapes that represent non-conforming meshes
 * require a single step of pre-processing ("noding") to first split edges
 * before coloring. The difference is described below:
 * 
 * <p style="margin-left: 40px">
 * <i>Conforming Meshes</i> : Consists of adjacent cells that not only share
 * edges, but every pair of shared edges are <b>identical</b> (having the same
 * coordinates) (such as a triangulation). <br>
 * <i>Non-Conforming Meshes</i> : Consists of adjacent cells that share edges
 * (i.e. edges may overlap) but adjacent edges do not necessarily have identical
 * start and end coordinates.
 * </p>
 * 
 * @author Michael Carleton
 * @since 1.2.0
 *
 */
public final class PGS_Coloring {

	private PGS_Coloring() {
	}

	/**
	 * Specifies the algorithm/heuristic used by the underlying graph coloring process to find
	 * a coloring for mesh faces. RLF, followed by DSATUR generally produce the
	 * "best" colorings (as measured by chromatic number, where lower is better).
	 */
	public enum ColoringAlgorithm {
		/**
		 * The greedy coloring algorithm with a random vertex ordering.
		 */
		RANDOM,
		/**
		 * The largest degree first greedy coloring algorithm.
		 * 
		 * <p>
		 * This algorithm orders the vertices in decreasing order of degree, the idea
		 * being that the large degree vertices can be colored more easily.
		 */
		LARGEST_DEGREE_FIRST, // aka Largest Degree Ordering ?
		/**
		 * The smallest degree last greedy coloring algorithm.
		 * 
		 * This is the greedy coloring algorithm with the smallest-last ordering of the
		 * vertices.
		 */
		SMALLEST_DEGREE_LAST,
		/**
		 * DSATUR (saturation degree ordering) is a variant on Largest Degree Ordering
		 * where the vertices are ordered in decreasing order by "saturation degree",
		 * defined as the number of distinct colors in the vertex neighborhood.
		 */
		DSATUR,
		/**
		 * Finds the coarsest coloring of a graph.
		 */
		COARSE,
		/**
		 * Recursive largest-first coloring (recommended).
		 */
		RLF,
		/**
		 * Finds a coloring using a genetic algorithm. Unlike all other algorithms this
		 * specifically targets a chromaticity of 4 (falls back to 5 if no solution is
		 * found).
		 */
		GENETIC
	}

	/**
	 * Computes a coloring of the given mesh shape, returning a color class for each
	 * mesh face.
	 * 
	 * @param meshShape         a GROUP PShape, whose children constitute the faces
	 *                          of a <b>conforming</b> mesh
	 * @param coloringAlgorithm coloring algorithm used to color the mesh
	 * @return a color-class map; a mapping of each mesh face (PShape) -> color
	 *         class (integer)
	 */
	public static Map<PShape, Integer> colorMesh(PShape meshShape, ColoringAlgorithm coloringAlgorithm) {
		return colorMesh(PGS_Conversion.getChildren(meshShape), coloringAlgorithm);
	}

	/**
	 * Computes a coloring of the given mesh shape, returning a color class for each
	 * mesh face.
	 * 
	 * @param shapes            a collection of shapes that constitute the faces of
	 *                          a <b>conforming</b> mesh
	 * @param coloringAlgorithm coloring algorithm used to color the mesh
	 * @return a color-class map; a mapping of each mesh face (PShape) -> color
	 *         class (integer)
	 */
	public static Map<PShape, Integer> colorMesh(Collection<PShape> shapes, ColoringAlgorithm coloringAlgorithm) {
		final Coloring<PShape> coloring = findColoring(shapes, coloringAlgorithm);
		return coloring.getColors();
	}

	/**
	 * Computes a coloring of the given mesh shape and colors its faces using the
	 * colors provided. This method mutates the fill colour of the input shape.
	 * 
	 * @param shape             a GROUP PShape, whose children constitute the faces
	 *                          of a <b>conforming</b> mesh
	 * @param coloringAlgorithm coloring algorithm used to color the mesh
	 * @param colorPalette      the (integer) colors with which to color the mesh
	 * @return the input shape (whose faces have now been colored)
	 */
	public static PShape colorMesh(PShape shape, ColoringAlgorithm coloringAlgorithm, int[] colorPalette) {
		final Coloring<PShape> coloring = findColoring(shape, coloringAlgorithm);
		if (coloring.getNumberColors() > colorPalette.length) {
			System.err.format("WARNING: Number of mesh colors (%s) exceeds those provided in palette (%s)%s", coloring.getNumberColors(),
					colorPalette.length, System.lineSeparator());
		}
		coloring.getColors().forEach((face, color) -> {
			int c = colorPalette[color % colorPalette.length]; // NOTE use modulo to avoid OOB exception
			face.setFill(c);
		});

		return shape; // return original shape (to keep method signature the same as colorNonMesh())
	}

	/**
	 * Computes a coloring of the given mesh shape and colors its faces using the
	 * colors provided. This method mutates the fill colour of the input shape.
	 * 
	 * @param shape             a GROUP PShape, whose children constitute the faces
	 *                          of a <b>conforming</b> mesh
	 * @param coloringAlgorithm coloring algorithm used to color the mesh
	 * @param colorPalette      the string colors (e.g. "#FFFFFF", or "cba5e8") with
	 *                          which to color the mesh
	 * @return the input shape (whose faces have now been colored)
	 */
	public static PShape colorMesh(PShape shape, ColoringAlgorithm coloringAlgorithm, String[] colorPalette) {
		return colorMesh(shape, coloringAlgorithm, RGB.hexToColor(colorPalette));
	}

	/**
	 * Computes a coloring of the given non-conforming mesh shape, returning a color
	 * class for each face of the pre-processed (noded) mesh.
	 * 
	 * @param shape             a GROUP PShape, whose children constitute the faces
	 *                          of a <b>non-conforming</b> mesh
	 * @param coloringAlgorithm coloring algorithm used to color the mesh
	 * @return a color-class map; a mapping of each noded mesh face (PShape) ->
	 *         color class (integer)
	 */
	public static Map<PShape, Integer> colorNonMesh(PShape shape, ColoringAlgorithm coloringAlgorithm) {
		final PShape mesh = nodeNonMesh(shape);
		return colorMesh(mesh, coloringAlgorithm);
	}

	/**
	 * Computes a coloring of the given non-conforming mesh shape and colors the
	 * faces of its noded representation using the colors provided.
	 * 
	 * @param shape             a GROUP PShape, whose children constitute the faces
	 *                          of a <b>non-conforming</b> mesh
	 * @param coloringAlgorithm coloring algorithm used to color the mesh
	 * @param colorPalette      the (integer) colors with which to color the mesh
	 * @return noded representation of the input shape (whose faces have now been
	 *         colored)
	 */
	public static PShape colorNonMesh(PShape shape, ColoringAlgorithm coloringAlgorithm, int[] colorPalette) {
		final PShape mesh = nodeNonMesh(shape);
		colorMesh(mesh, coloringAlgorithm, colorPalette);
		PGS_Conversion.setAllStrokeColor(mesh, RGB.WHITE, 2);
		return mesh;
	}

	/**
	 * Computes a coloring of the given non-conforming mesh shape and colors the
	 * faces of its noded representation using the colors provided.
	 * 
	 * @param shape             a GROUP PShape, whose children constitute the faces
	 *                          of a <b>non-conforming</b> mesh
	 * @param coloringAlgorithm coloring algorithm used to color the mesh
	 * @param colorPalette      the string colors (e.g. "#FFFFFF", or "cba5e8") with
	 *                          which to color the mesh
	 * @return noded representation of the input shape (whose faces have now been
	 *         colored)
	 */
	public static PShape colorNonMesh(PShape shape, ColoringAlgorithm coloringAlgorithm, String[] colorPalette) {
		final PShape mesh = nodeNonMesh(shape);
		colorMesh(mesh, coloringAlgorithm, colorPalette);
		PGS_Conversion.setAllStrokeColor(mesh, RGB.WHITE, 2);
		return mesh;
	}

	private static Coloring<PShape> findColoring(PShape meshShape, ColoringAlgorithm coloringAlgorithm) {
		return findColoring(PGS_Conversion.getChildren(meshShape), coloringAlgorithm);
	}

	/**
	 * Finds a coloring for the graphable/mesh-like shape (as given by a collection
	 * of faces) using the coloring algorithm specified.
	 */
	private static Coloring<PShape> findColoring(Collection<PShape> shapes, ColoringAlgorithm coloringAlgorithm) {
		final AbstractBaseGraph<PShape, DefaultEdge> graph = PGS_Conversion.toDualGraph(shapes);
		final Coloring<PShape> coloring;

		switch (coloringAlgorithm) {
			case RANDOM : // randomly ordered sequential
				coloring = new RandomGreedyColoring<>(graph, new Random()).getColoring();
				break;
			case SMALLEST_DEGREE_LAST :
				coloring = new SmallestDegreeLastColoring<>(graph).getColoring();
				break;
			case LARGEST_DEGREE_FIRST :
				coloring = new LargestDegreeFirstColoring<>(graph).getColoring(); // Greedy Welsh-Powell
				break;
			case DSATUR :
				coloring = new SaturationDegreeColoring<>(graph).getColoring();
				break;
			case COARSE :
				coloring = new ColorRefinementAlgorithm<>(graph).getColoring();
				break;
			case GENETIC :
				coloring = new GeneticColoring<>(graph).getColoring();
				break;
			case RLF :
			default :
				coloring = new RLFColoring<>(graph).getColoring();
		}
		return coloring;
	}

	/**
	 * Converts a non-conforming mesh shape into a conforming mesh by "noding" it.
	 * This essentially means splitting edges into two at points where they
	 * intersect (touch) another edge.
	 * 
	 * @param shape a GROUP PShape
	 * @return the input shape, having been noded and polygonized
	 */
	private static PShape nodeNonMesh(PShape shape) {
		final List<SegmentString> segmentStrings = new ArrayList<>(shape.getChildCount() * 3);

		for (PShape face : shape.getChildren()) {
			for (int i = 0; i < face.getVertexCount(); i++) {
				final PVector a = face.getVertex(i);
				final PVector b = face.getVertex((i + 1) % face.getVertexCount());
				if (!a.equals(b)) {
					segmentStrings.add(PGS.createSegmentString(a, b));
				}
			}
		}
		return PGS.polygonizeSegments(segmentStrings, true);
	}

}
