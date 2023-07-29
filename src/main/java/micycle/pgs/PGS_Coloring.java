package micycle.pgs;

import java.util.Collection;
import java.util.Map;

import org.jgrapht.alg.color.ColorRefinementAlgorithm;
import org.jgrapht.alg.color.LargestDegreeFirstColoring;
import org.jgrapht.alg.color.RandomGreedyColoring;
import org.jgrapht.alg.color.SaturationDegreeColoring;
import org.jgrapht.alg.color.SmallestDegreeLastColoring;
import org.jgrapht.alg.interfaces.VertexColoringAlgorithm.Coloring;
import org.jgrapht.graph.AbstractBaseGraph;
import org.jgrapht.graph.DefaultEdge;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import micycle.pgs.color.Colors;
import micycle.pgs.color.ColorUtils;
import micycle.pgs.commons.GeneticColoring;
import micycle.pgs.commons.RLFColoring;
import processing.core.PShape;

/**
 * This class provides methods to color meshes and mesh-like shapes. It ensures
 * that no two adjacent faces share the same color while also minimizing the
 * total number of colors used.
 * <p>
 * This class differentiates between "conforming meshes" and "non-conforming
 * meshes".
 * 
 * <p style="margin-left: 40px">
 * <i>Conforming Meshes</i> : Consist of adjacent cells that share edges and
 * every pair of shared edges are identical, meaning they have the same
 * coordinates. An example of a conforming mesh is a triangulation. <br>
 * <i>Non-Conforming Meshes</i> : Consist of adjacent cells that share edges
 * (i.e. edges may overlap) but adjacent edges do not necessarily have identical
 * start and end coordinates.
 * </p>
 * For a non-conforming mesh, a pre-processing step called "noding" is required
 * to split edges before coloring.
 * 
 * @author Michael Carleton
 * @since 1.2.0
 */
public final class PGS_Coloring {

	private PGS_Coloring() {
	}

	/**
	 * Specifies the algorithm/heuristic used by the underlying graph coloring
	 * process to find a coloring for mesh faces. RLF, followed by DSATUR generally
	 * produce the "best" colorings (as measured by chromatic number, where lower is
	 * better).
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
		 * Repeatedly calls the recursive largest-first (RLF) algorithm until a
		 * 4-coloring is found. The operation will break after 250 attempts if a
		 * 4-coloring is still not found; in this case, the result from the final
		 * attempt is returned.
		 */
		RLF_BRUTE_FORCE_4COLOR,
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
			face.setFill(true); // just in case
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
		return colorMesh(shape, coloringAlgorithm, ColorUtils.hexToColor(colorPalette));
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
		final PShape mesh = PGS_Meshing.nodeNonMesh(shape);
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
		final PShape mesh = PGS_Meshing.nodeNonMesh(shape);
		colorMesh(mesh, coloringAlgorithm, colorPalette);
		PGS_Conversion.setAllStrokeColor(mesh, Colors.WHITE, 2);
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
		final PShape mesh = PGS_Meshing.nodeNonMesh(shape);
		colorMesh(mesh, coloringAlgorithm, colorPalette);
		PGS_Conversion.setAllStrokeColor(mesh, Colors.WHITE, 2);
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
		Coloring<PShape> coloring;

		switch (coloringAlgorithm) {
			case RANDOM : // randomly ordered sequential
				coloring = new RandomGreedyColoring<>(graph, new XoRoShiRo128PlusRandom()).getColoring();
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
			case RLF_BRUTE_FORCE_4COLOR :
				int iterations = 0;
				do {
					coloring = new RLFColoring<>(graph).getColoring();
					iterations++;
				} while (coloring.getNumberColors() > 4 && iterations < 250);
				break;
			case RLF :
			default :
				coloring = new RLFColoring<>(graph, 1337).getColoring(); // NOTE fixed seed of 1337
		}
		return coloring;
	}

}
