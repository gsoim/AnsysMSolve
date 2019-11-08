using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Postprocessing
{
    /// <summary>
	/// Paraview file  generator for 3D NURBS geometries.
	/// </summary>
	public class ParaviewNurbs3D
	{
		private readonly string _filename;
		private readonly Model _model;
		private readonly IVectorView _solution;

		/// <summary>
		/// Defines a Paraview FileWriter.
		/// </summary>
		/// <param name="model">An isogeometric <see cref="Model"/>.</param>
		/// <param name="solution">An <see cref="IVectorView"/> containing the solution of the linear system.</param>
		/// <param name="filename">The name of the paraview file to be generated.</param>
		public ParaviewNurbs3D(Model model, IVectorView solution, string filename)
		{
			_model = model;
			_solution = solution;
			_filename = filename;
		}

		/// <summary>
		/// Calculates the cartesian coordinate of a parametric point.
		/// </summary>
		/// <param name="numberOfCPKsi">Number of <see cref="ControlPoint"/> in the parametric axis Ksi.</param>
		/// <param name="degreeKsi">Degree of splines for the parametric axis Ksi.</param>
		/// <param name="knotValueVectorKsi"> Knot Value Vector for the parametric axis Ksi.</param>
		/// <param name="numberOfCPHeta">Number of <see cref="ControlPoint"/> in the parametric axis Heta.</param>
		/// <param name="degreeHeta">Degree of splines for the parametric axis Heta.</param>
		/// <param name="knotValueVectorHeta"> Knot Value Vector for the parametric axis Heta.</param>
		/// <param name="numberOfCPZeta">Number of <see cref="ControlPoint"/> in the parametric axis Zeta.</param>
		/// <param name="degreeZeta">Degree of splines for the parametric axis Zeta.</param>
		/// <param name="knotValueVectorZeta"> Knot Value Vector for the parametric axis Zeta.</param>
		/// <param name="projectiveControlPointCoordinates">A <see cref="double"/> two dimensional array containing the projective coordinates of the control points.</param>
		/// <param name="ksiCoordinate">Parametric coordinate Ksi of the physical point to be evaluated.</param>
		/// <param name="hetaCoordinate">Parametric coordinate Heta of the physical point to be evaluated.</param>
		/// <param name="zetaCoordinate">Parametric coordinate Zeta of the physical point to be evaluated.</param>
		/// <returns>A <see cref="Vector"/> containing the physical coordinates of the solid point calculated.</returns>
		public static Vector SolidPoint3D(int numberOfCPKsi, int degreeKsi, Vector knotValueVectorKsi,
			int numberOfCPHeta, int degreeHeta, Vector knotValueVectorHeta, int numberOfCPZeta, int degreeZeta,
			Vector knotValueVectorZeta, double[,] projectiveControlPointCoordinates, double ksiCoordinate, double hetaCoordinate,
			double zetaCoordinate)
		{
			var spanKsi = ParaviewNurbs2D.FindSpan(numberOfCPKsi, degreeKsi, ksiCoordinate, knotValueVectorKsi);
			var spanHeta = ParaviewNurbs2D.FindSpan(numberOfCPHeta, degreeHeta, hetaCoordinate, knotValueVectorHeta);
			var spanZeta = ParaviewNurbs2D.FindSpan(numberOfCPZeta, degreeZeta, zetaCoordinate, knotValueVectorZeta);

			var pointFunctionsKsi = ParaviewNurbs2D.BasisFunctions(spanKsi, ksiCoordinate, degreeKsi, knotValueVectorKsi);
			var pointFunctionsHeta = ParaviewNurbs2D.BasisFunctions(spanHeta, hetaCoordinate, degreeHeta, knotValueVectorHeta);
			var pointFunctionsZeta = ParaviewNurbs2D.BasisFunctions(spanZeta, zetaCoordinate, degreeZeta, knotValueVectorZeta);

			var cartesianPoint = Vector.CreateZero(4);
			var indexKsi = spanKsi - degreeKsi;

			for (int k = 0; k <= degreeZeta; k++)
			{
				var indexZeta = spanZeta - degreeZeta + k;
				for (int j = 0; j <= degreeHeta; j++)
				{
					var indexHeta = spanHeta - degreeHeta + j;
					for (int i = 0; i <= degreeKsi; i++)
					{
						var cpIndex = (indexKsi + i) * (numberOfCPHeta + 1) * (numberOfCPZeta + 1) +
									  indexHeta * (numberOfCPZeta + 1) + indexZeta;

						var cpCoordinates = Vector.CreateFromArray(new double[]
						{
							projectiveControlPointCoordinates[cpIndex, 0],
							projectiveControlPointCoordinates[cpIndex, 1],
							projectiveControlPointCoordinates[cpIndex, 2],
							projectiveControlPointCoordinates[cpIndex, 3]
						});
						cpCoordinates.Scale(pointFunctionsKsi[i]);
						cpCoordinates.Scale(pointFunctionsHeta[j]);
						cpCoordinates.Scale(pointFunctionsZeta[k]);
						cartesianPoint = cartesianPoint + cpCoordinates;
					}
				}
			}

			return cartesianPoint;
		}

		/// <summary>
		/// Creates Paraview File of the 3D NURBS geometry.
		/// </summary>
		public void CreateParaviewFile()
		{
			var uniqueKnotsKsi = _model.PatchesDictionary[0].KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity();
			var uniqueKnotsHeta = _model.PatchesDictionary[0].KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity();
			var uniqueKnotsZeta = _model.PatchesDictionary[0].KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity();

			var numberOfKnotsKsi = uniqueKnotsKsi[0].Length;
			var numberOfKnotsHeta = uniqueKnotsHeta[0].Length;
			var numberOfKnotsZeta = uniqueKnotsZeta[0].Length;

			var knots = new double[numberOfKnotsKsi * numberOfKnotsHeta * numberOfKnotsZeta, 3];
			var count = 0;
			var patch = _model.PatchesDictionary[0];

			var projectiveControlPoints = CalculateProjectiveControlPoints();

			for (int knotKsiIndex = 0; knotKsiIndex < numberOfKnotsKsi; knotKsiIndex++)
			{
				for (int knotHetaIndex = 0; knotHetaIndex < numberOfKnotsHeta; knotHetaIndex++)
				{
					for (int knotZetaIndex = 0; knotZetaIndex < numberOfKnotsZeta; knotZetaIndex++)
					{
						var ksiCoordinate = uniqueKnotsKsi[0][knotKsiIndex];
						var hetaCoordinate = uniqueKnotsHeta[0][knotHetaIndex];
						var zetaCoordinate = uniqueKnotsZeta[0][knotZetaIndex];

						var point3D = SolidPoint3D(
							patch.NumberOfControlPointsKsi - 1, patch.DegreeKsi, patch.KnotValueVectorKsi,
							patch.NumberOfControlPointsHeta - 1, patch.DegreeHeta, patch.KnotValueVectorHeta,
							patch.NumberOfControlPointsZeta - 1, patch.DegreeZeta, patch.KnotValueVectorZeta,
							projectiveControlPoints, ksiCoordinate, hetaCoordinate, zetaCoordinate);

						knots[count, 0] = point3D[0] / point3D[3];
						knots[count, 1] = point3D[1] / point3D[3];
						knots[count++, 2] = point3D[2] / point3D[3];
					}
				}
			}

			var numberOfElements = (numberOfKnotsKsi - 1) * (numberOfKnotsHeta - 1) * (numberOfKnotsZeta - 1);
			var elementConnectivity = new int[numberOfElements, 8];

			foreach (var element in _model.Elements)
				for (int i = 0; i < element.KnotsDictionary.Count; i++)
					elementConnectivity[element.ID, i] = element.Knots.ToArray()[i].ID;

			var knotDisplacements = new double[knots.GetLength(0), 3];
			foreach (var element in _model.Elements)
			{
				var localDisplacements = Matrix.CreateZero(element.ControlPointsDictionary.Count, 3);
				var counterCP = 0;
				foreach (var controlPoint in element.ControlPoints)
				{
					localDisplacements[counterCP, 0] =
						(!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationX))
							? 0.0
							: _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationX]];
					localDisplacements[counterCP, 1] =
						(!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationY))
							? 0.0
							: _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationY]];
					localDisplacements[counterCP, 2] =
						(!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationZ))
							? 0.0
							: _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationZ]];
					counterCP++;
				}
				var elementKnotDisplacements = element.ElementType.CalculateDisplacementsForPostProcessing(element, localDisplacements);
				for (int i = 0; i < elementConnectivity.GetLength(1); i++)
				{
					var knotConnectivity = elementConnectivity[element.ID, i];
					knotDisplacements[knotConnectivity, 0] = elementKnotDisplacements[i, 0];
					knotDisplacements[knotConnectivity, 1] = elementKnotDisplacements[i, 1];
					knotDisplacements[knotConnectivity, 2] = elementKnotDisplacements[i, 2];
				}
			}

			WriteParaviewFile3D(knots, knotDisplacements);
		}

		private double[,] CalculateProjectiveControlPoints()
		{
			var projectiveCPs = new double[_model.PatchesDictionary[0].ControlPoints.Count, 4];
			foreach (var controlPoint in _model.PatchesDictionary[0].ControlPoints)
			{
				var weight = controlPoint.WeightFactor;
				projectiveCPs[controlPoint.ID, 0] = controlPoint.X * weight;
				projectiveCPs[controlPoint.ID, 1] = controlPoint.Y * weight;
				projectiveCPs[controlPoint.ID, 2] = controlPoint.Z * weight;
				projectiveCPs[controlPoint.ID, 3] = weight;
			}

			return projectiveCPs;
		}

		private void WriteParaviewFile3D(double[,] nodeCoordinates, double[,] displacements)
		{
			var patch = _model.PatchesDictionary[0];

			var numberOfKnotsKsi = patch.KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0].Length;
			var numberOfKnotsHeta = patch.KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0].Length;
			var numberOfKnotZeta = patch.KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0].Length;

			var numberOfNodes = nodeCoordinates.GetLength(0);

			using (StreamWriter outputFile = new StreamWriter($"..\\..\\..\\MGroup.IGA.Tests\\OutputFiles\\{_filename}Paraview.vts"))
			{
				outputFile.WriteLine("<?xml version=\"1.0\"?>");
				outputFile.WriteLine("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\" >");
				outputFile.WriteLine($"<StructuredGrid  WholeExtent=\"{0} {numberOfKnotZeta - 1} {0} {numberOfKnotsHeta - 1} {0} {numberOfKnotsKsi - 1}\">");
				outputFile.WriteLine($"<Piece Extent=\"{0} {numberOfKnotZeta - 1} {0} {numberOfKnotsHeta - 1} {0} {numberOfKnotsKsi - 1}\">");

				outputFile.WriteLine("<PointData Vectors=\"Disp\"  >");
				outputFile.WriteLine("<DataArray type=\"Float32\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">");
				for (int i = 0; i < numberOfNodes; i++)
					outputFile.WriteLine($"{displacements[i, 0]} {displacements[i, 1]} {displacements[i, 2]}");
				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</PointData>");

				outputFile.WriteLine("<Celldata>");
				outputFile.WriteLine("</Celldata>");
				outputFile.WriteLine("<Points>");
				outputFile.WriteLine("<DataArray type=\"Float32\" Name=\"Array\" NumberOfComponents=\"3\" format=\"ascii\">");
				for (int i = 0; i < numberOfNodes; i++)
					outputFile.WriteLine($"{nodeCoordinates[i, 0]} {nodeCoordinates[i, 1]} {nodeCoordinates[i, 2]}");
				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Points>");

				outputFile.WriteLine("</Piece>");
				outputFile.WriteLine("</StructuredGrid>");
				outputFile.WriteLine("</VTKFile>");
			}
		}
	}
}
