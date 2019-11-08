using System.IO;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Postprocessing
{
    /// <summary>
	/// Enumeration containing different types of T-Spline shell elements.
	/// </summary>
	public enum TSplineShellType
	{
		Linear, Section, Thickness
	}

	/// <summary>
	/// Paraview file  generator for T-Spline shell geometries.
	/// </summary>
	public class ParaviewTsplineShells
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
		public ParaviewTsplineShells(Model model, IVectorView solution, string filename)
		{
			_model = model;
			_solution = solution;
			_filename = filename;
		}

		/// <summary>
		/// Creates Paraview File of the T-Splines shells geometry.
		/// </summary>
		public void CreateParaviewFile(TSplineShellType shellType = TSplineShellType.Linear)
		{
			var projectiveControlPoints = CalculateProjectiveControlPoints();
			var numberOfPointsPerElement = 4;
			var nodes = new double[_model.Elements.Count * numberOfPointsPerElement, 3];
			var pointIndex = 0;
			foreach (var element in _model.Elements)
			{
				var tsplineElement = element as TSplineKirchhoffLoveShellElement;
				var elementPoints = tsplineElement.CalculatePointsForPostProcessing(tsplineElement);

				for (int i = 0; i < elementPoints.GetLength(0); i++)
				{
					nodes[pointIndex, 0] = elementPoints[i, 0];
					nodes[pointIndex, 1] = elementPoints[i, 1];
					nodes[pointIndex++, 2] = elementPoints[i, 2];
				}
			}

			var elementConnectivity = CreateTsplineConnectivity();

			var pointDisplacements = new double[nodes.GetLength(0), 3];
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
					localDisplacements[counterCP++, 2] =
						(!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, StructuralDof.TranslationZ))
							? 0.0
							: _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, StructuralDof.TranslationZ]];
				}
				var elementKnotDisplacements = element.ElementType.CalculateDisplacementsForPostProcessing(element, localDisplacements);
				for (int i = 0; i < elementConnectivity.GetLength(1); i++)
				{
					var knotConnectivity = elementConnectivity[element.ID, i];
					pointDisplacements[knotConnectivity, 0] = elementKnotDisplacements[i, 0];
					pointDisplacements[knotConnectivity, 1] = elementKnotDisplacements[i, 1];
					pointDisplacements[knotConnectivity, 2] = elementKnotDisplacements[i, 2];
				}
			}

			WriteTSplineShellsFile(nodes, elementConnectivity, pointDisplacements);
		}

		private void WriteTSplineShellsFile(double[,] nodeCoordinates, int[,] elementConnectivity, double[,] displacements)
		{
			var numberOfPoints = nodeCoordinates.GetLength(0);
			var numberOfCells = elementConnectivity.GetLength(0);

			const int numberOfVerticesPerCell = 4;
			const int paraviewCellCode = 9;

			using (StreamWriter outputFile = new StreamWriter($"..\\..\\..\\MGroup.IGA.Tests\\OutputFiles\\{_filename}Paraview.vtu"))
			{
				outputFile.WriteLine("<?xml version=\"1.0\"?>");
				outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">");
				outputFile.WriteLine("<UnstructuredGrid>");

				outputFile.WriteLine($"<Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
				outputFile.WriteLine($"<PointData Vectors=\"U\">");
				outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"U\" format=\"ascii\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{displacements[i, 0]} {displacements[i, 1]} {displacements[i, 2]}");

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</PointData>");
				outputFile.WriteLine("<Points>");
				outputFile.WriteLine("<DataArray type=\"Float32\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{nodeCoordinates[i, 0]} {nodeCoordinates[i, 1]} {nodeCoordinates[i, 2]}");

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Points>");
				outputFile.WriteLine("<Cells>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"connectivity\">");

				for (int i = 0; i < numberOfCells; i++)
				{
					for (int j = 0; j < elementConnectivity.GetLength(1); j++)
						outputFile.Write($"{elementConnectivity[i, j]} ");
					outputFile.WriteLine("");
				}

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"offsets\">");

				var offset = 0;
				for (int i = 0; i < numberOfCells; i++)
				{
					offset += numberOfVerticesPerCell;
					outputFile.WriteLine(offset);
				}

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"types\">");

				for (int i = 0; i < numberOfCells; i++)
					outputFile.WriteLine(paraviewCellCode);

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Cells>");
				outputFile.WriteLine("</Piece>");
				outputFile.WriteLine("</UnstructuredGrid>");
				outputFile.WriteLine("</VTKFile>");
			}
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

		private int[,] CreateTsplineConnectivity()
		{
			var elementConnectivity = new int[_model.Elements.Count, 4];
			var pointCounter = 0;
			for (int i = 0; i < _model.Elements.Count; i++)
			{
				elementConnectivity[i, 0] = pointCounter++;
				elementConnectivity[i, 1] = pointCounter++;
				elementConnectivity[i, 2] = pointCounter++;
				elementConnectivity[i, 3] = pointCounter++;
			}

			return elementConnectivity;
		}
	}
}
