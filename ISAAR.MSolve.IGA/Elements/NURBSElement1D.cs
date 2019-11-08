using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using Element = ISAAR.MSolve.IGA.Entities.Element;

namespace ISAAR.MSolve.IGA.Elements
{
    /// <summary>
	/// An one-dimensional continuum element that utilizes Non-Uniform Rational B-Splines for shape functions.
	/// Authors: Dimitris Tsapetis.
	/// </summary>
	public class NurbsElement1D : Element, IStructuralIsogeometricElement
	{
		protected static readonly IDofType[] ControlPointDofTypes = { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private IDofType[][] _dofTypes;

		/// <summary>
		/// Retrieves the type of Finite Element used. Since the element is Isogeometric its type is defined as unknown.
		/// </summary>
		public CellType CellType { get; } = CellType.Unknown;

		/// <summary>
		/// Property that Polynomial degree of the NURBS shape functions.
		/// </summary>
		public int Degree { get; set; }

		#region IStructuralIsogeometricElement

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		/// <summary>
		/// Retrieves the number of Dimensions of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.OneD;

		/// <summary>
		/// Boolean property that determines whether the material used for this elements has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Calculates the forces applies to an <see cref="NurbsElement1D"/> due to <see cref="FEM.Entities.MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="loads">A list of <see cref="FEM.Entities.MassAccelerationLoad"/>. For more info see <seealso cref="FEM.Entities.MassAccelerationLoad"/></param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => Array.Empty<double>();

		/// <summary>
		/// Calculates displacements of knots for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements) =>
			Array.Empty<double>();

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements) =>
			Array.Empty<double>();

		/// <summary>
		/// This method calculates the Neumann boundary condition when applied to a one dimensional NURBS element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="neumann"><inheritdoc cref="NeumannBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann)
		{
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(element);
			Dictionary<int, double> neumannLoad = new Dictionary<int, double>();
			IList<ControlPoint> controlPoints = new List<ControlPoint>();

			CalculateEdgeControlPoints(element, edge, controlPoints);

			CalculatePressure1D(element, edge, neumann, controlPoints, gaussPoints, neumannLoad);
			return neumannLoad;
		}

		/// <summary>
		/// This method calculates the Neumann boundary condition when applied to a one dimensional NURBS element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure)
		{
			Contract.Requires(element != null, "The element cannot be null");
			Contract.Requires(edge != null, "The edge cannot be null");
			Contract.Requires(pressure != null, "The pressure BC cannot be null");

			var gaussPoints = CreateElementGaussPoints(element);
			var pressureLoad = new Dictionary<int, double>();
			IList<ControlPoint> controlPoints = new List<ControlPoint>();

			CalculateEdgeControlPoints(element, edge, controlPoints);

			CalculatePressure1D(element, edge, pressure, controlPoints, gaussPoints, pressureLoad);
			return pressureLoad;
		}

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsElement1D"/> as it refers to two-dimensional loads.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="NeumannBoundaryCondition"/> was applied to.</param>
		/// <param name="neumann">The <see cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann) => throw new NotSupportedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsElement1D"/> as it refers to two-dimensional loads.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="PressureBoundaryCondition"/> was applied to.</param>
		/// <param name="pressure">The <see cref="PressureBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure) => throw new NotSupportedException();

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState() => throw new NotImplementedException();

		/// <summary>
		/// Clear any saved material states of the element.
		/// </summary>
		public void ClearMaterialStresses() => throw new NotImplementedException();

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="NurbsElement1D"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="ControlPoint"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			Contract.Requires(element != null, "The element cannot be null");

			_dofTypes = new IDofType[element.Nodes.Count][];
			for (var i = 0; i < element.Nodes.Count; i++)
			{
				_dofTypes[i] = ControlPointDofTypes;
			}

			return _dofTypes;
		}

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="NurbsElement1D"/>.</returns>
		public IMatrix MassMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified() => throw new NotImplementedException();

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState() => throw new NotImplementedException();

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="NurbsElement1D"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element) => throw new NotImplementedException();

		private static void CalculateEdgeControlPoints(Element element, Edge edge, IList<ControlPoint> controlPoints)
		{
			foreach (var controlPoint in element.ControlPoints)
			{
				if (element.Patch.NumberOfDimensions == 2)
				{
					controlPoints.Add(new ControlPoint()
					{
						ID = (edge.ID < 2)
							? controlPoint.ID % element.Patch.NumberOfControlPointsHeta
							: controlPoint.ID / element.Patch.NumberOfControlPointsHeta,
						Ksi = controlPoint.Ksi,
						Heta = controlPoint.Heta,
						Zeta = controlPoint.Zeta,
						X = controlPoint.X,
						Y = controlPoint.Y,
						Z = controlPoint.Z,
						WeightFactor = controlPoint.WeightFactor,
					});
				}
				else
				{
					var id = FindAxisControlPointId3D(element, edge, controlPoint);
					controlPoints.Add(new ControlPoint()
					{
						ID = id,
						Ksi = controlPoint.Ksi,
						Heta = controlPoint.Heta,
						Zeta = controlPoint.Zeta,
						X = controlPoint.X,
						Y = controlPoint.Y,
						Z = controlPoint.Z,
						WeightFactor = controlPoint.WeightFactor,
					});
				}
			}
		}

		private static void CalculatePressure1D(
			Element element,
			Edge edge,
			NeumannBoundaryCondition neumann,
			IList<ControlPoint> controlPoints,
			IList<GaussLegendrePoint3D> gaussPoints,
			IDictionary<int, double> neumannLoad)
		{
			var nurbs = new Nurbs1D(element, controlPoints, edge);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				double xGaussPoint = 0;
				double yGaussPoint = 0;
				double zGaussPoint = 0;
				double jacobian1 = 0.0;
				double jacobian2 = 0.0;
				var elementControlPoints = element.ControlPointsDictionary.Values.ToArray();
				for (int k = 0; k < elementControlPoints.Length; k++)
				{
					xGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].X;
					yGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].Y;
					zGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].Z;
					jacobian1 += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
					jacobian2 += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				}

				double jacdet = Math.Sqrt(Math.Pow(jacobian1, 2) + Math.Pow(jacobian2, 2));
				var loadGaussPoint = neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint);

				for (int k = 0; k < element.ControlPointsDictionary.Count; k++)
				{
					if (element.Model.GlobalDofOrdering.GlobalFreeDofs.Contains(elementControlPoints[k],
						StructuralDof.TranslationX))
					{
						int dofIDX =
							element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[k],
								StructuralDof.TranslationX];
						if (neumannLoad.ContainsKey(dofIDX))
						{
							neumannLoad[dofIDX] += jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] *
												   loadGaussPoint[0];
						}
						else
						{
							neumannLoad.Add(
								dofIDX,
								jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPoint[0]);
						}

					}

					if (element.Model.GlobalDofOrdering.GlobalFreeDofs.Contains(
						elementControlPoints[k],
						StructuralDof.TranslationY))
					{
						var dofIDY =
							element.Model.GlobalDofOrdering.GlobalFreeDofs[
								elementControlPoints[k],
								StructuralDof.TranslationY];
						if (neumannLoad.ContainsKey(dofIDY))
						{
							neumannLoad[dofIDY] += jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] *
							                       loadGaussPoint[1];
						}
						else
						{
							neumannLoad.Add(dofIDY,
								jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPoint[1]);
						}
					}


					if (element.Model.GlobalDofOrdering.GlobalFreeDofs.Contains(
						elementControlPoints[k],
						StructuralDof.TranslationZ))
					{
						var dofIDZ =
							element.Model.GlobalDofOrdering.GlobalFreeDofs[
								elementControlPoints[k],
								StructuralDof.TranslationZ];
						if (neumannLoad.ContainsKey(dofIDZ))
						{
							neumannLoad[dofIDZ] += jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] *
							                       loadGaussPoint[2];
						}
						else
						{
							neumannLoad.Add(dofIDZ,
								jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPoint[2]);
						}
					}


				}
			}
		}

		private static void CalculatePressure1D(
			Element element,
			Edge edge,
			PressureBoundaryCondition pressure,
			IList<ControlPoint> controlPoints,
			IList<GaussLegendrePoint3D> gaussPoints,
			IDictionary<int, double> pressureLoad)
		{
			var nurbs = new Nurbs1D(element, controlPoints, edge);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				double xGaussPoint = 0;
				double yGaussPoint = 0;
				double jacobian1 = 0.0;
				double jacobian2 = 0.0;
				var elementControlPoints = element.ControlPointsDictionary.Values.ToArray();
				for (int k = 0; k < elementControlPoints.Length; k++)
				{
					xGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].X;
					yGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].Y;
					jacobian1 += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
					jacobian2 += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				}

				double jacdet = Math.Sqrt(Math.Pow(jacobian1, 2) + Math.Pow(jacobian2, 2));

				double norm = Math.Sqrt(Math.Pow(xGaussPoint, 2) + Math.Pow(yGaussPoint, 2));
				var loadGaussPointX = pressure.Value * xGaussPoint / norm;
				var loadGaussPointY = pressure.Value * yGaussPoint / norm;

				for (int k = 0; k < elementControlPoints.Length; k++)
				{
					int dofIDX =
						element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[k], StructuralDof.TranslationX];
					int dofIDY =
						element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[k], StructuralDof.TranslationY];
					if (pressureLoad.ContainsKey(dofIDX))
					{
						pressureLoad[dofIDX] +=
							jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPointX;
					}
					else
					{
						pressureLoad.Add(dofIDX, jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPointX);
					}

					if (pressureLoad.ContainsKey(dofIDY))
					{
						pressureLoad[dofIDY] +=
							jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPointY;
					}
					else
					{
						pressureLoad.Add(dofIDY, jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPointY);
					}
				}
			}
		}

		private static int FindAxisControlPointId3D(Element element, Edge edge, ControlPoint controlPoint)
		{
			int ID = -1;
			switch (edge.ID)
			{
				case 1:
				case 2:
				case 5:
				case 6:
					ID = controlPoint.ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) %
						 element.Patch.NumberOfControlPointsZeta;
					break;

				case 3:
				case 4:
				case 7:
				case 8:
					ID = controlPoint.ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) /
						 element.Patch.NumberOfControlPointsZeta;
					break;

				case 9:
				case 10:
				case 11:
				case 12:
					ID = controlPoint.ID / (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta);
					break;
			}

			return ID;
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(((NurbsElement1D)element).Degree, element.Knots.ToList());
		}

		#endregion IStructuralIsogeometricElement
	}
}
