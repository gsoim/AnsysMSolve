using System.Diagnostics.Contracts;

namespace MGroup.IGA.Elements
{
	using System;
	using System.Collections.Generic;
	using System.Linq;

	using MGroup.IGA.Entities;
	using MGroup.IGA.Entities.Loads;
	using MGroup.IGA.Interfaces;
	using MGroup.IGA.SupportiveClasses;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Materials.Interfaces;
	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.FreedomDegrees;
	using MGroup.MSolve.Discretization.Interfaces;
	using MGroup.MSolve.Discretization.Loads;
	using MGroup.MSolve.Discretization.Mesh;

	/// <summary>
	/// A two-dimensional continuum element that utilizes NURBS for shape functions.
	/// Authors: Dimitris Tsapetis.
	/// </summary>
	public class NurbsElement2D : Element, IStructuralIsogeometricElement
	{
		protected static readonly IDofType[] ControlPointDofTypes = { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private IDofType[][] _dofTypes;

		/// <summary>
		/// Retrieves the type of Finite Element used. Since the element is Isogeometric its type is defined as unknown.
		/// </summary>
		public CellType CellType { get; } = CellType.Unknown;

		#region IStructuralIsogeometricElement

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		/// <summary>
		/// Retrieves the number of Dimensions of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		/// <summary>
		/// Boolean property that determines whether the material used for this elements has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Calculates displacements of knots for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
		{
			Contract.Requires(element != null, "The element cannot be null");
			var nurbsElement = (NurbsElement2D)element;
			var elementControlPoints = nurbsElement.ControlPoints.ToArray();
			var elemenetKnots = nurbsElement.Knots.ToArray();
			var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] { elemenetKnots[0].Ksi, elemenetKnots[2].Ksi });
			var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { elemenetKnots[0].Heta, elemenetKnots[1].Heta });
			Nurbs2D nurbs = new Nurbs2D(nurbsElement, elementControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);
			var knotDisplacements = new double[4, 2];
			var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			for (int j = 0; j < elemenetKnots.Length; j++)
			{
				for (int i = 0; i < elementControlPoints.Length; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] += nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
				}
			}

			return knotDisplacements;
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsElement2D"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="neumann"><inheritdoc cref="NeumannBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann) => throw new NotSupportedException();

		/// <summary>
		/// This method calculates the Neumann boundary condition when applied to a two-dimensional NURBS element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/></param>
		/// <param name="face">An two dimensional boundary entity. For more info see <see cref="Face"/>.</param>
		/// <param name="neumann"><inheritdoc cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann)
		{
			Contract.Requires(element != null, "The element cannot be null");
			Contract.Requires(face != null, "The face cannot be null");
			Contract.Requires(neumann != null, "The Neumann Boundary condition cannot be null");

			IList<GaussLegendrePoint3D> gaussPoints =
				CreateElementGaussPoints(element, face.Degrees[0], face.Degrees[1]);
			Dictionary<int, double> neumannLoad = new Dictionary<int, double>();
			var elementControlPoints = element.ControlPoints.ToArray();
			Nurbs2D nurbs = new Nurbs2D(element, elementControlPoints, face);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = JacobianMatrixForLoadCalculation(element, nurbs, j, out var xGaussPoint, out var yGaussPoint, out var zGaussPoint);

				Vector surfaceBasisVector1 = Vector.CreateZero(3);
				surfaceBasisVector1[0] = jacobianMatrix[0, 0];
				surfaceBasisVector1[1] = jacobianMatrix[0, 1];
				surfaceBasisVector1[2] = jacobianMatrix[0, 2];

				Vector surfaceBasisVector2 = Vector.CreateZero(3);
				surfaceBasisVector2[0] = jacobianMatrix[1, 0];
				surfaceBasisVector2[1] = jacobianMatrix[1, 1];
				surfaceBasisVector2[2] = jacobianMatrix[1, 2];

				Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);

				double jacdet = jacobianMatrix[0, 0] * jacobianMatrix[1, 1]
								- jacobianMatrix[1, 0] * jacobianMatrix[0, 1];

				CalculateNeumannLoad2D(element, neumann, neumannLoad, nurbs, j, jacdet, gaussPoints, xGaussPoint, yGaussPoint, zGaussPoint, surfaceBasisVector3);
			}

			return neumannLoad;
		}

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState() => throw new NotImplementedException();

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="NurbsElement2D"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="ControlPoint"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			Contract.Requires(element != null, "The element cannot be null");
			var nurbsElement = (NurbsElement2D)element;
			_dofTypes = new IDofType[nurbsElement.ControlPointsDictionary.Count][];
			for (var i = 0; i < nurbsElement.ControlPointsDictionary.Count; i++)
			{
				_dofTypes[i] = ControlPointDofTypes;
			}

			return _dofTypes;
		}

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="NurbsElement2D"/>.</returns>
		public IMatrix MassMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="NurbsElement2D"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element)
		{
			Contract.Requires(element != null, "The element cannot be null");
			var nurbsElement = (NurbsElement2D)element;
			var gaussPoints = CreateElementGaussPoints(nurbsElement);
			var stiffnessMatrixElement = Matrix.CreateZero(
				nurbsElement.ControlPointsDictionary.Count * 2,
				nurbsElement.ControlPointsDictionary.Count * 2);
			var elementControlPoints = nurbsElement.ControlPoints.ToArray();
			var nurbs = new Nurbs2D(nurbsElement, elementControlPoints);

			for (var j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = Matrix.CreateZero(2, 2);

				for (int k = 0; k < elementControlPoints.Length; k++)
				{
					jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
					jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
					jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].X;
					jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				}

				double jacdet = (jacobianMatrix[0, 0] * jacobianMatrix[1, 1])
								- (jacobianMatrix[1, 0] * jacobianMatrix[0, 1]);

				Matrix B1 = Matrix.CreateZero(3, 4);

				B1[0, 0] += jacobianMatrix[1, 1] / jacdet;
				B1[0, 1] += -jacobianMatrix[0, 1] / jacdet;
				B1[1, 2] += -jacobianMatrix[1, 0] / jacdet;
				B1[1, 3] += jacobianMatrix[0, 0] / jacdet;
				B1[2, 0] += -jacobianMatrix[1, 0] / jacdet;
				B1[2, 1] += jacobianMatrix[0, 0] / jacdet;
				B1[2, 2] += jacobianMatrix[1, 1] / jacdet;
				B1[2, 3] += -jacobianMatrix[0, 1] / jacdet;

				Matrix B2 = Matrix.CreateZero(4, 2 * elementControlPoints.Length);
				for (int column = 0; column < 2 * elementControlPoints.Length; column += 2)
				{
					B2[0, column] += nurbs.NurbsDerivativeValuesKsi[column / 2, j];
					B2[1, column] += nurbs.NurbsDerivativeValuesHeta[column / 2, j];
					B2[2, column + 1] += nurbs.NurbsDerivativeValuesKsi[column / 2, j];
					B2[3, column + 1] += nurbs.NurbsDerivativeValuesHeta[column / 2, j];
				}

				Matrix B = B1 * B2;
				IMatrixView elasticityMatrix = ((IContinuumMaterial2D)nurbsElement.Patch.Material).ConstitutiveMatrix;
				Matrix stiffnessMatrixGaussPoint = B.ThisTransposeTimesOtherTimesThis(elasticityMatrix);
				stiffnessMatrixGaussPoint *= jacdet * gaussPoints[j].WeightFactor * nurbsElement.Patch.Thickness;

				for (int m = 0; m < elementControlPoints.Length * 2; m++)
				{
					for (int n = 0; n < elementControlPoints.Length * 2; n++)
					{
						stiffnessMatrixElement[m, n] += stiffnessMatrixGaussPoint[m, n];
					}
				}
			}

			return stiffnessMatrixElement;
		}

		private static void CalculateNeumannLoad2D(Element element, NeumannBoundaryCondition neumann, Dictionary<int, double> neumannLoad,
			Nurbs2D nurbs, int j, double jacdet, IList<GaussLegendrePoint3D> gaussPoints, double xGaussPoint, double yGaussPoint, double zGaussPoint,
			Vector surfaceBasisVector3)
		{
			var elementControlPoints = element.ControlPoints.ToArray();
			for (int k = 0; k < elementControlPoints.Length; k++)
			{
				int dofIDX =
					element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[k], StructuralDof.TranslationX];
				int dofIDY =
					element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[k], StructuralDof.TranslationY];
				int dofIDZ =
					element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[k], StructuralDof.TranslationZ];

				if (neumannLoad.ContainsKey(dofIDX))
				{
					neumannLoad[dofIDX] += nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
										   neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[0] *
										   surfaceBasisVector3[0];
				}
				else
				{
					neumannLoad.Add(
						dofIDX,
						nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
						neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[0] * surfaceBasisVector3[0]);
				}

				if (neumannLoad.ContainsKey(dofIDY))
				{
					neumannLoad[dofIDY] += nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
										   neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[1] *
										   surfaceBasisVector3[1];
				}
				else
				{
					neumannLoad.Add(
						dofIDY,
						nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor * neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[1] * surfaceBasisVector3[1]);
				}

				if (neumannLoad.ContainsKey(dofIDZ))
				{
					neumannLoad[dofIDZ] += nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor * neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[2] * surfaceBasisVector3[2];
				}
				else
				{
					neumannLoad.Add(dofIDZ,
						nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor * neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint)[2] * surfaceBasisVector3[2]);
				}
			}
		}

		private static Matrix JacobianMatrixForLoadCalculation(Element element, Nurbs2D nurbs, int j, out double xGaussPoint,
			out double yGaussPoint, out double zGaussPoint)
		{
			var elementControlPoints = element.ControlPoints.ToArray();
			Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
			xGaussPoint = 0;
			yGaussPoint = 0;
			zGaussPoint = 0;
			for (int k = 0; k < elementControlPoints.Length; k++)
			{
				xGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].X;
				yGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].Y;
				zGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		#endregion IStructuralIsogeometricElement

		/// <summary>
		/// Calculates the forces applies to an <see cref="NurbsElement2D"/> due to <see cref="MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => throw new NotImplementedException();

		/// <summary>
		/// This method calculates the Pressure boundary condition when applied to a two-dimensional NURBS element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="face">An two dimensional boundary entity. For more info see <see cref="Face"/>.</param>
		/// <param name="pressure"><inheritdoc cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure)
		{
			Contract.Requires(element != null, "The element cannot be null");
			Contract.Requires(face != null, "The face cannot be null");
			Contract.Requires(pressure != null, "The pressure boundary condition cannot be null");

			var dofs = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

			IList<GaussLegendrePoint3D> gaussPoints =
				CreateElementGaussPoints(element, face.Degrees[0], face.Degrees[1]);
			Dictionary<int, double> pressureLoad = new Dictionary<int, double>();
			var elementControlPoints = element.ControlPoints.ToArray();
			Nurbs2D nurbs = new Nurbs2D(element, elementControlPoints, face);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
				double xGaussPoint = 0;
				double yGaussPoint = 0;
				double zGaussPoint = 0;
				for (int k = 0; k < elementControlPoints.Length; k++)
				{
					xGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].X;
					yGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].Y;
					zGaussPoint += nurbs.NurbsValues[k, j] * elementControlPoints[k].Z;
					jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
					jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
					jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
					jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].X;
					jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
					jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
				}

				Vector surfaceBasisVector1 = Vector.CreateZero(3);
				surfaceBasisVector1[0] = jacobianMatrix[0, 0];
				surfaceBasisVector1[1] = jacobianMatrix[0, 1];
				surfaceBasisVector1[2] = jacobianMatrix[0, 2];

				Vector surfaceBasisVector2 = Vector.CreateZero(3);
				surfaceBasisVector2[0] = jacobianMatrix[1, 0];
				surfaceBasisVector2[1] = jacobianMatrix[1, 1];
				surfaceBasisVector2[2] = jacobianMatrix[1, 2];

				Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);

				double jacdet = (jacobianMatrix[0, 0] * jacobianMatrix[1, 1])
								- (jacobianMatrix[1, 0] * jacobianMatrix[0, 1]);

				for (int k = 0; k < elementControlPoints.Length; k++)
				{
					for (int m = 0; m < 3; m++)
					{
						int dofID = element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[k], dofs[m]];
						if (pressureLoad.ContainsKey(dofID))
						{
							pressureLoad[dofID] += nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor *
												   pressure.Value * surfaceBasisVector3[m];
						}
						else
						{
							pressureLoad.Add(dofID,
								nurbs.NurbsValues[k, j] * jacdet * gaussPoints[j].WeightFactor * pressure.Value * surfaceBasisVector3[m]);
						}
					}
				}
			}

			return pressureLoad;
		}

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsElement2D"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure) => throw new NotSupportedException();

		/// <summary>
		/// Clear any saved material states of the element.
		/// </summary>
		public void ClearMaterialStresses() => throw new NotImplementedException();

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified() => throw new NotImplementedException();

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState() => throw new NotImplementedException();

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Knots.ToArray());
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element, int degreeKsi, int degreeHeta)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, element.Knots.ToArray());
		}
	}
}
