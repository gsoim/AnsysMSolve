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
	/// A three-dimensional continuum element that utilizes NURBS for shape functions.
	/// Authors: Dimitris Tsapetis.
	/// </summary>
	public class NurbsElement3D : Element, IStructuralIsogeometricElement
	{
		protected static readonly IDofType[] ControlPointDofTypes = { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private IDofType[][] _dofTypes;

		/// <summary>
		/// Retrieves the type of Finite Element used. Since the element is Isogeometric its type is defined as unknown.
		/// </summary>
		public CellType CellType { get; } = CellType.Unknown;

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		/// <summary>
		/// Retrieves the number of Dimensions of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		/// <summary>
		/// Boolean property that determines whether the material used for this elements has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Calculates the forces applies to an <see cref="NurbsElement3D"/> due to <see cref="MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => throw new NotImplementedException();

		/// <summary>
		/// Calculates displacements of knots for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
		{
			Contract.Requires(element != null, "The element cannot be null");

			var nurbsElement = (NurbsElement3D)element;
			var elementControlPoints = nurbsElement.ControlPoints.ToArray();
			var elementKnots = nurbsElement.Knots.ToArray();
			var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] { elementKnots[0].Ksi, elementKnots[4].Ksi });
			var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { elementKnots[0].Heta, elementKnots[2].Heta });
			var knotParametricCoordinatesΖeta = Vector.CreateFromArray(new double[] { elementKnots[0].Zeta, elementKnots[1].Zeta });
			Nurbs3D nurbs = new Nurbs3D(nurbsElement, elementControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta, knotParametricCoordinatesΖeta);
			var knotDisplacements = new double[8, 3];
			var paraviewKnotRenumbering = new int[] { 0, 4, 2, 6, 1, 5, 3, 7 };
			for (int j = 0; j < elementKnots.Length; j++)
			{
				for (int i = 0; i < elementControlPoints.Length; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] += nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
					knotDisplacements[paraviewKnotRenumbering[j], 2] += nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
				}
			}

			return knotDisplacements;
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsElement3D"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="neumann"><inheritdoc cref="NeumannBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann) => throw new NotSupportedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsElement3D"/> as it refers to two-dimensional loads.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement3D"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="NeumannBoundaryCondition"/> was applied to.</param>
		/// <param name="neumann">The <see cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann) => throw new NotSupportedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsElement3D"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure) => throw new NotSupportedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsElement3D"/> as it refers to two-dimensional loads.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement3D"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="PressureBoundaryCondition"/> was applied to.</param>
		/// <param name="pressure">The <see cref="PressureBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure) => throw new NotSupportedException();

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
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
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="NurbsElement3D"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="ControlPoint"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			Contract.Requires(element != null, "The element cannot be null");

			var nurbsElement = (NurbsElement3D)element;
			_dofTypes = new IDofType[nurbsElement.ControlPointsDictionary.Count][];
			for (int i = 0; i < nurbsElement.ControlPointsDictionary.Count; i++)
			{
				_dofTypes[i] = ControlPointDofTypes;
			}

			return _dofTypes;
		}

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="NurbsElement3D"/>.</returns>
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
		/// <param name="element">An element of type <see cref="NurbsElement3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="NurbsElement3D"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element)
		{
			var nurbsElement = (NurbsElement3D)element;
			var elementControlPoints = nurbsElement.ControlPoints.ToArray();
			IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(nurbsElement);
			Matrix stiffnessMatrixElement = Matrix.CreateZero(nurbsElement.ControlPointsDictionary.Count * 3, nurbsElement.ControlPointsDictionary.Count * 3);

			Nurbs3D nurbs = new Nurbs3D(nurbsElement, elementControlPoints);

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				Matrix jacobianMatrix = CalculateJacobian(elementControlPoints, nurbs, j);

				double jacdet = CalculateJacobianDeterminant(jacobianMatrix);

				Matrix inverseJacobian = CalculateInverseJacobian(jacobianMatrix, jacdet);

				Matrix B1 = CalculateDeformationMatrix1(inverseJacobian);

				Matrix B2 = CalculateDeformationMatrix2(elementControlPoints, nurbs, j);

				Matrix B = B1 * B2;

				IMatrixView E = ((IContinuumMaterial3D)nurbsElement.Patch.Material).ConstitutiveMatrix;
				Matrix stiffnessMatrixGaussPoint = B.ThisTransposeTimesOtherTimesThis(E) * jacdet * gaussPoints[j].WeightFactor;

				for (int m = 0; m < elementControlPoints.Length * 3; m++)
				{
					for (int n = 0; n < elementControlPoints.Length * 3; n++)
					{
						stiffnessMatrixElement[m, n] += stiffnessMatrixGaussPoint[m, n];
					}
				}
			}
			return stiffnessMatrixElement;
		}

		private static Matrix CalculateDeformationMatrix1(Matrix inverseJacobian)
		{
			Matrix B1 = Matrix.CreateZero(6, 9);

			B1[0, 0] += inverseJacobian[0, 0];
			B1[0, 1] += inverseJacobian[0, 1];
			B1[0, 2] += inverseJacobian[0, 2];

			B1[1, 3] += inverseJacobian[1, 0];
			B1[1, 4] += inverseJacobian[1, 1];
			B1[1, 5] += inverseJacobian[1, 2];

			B1[2, 6] += inverseJacobian[2, 0];
			B1[2, 7] += inverseJacobian[2, 1];
			B1[2, 8] += inverseJacobian[2, 2];

			B1[3, 0] += inverseJacobian[1, 0];
			B1[3, 1] += inverseJacobian[1, 1];
			B1[3, 2] += inverseJacobian[1, 2];
			B1[3, 3] += inverseJacobian[0, 0];
			B1[3, 4] += inverseJacobian[0, 1];
			B1[3, 5] += inverseJacobian[0, 2];

			B1[4, 3] += inverseJacobian[2, 0];
			B1[4, 4] += inverseJacobian[2, 1];
			B1[4, 5] += inverseJacobian[2, 2];
			B1[4, 6] += inverseJacobian[1, 0];
			B1[4, 7] += inverseJacobian[1, 1];
			B1[4, 8] += inverseJacobian[1, 2];

			B1[5, 0] += inverseJacobian[2, 0];
			B1[5, 1] += inverseJacobian[2, 1];
			B1[5, 2] += inverseJacobian[2, 2];
			B1[5, 6] += inverseJacobian[0, 0];
			B1[5, 7] += inverseJacobian[0, 1];
			B1[5, 8] += inverseJacobian[0, 2];
			return B1;
		}

		private static Matrix CalculateDeformationMatrix2(IList<ControlPoint> elementControlPoints, Nurbs3D nurbs, int j)
		{
			Matrix B2 = Matrix.CreateZero(9, 3 * elementControlPoints.Count);
			for (int column = 0; column < 3 * elementControlPoints.Count; column += 3)
			{
				B2[0, column] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[1, column] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[2, column] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];

				B2[3, column + 1] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[4, column + 1] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[5, column + 1] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];

				B2[6, column + 2] += nurbs.NurbsDerivativeValuesKsi[column / 3, j];
				B2[7, column + 2] += nurbs.NurbsDerivativeValuesHeta[column / 3, j];
				B2[8, column + 2] += nurbs.NurbsDerivativeValuesZeta[column / 3, j];
			}

			return B2;
		}

		private static Matrix CalculateInverseJacobian(Matrix jacobianMatrix, double jacdet)
		{
			Matrix inverseJacobian = Matrix.CreateZero(3, 3);

			inverseJacobian[0, 0] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
			inverseJacobian[0, 1] = jacobianMatrix[0, 2] * jacobianMatrix[2, 1] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2];
			inverseJacobian[0, 2] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1];

			inverseJacobian[1, 0] = jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2];
			inverseJacobian[1, 1] = jacobianMatrix[0, 0] * jacobianMatrix[2, 2] - jacobianMatrix[0, 2] * jacobianMatrix[2, 0];
			inverseJacobian[1, 2] = jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[0, 0] * jacobianMatrix[1, 2];

			inverseJacobian[2, 0] = jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[1, 1] * jacobianMatrix[2, 0];
			inverseJacobian[2, 1] = jacobianMatrix[0, 1] * jacobianMatrix[2, 0] - jacobianMatrix[0, 0] * jacobianMatrix[2, 1];
			inverseJacobian[2, 2] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

			inverseJacobian = inverseJacobian * (1 / jacdet);
			return inverseJacobian;
		}

		private static Matrix CalculateJacobian(IList<ControlPoint> elementControlPoints, Nurbs3D nurbs, int j)
		{
			Matrix jacobianMatrix = Matrix.CreateZero(3, 3);

			for (int k = 0; k < elementControlPoints.Count; k++)
			{
				jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[2, 0] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[2, 1] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[2, 2] += nurbs.NurbsDerivativeValuesZeta[k, j] * elementControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		private static double CalculateJacobianDeterminant(Matrix jacobianMatrix)
		{
			return jacobianMatrix[0, 0] * (jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[2, 1] * jacobianMatrix[1, 2])
								- jacobianMatrix[0, 1] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 2] - jacobianMatrix[2, 0] * jacobianMatrix[1, 2])
								+ jacobianMatrix[0, 2] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[2, 0] * jacobianMatrix[1, 1]);
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Patch.DegreeZeta, element.Knots.ToArray());
		}
	}
}
