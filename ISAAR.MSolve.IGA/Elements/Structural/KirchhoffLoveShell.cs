using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;
using Element = ISAAR.MSolve.IGA.Entities.Element;

namespace ISAAR.MSolve.IGA.Elements.Structural
{
  //  public enum GeometricalFormulation
  //  {
		//Linear, 
		//NonLinear
  //  }

  //  public class KirchhoffLoveShellFactory
  //  {
  //      public KirchhoffLoveShell CreateKirchhoffLoveShellElement(GeometricalFormulation formulation,
  //          IShellSectionMaterial material, IShapeFunction2D shapeFunctions,
  //          GaussLegendrePoint3D[] gaussPoints, double thickness)
  //      {

		//	//var geometricalCalculation= formulation==GeometricalFormulation.Linear? new LinearKLFormulation():new NonLinearKLFormulation;
		//	return new KirchhoffLoveShell(material,shapeFunctions,gaussPoints,thickness);
  //      }


  //  }


    /// <summary>
	/// An shell element that utilizes Non-Uniform Rational B-Splines for shape functions.
	/// It is based on Kirchhoff-Love theory. Geometrically linear formulation.
	/// For more information please refer to <see href="https://www.sciencedirect.com/science/article/pii/S0045782509002680"/>
	/// Authors: Dimitris Tsapetis.
	/// </summary>
	public class KirchhoffLoveShell : Element, IStructuralIsogeometricElement, ISurfaceLoadedElement
	{
        public KirchhoffLoveShell(IShellSectionMaterial material,
            IShapeFunction2D shapeFunctions, GaussLegendrePoint3D[] gaussPoints, 
            double thickness )
        {
            _material = material;
            _shapeFunctions = shapeFunctions;
            _gaussPoints = gaussPoints;
            _thickness = thickness;

            foreach (var gaussPoint in _gaussPoints)
                materialsAtMidsurfaceGP.Add(gaussPoint, _material.Clone());
        }

        private readonly Dictionary<GaussLegendrePoint3D, IShellSectionMaterial> materialsAtMidsurfaceGP =
            new Dictionary<GaussLegendrePoint3D, IShellSectionMaterial>();
		protected static readonly IDofType[] ControlPointDofTypes = { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        private readonly IShellSectionMaterial _material;
        internal readonly IShapeFunction2D _shapeFunctions;
        private readonly GaussLegendrePoint3D[] _gaussPoints;
        private readonly double _thickness;
        private IDofType[][] dofTypes;

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
		public bool MaterialModified => throw new NotImplementedException();

		/// <summary>
		/// Calculates the forces applies to an <see cref="NurbsKirchhoffLoveShellElement"/> due to <see cref="FEM.Entities.MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <param name="loads">A list of <see cref="FEM.Entities.MassAccelerationLoad"/>. For more info see <seealso cref="FEM.Entities.MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => throw new NotImplementedException();

		/// <summary>
		/// Calculates displacements of knots for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
        {
            throw new NotImplementedException();
			//var nurbsElement = (NurbsKirchhoffLoveShellElement)element;
			//var elementControlPoints = nurbsElement.ControlPoints.ToArray();
			//var elementKnots = nurbsElement.Knots.ToArray();
			//var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] { elementKnots[0].Ksi, elementKnots[2].Ksi });
			//var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { elementKnots[0].Heta, elementKnots[1].Heta });
			//var knotDisplacements = new double[4, 3];
			//var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			//for (int j = 0; j < elementKnots.Length; j++)
			//{
			//	for (int i = 0; i < elementControlPoints.Length; i++)
			//	{
			//		knotDisplacements[paraviewKnotRenumbering[j], 0] += nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
			//		knotDisplacements[paraviewKnotRenumbering[j], 1] += nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
			//		knotDisplacements[paraviewKnotRenumbering[j], 2] += nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
			//	}
			//}

			//return knotDisplacements;
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsKirchhoffLoveShellElement"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="neumann"><inheritdoc cref="NeumannBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsKirchhoffLoveShellElement"/> as it refers to two-dimensional loads.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="NeumannBoundaryCondition"/> was applied to.</param>
		/// <param name="neumann">The <see cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsKirchhoffLoveShellElement"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="NurbsKirchhoffLoveShellElement"/> as it refers to two-dimensional loads.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="PressureBoundaryCondition"/> was applied to.</param>
		/// <param name="pressure">The <see cref="PressureBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure) => throw new NotImplementedException();

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
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
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of a <see cref="NurbsKirchhoffLoveShellElement"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="ControlPoint"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			dofTypes = new IDofType[element.Nodes.Count][];
			for (int i = 0; i < element.Nodes.Count; i++)
			{
				dofTypes[i] = ControlPointDofTypes;
			}

			return dofTypes;
		}

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="NurbsKirchhoffLoveShellElement"/>.</returns>
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
		/// <param name="element">An element of type <see cref="NurbsKirchhoffLoveShellElement"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="NurbsKirchhoffLoveShellElement"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element)
        {
            var numberOfCP = element.Nodes.Count;

			Matrix stiffnessMatrixElement = Matrix.CreateZero(numberOfCP * 3,
                numberOfCP * 3);


			for (int j = 0; j < _gaussPoints.Length; j++)
			{
				var jacobianMatrix = CalculateJacobian(element, _shapeFunctions, j);

				var hessianMatrix = CalculateHessian(element, _shapeFunctions, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
				var J1 = surfaceBasisVector3.Norm2();
				surfaceBasisVector3.ScaleIntoThis(1 / J1);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                var material = materialsAtMidsurfaceGP[_gaussPoints[j]];
                material.TangentVectorV1 = surfaceBasisVector1.CopyToArray();
                material.TangentVectorV2 = surfaceBasisVector2.CopyToArray();
                material.NormalVectorV3 = surfaceBasisVector3.CopyToArray();
                material.Thickness = this._thickness;

                var membraneConstitutiveMatrix = material.MembraneConstitutiveMatrix;
                var bendingConstitutiveMatrix = material.BendingConstitutiveMatrix;

                var Bmembrane = CalculateMembraneDeformationMatrix(_shapeFunctions, j, surfaceBasisVector1,
                    surfaceBasisVector2, numberOfCP);

                var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, _shapeFunctions, j,
                    surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1,
                    surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, numberOfCP);
				
				var Kmembrane = Bmembrane.Transpose() * membraneConstitutiveMatrix * Bmembrane *  J1 *
								_gaussPoints[j].WeightFactor;
				
				var Kbending = Bbending.Transpose() * bendingConstitutiveMatrix * Bbending  * J1 *
							   _gaussPoints[j].WeightFactor;

				stiffnessMatrixElement.AddIntoThis(Kmembrane);
				stiffnessMatrixElement.AddIntoThis(Kbending);
			}
			return stiffnessMatrixElement;
		}

		private static Matrix CalculateHessian(IElement shellElement, IShapeFunction2D nurbs, int j)
		{
			var elementControlPoints = shellElement.Nodes.ToArray();
			Matrix hessianMatrix = Matrix.CreateZero(3, 3);
			for (int k = 0; k < elementControlPoints.Length; k++)
			{
				hessianMatrix[0, 0] += nurbs.SecondDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
				hessianMatrix[0, 1] += nurbs.SecondDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				hessianMatrix[0, 2] += nurbs.SecondDerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
				hessianMatrix[1, 0] += nurbs.SecondDerivativeValuesHeta[k, j] * elementControlPoints[k].X;
				hessianMatrix[1, 1] += nurbs.SecondDerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				hessianMatrix[1, 2] += nurbs.SecondDerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
				hessianMatrix[2, 0] += nurbs.SecondDerivativeValuesKsiHeta[k, j] * elementControlPoints[k].X;
				hessianMatrix[2, 1] += nurbs.SecondDerivativeValuesKsiHeta[k, j] * elementControlPoints[k].Y;
				hessianMatrix[2, 2] += nurbs.SecondDerivativeValuesKsiHeta[k, j] * elementControlPoints[k].Z;
			}

			return hessianMatrix;
		}

		private static Matrix CalculateJacobian(IElement shellElement, IShapeFunction2D nurbs, int j)
		{
			var elementControlPoints = shellElement.Nodes.ToArray();
			Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
			for (int k = 0; k < elementControlPoints.Length; k++)
			{
				jacobianMatrix[0, 0] += nurbs.DerivativeValuesKsi[k, j] * elementControlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.DerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[0, 2] += nurbs.DerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[1, 0] += nurbs.DerivativeValuesHeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.DerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[1, 2] += nurbs.DerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		private static Vector CalculateSurfaceBasisVector1(Matrix Matrix, int row)
		{
			Vector surfaceBasisVector1 = Vector.CreateZero(3);
			surfaceBasisVector1[0] = Matrix[row, 0];
			surfaceBasisVector1[1] = Matrix[row, 1];
			surfaceBasisVector1[2] = Matrix[row, 2];
			return surfaceBasisVector1;
		}

		private Matrix CalculateBendingDeformationMatrix(Vector surfaceBasisVector3, IShapeFunction2D nurbs, int j,
			Vector surfaceBasisVector2, Vector surfaceBasisVectorDerivative1, Vector surfaceBasisVector1, double J1,
			Vector surfaceBasisVectorDerivative2, Vector surfaceBasisVectorDerivative12, int numberOfControlPoints)
		{
			Matrix Bbending = Matrix.CreateZero(3, numberOfControlPoints * 3);
			for (int column = 0; column < numberOfControlPoints * 3; column += 3)
			{
				#region BI1

				var BI1 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
				BI1.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				var auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI1.AddIntoThis(auxVector);
				
				BI1.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative1));
				auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative1);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				BI1.AddIntoThis(auxVector);

				auxVector = surfaceBasisVectorDerivative1.CrossProduct(surfaceBasisVector2);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI1.AddIntoThis(auxVector);

				BI1.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-nurbs.SecondDerivativeValuesKsi[column / 3, j]);
				BI1.AddIntoThis(auxVector);

				#endregion BI1

				#region BI2

				Vector BI2 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
				BI2.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				BI2.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative2));
				auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative2);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				auxVector = surfaceBasisVectorDerivative2.CrossProduct(surfaceBasisVector2);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				BI2.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-nurbs.SecondDerivativeValuesHeta[column / 3, j]);
				BI2.AddIntoThis(auxVector);

				#endregion BI2

				#region BI3

				Vector BI3 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
				BI3.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				BI3.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative12));
				auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative12);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				auxVector = surfaceBasisVectorDerivative2.CrossProduct(surfaceBasisVector2);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				BI3.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-nurbs.SecondDerivativeValuesKsiHeta[column / 3, j]);
				BI3.AddIntoThis(auxVector);

				#endregion BI3

				Bbending[0, column] = BI1[0];
				Bbending[0, column + 1] = BI1[1];
				Bbending[0, column + 2] = BI1[2];

				Bbending[1, column] = BI2[0];
				Bbending[1, column + 1] = BI2[1];
				Bbending[1, column + 2] = BI2[2];

				Bbending[2, column] = 2 * BI3[0];
				Bbending[2, column + 1] = 2 * BI3[1];
				Bbending[2, column + 2] = 2 * BI3[2];
			}

			return Bbending;
		}

		private Matrix CalculateConstitutiveMatrix(KirchhoffLoveShell element, Vector surfaceBasisVector1, Vector surfaceBasisVector2)
		{
			var auxMatrix1 = Matrix.CreateZero(2, 2);
			auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
			auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
			auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
			auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
			(Matrix inverse, double det) = auxMatrix1.InvertAndDeterminant();

			var material = _material;
			var constitutiveMatrix = Matrix.CreateFromArray(new double[3, 3]
			{
				{
					inverse[0,0]*inverse[0,0],
					material.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-material.PoissonRatio)*inverse[1,0]*inverse[1,0],
					inverse[0,0]*inverse[1,0]
				},
				{
					material.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-material.PoissonRatio)*inverse[1,0]*inverse[1,0],
					inverse[1,1]*inverse[1,1],
					inverse[1,1]*inverse[1,0]
				},
				{
					inverse[0,0]*inverse[1,0],
					inverse[1,1]*inverse[1,0],
					0.5*(1-material.PoissonRatio)*inverse[0,0]*inverse[1,1]+(1+material.PoissonRatio)*inverse[1,0]*inverse[1,0]
				},
			});
			return constitutiveMatrix;
		}

		private Matrix CalculateMembraneDeformationMatrix(IShapeFunction2D nurbs, int j, Vector surfaceBasisVector1,
			Vector surfaceBasisVector2, int numberOfControlPoints)
		{
			Matrix dRIa = Matrix.CreateZero(3, numberOfControlPoints);
			for (int i = 0; i < numberOfControlPoints; i++)
			{
				for (int m = 0; m < 3; m++)
				{
					dRIa[m, i] = nurbs.DerivativeValuesHeta[i, j] * surfaceBasisVector1[m] +
								 nurbs.DerivativeValuesKsi[i, j] * surfaceBasisVector2[m];
				}
			}

			Matrix Bmembrane = Matrix.CreateZero(3, numberOfControlPoints * 3);
			for (int column = 0; column < numberOfControlPoints * 3; column += 3)
			{
				Bmembrane[0, column] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[0];
				Bmembrane[0, column + 1] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[1];
				Bmembrane[0, column + 2] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[2];

				Bmembrane[1, column] = nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[0];
				Bmembrane[1, column + 1] = nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[1];
				Bmembrane[1, column + 2] = nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[2];

				Bmembrane[2, column] = dRIa[0, column / 3];
				Bmembrane[2, column + 1] = dRIa[1, column / 3];
				Bmembrane[2, column + 2] = dRIa[2, column / 3];
			}

			return Bmembrane;
		}

		public Dictionary<int, double> CalculateSurfacePressure(Element element, double pressureMagnitude)
		{
			var shellElement = (KirchhoffLoveShell)element;
			var elementControlPoints = shellElement.ControlPoints.ToArray();
			var pressureLoad = new Dictionary<int, double>();

			for (var j = 0; j < _gaussPoints.Length; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, _shapeFunctions, j);
				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
				var surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
				var J1 = surfaceBasisVector3.Norm2();
				surfaceBasisVector3.ScaleIntoThis(1 / J1);

				for (int i = 0; i < elementControlPoints.Length; i++)
				{
					for (int k = 0; k < ControlPointDofTypes.Length; k++)
					{
						int dofId = element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[i], ControlPointDofTypes[k]];

						if (pressureLoad.ContainsKey(dofId))
						{
							pressureLoad[dofId] += pressureMagnitude * surfaceBasisVector3[k] *
                                                   _shapeFunctions.Values[i, j] * _gaussPoints[j].WeightFactor;
						}
						else
						{
							pressureLoad.Add(dofId, pressureMagnitude * surfaceBasisVector3[k] * _shapeFunctions.Values[i, j] * _gaussPoints[j].WeightFactor);
						}
					}
				}
			}

			return pressureLoad;
		}

		public Dictionary<int, double> CalculateSurfaceDistributedLoad(Element element, IDofType loadedDof, double loadMagnitude)
		{
			var shellElement = (KirchhoffLoveShell)element;
			var elementControlPoints = shellElement.ControlPoints.ToArray();
			var distributedLoad = new Dictionary<int, double>();

			for (var j = 0; j < _gaussPoints.Length; j++)
			{
				var jacobianMatrix = CalculateJacobian(shellElement, _shapeFunctions, j);
				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
				var surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
				var J1 = surfaceBasisVector3.Norm2();
				surfaceBasisVector3.ScaleIntoThis(1 / J1);

				for (int i = 0; i < elementControlPoints.Length; i++)
				{
					var loadedDofIndex = ControlPointDofTypes.FindFirstIndex(loadedDof);
					if (!element.Model.GlobalDofOrdering.GlobalFreeDofs.Contains(elementControlPoints[i], loadedDof))
						continue;
					var dofId = element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[i], loadedDof];

					if (distributedLoad.ContainsKey(dofId))
					{
						distributedLoad[dofId] += loadMagnitude * J1 *
                                                  _shapeFunctions.Values[i, j] * _gaussPoints[j].WeightFactor;
					}
					else
					{
						distributedLoad.Add(dofId, loadMagnitude * _shapeFunctions.Values[i, j] * J1 * _gaussPoints[j].WeightFactor);
					}
				}
			}

			return distributedLoad;
		}
	}
}
