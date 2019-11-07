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
	/// An shell element that utilizes T-Splines for shape functions.
	/// It is based on Kirchhoff-Love theory. Geometrically linear formulation.
	/// The material constitutive laws are integrated through the thickness at each midsurface integration point to take into account any material non-linearities.
	/// Authors: Dimitris Tsapetis.
	/// </summary>
	public class TSplineKirchhoffLoveShellElementMaterial : Element, IStructuralIsogeometricElement
	{
		protected static readonly IDofType[] ControlPointDofTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private IDofType[][] _dofTypes;
		private const int ThicknessIntegrationDegree = 2;

		private readonly Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>
			_materialsAtThicknessGp =
				new Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>();

		private readonly Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>> _thicknessIntegrationPoints = new Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>>();

		/// <summary>
		/// Creates a <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.
		/// </summary>
		/// <param name="id">The element id.</param>
		/// <param name="patch">The patch that contains the element.</param>
		/// <param name="degreeKsi">Degree for parametric axis Ksi.</param>
		/// <param name="degreeHeta">Degree for parametric axis Heta.</param>
		/// <param name="thickness">Shell thickness. Constant throughout the shell.</param>
		/// <param name="extractionOperator">Bezier extraction operation from TSplines to Bezier elements.</param>
		/// <param name="shellMaterial">Material of the shell element.</param>
		public TSplineKirchhoffLoveShellElementMaterial(int id, Patch patch, int degreeKsi, int degreeHeta,
			double thickness, Matrix extractionOperator, IShellMaterial shellMaterial)
		{
			this.ID = id;
			this.Patch = patch;
			this.DegreeKsi = degreeKsi;
			this.DegreeHeta = degreeHeta;
			this.Thickness = thickness;
			this.ExtractionOperator = extractionOperator;

			CreateElementGaussPoints(this);
			foreach (var medianSurfaceGP in _thicknessIntegrationPoints.Keys)
			{
				_materialsAtThicknessGp.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IShellMaterial>());
				foreach (var point in _thicknessIntegrationPoints[medianSurfaceGP])
				{
					_materialsAtThicknessGp[medianSurfaceGP].Add(point, shellMaterial.Clone());
				}
			}
		}

		/// <summary>
		/// Retrieves the type of Finite Element used. Since the element is Isogeometric its type is defined as unknown.
		/// </summary>
		public CellType CellType { get; } = CellType.Unknown;

		/// <summary>
		/// Property that Polynomial degree of the T-Splines shape functions per axis Heta.
		/// </summary>
		public int DegreeHeta { get; set; }

		/// <summary>
		/// Property that Polynomial degree of the T-Splines shape functions per axis Ksi.
		/// </summary>
		public int DegreeKsi { get; set; }

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
		/// A <see cref="Matrix"/> that contains the Bezier extraction operator.
		/// </summary>
		public Matrix ExtractionOperator { get; set; }

		/// <summary>
		/// Boolean property that determines whether the material used for this elements has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Shell thickness. Constant throughout the element.
		/// </summary>
		public double Thickness { get; set; }

		/// <summary>
		/// Calculates the forces applies to an <see cref="TSplineKirchhoffLoveShellElementMaterial"/> due to <see cref="MassAccelerationLoad"/>
		/// </summary>
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/></param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => throw new NotImplementedException();

		/// <summary>
		/// Calculates displacements of knots for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
		{
			var tsplineElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			var elementControlPoints = tsplineElement.ControlPoints.ToArray();
			var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] { -1, 1 });
			var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { -1, 1 });

			var tsplines = new ShapeTSplines2DFromBezierExtraction(tsplineElement, elementControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

			var knotDisplacements = new double[4, 3];
			var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			for (var j = 0; j < knotDisplacements.GetLength(0); j++)
			{
				for (int i = 0; i < elementControlPoints.Length; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += tsplines.TSplineValues[i, j] * localDisplacements[i, 0];
					knotDisplacements[paraviewKnotRenumbering[j], 1] += tsplines.TSplineValues[i, j] * localDisplacements[i, 1];
					knotDisplacements[paraviewKnotRenumbering[j], 2] += tsplines.TSplineValues[i, j] * localDisplacements[i, 2];
				}
			}

			return knotDisplacements;
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			var gaussPoints = CreateElementGaussPoints(shellElement);
			var ElementNodalForces = new double[shellElement.ControlPointsDictionary.Count * 3];
			var tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints.ToArray());
			var elementControlPoints = shellElement.ControlPoints.ToArray();

			for (var j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(elementControlPoints, tsplines, j);

				var hessianMatrix = CalculateHessian(elementControlPoints, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
				};

				var norm = surfaceBasisVector3.Sum(t => t * t);
				var J1 = Math.Sqrt(norm);

				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1,
					surfaceBasisVector2, elementControlPoints);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, elementControlPoints);

				var (MembraneForces, BendingMoments) =
					IntegratedStressesOverThickness(gaussPoints, j);

				var wfactor = J1 * gaussPoints[j].WeightFactor;
				for (var i = 0; i < Bmembrane.GetLength(1); i++)
				{
					for (var k = 0; k < Bmembrane.GetLength(0); k++)
					{
						ElementNodalForces[i] += Bmembrane[k, i] * MembraneForces[k] * wfactor +
												 Bbending[k, i] * BendingMoments[k] * wfactor;
					}
				}
			}

			return ElementNodalForces;
		}

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="TSplineKirchhoffLoveShellElementMaterial"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="neumann"><inheritdoc cref="NeumannBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="TSplineKirchhoffLoveShellElementMaterial"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="NeumannBoundaryCondition"/> was applied to.</param>
		/// <param name="neumann">The <see cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="TSplineKirchhoffLoveShellElementMaterial"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure) => throw new NotImplementedException();

		/// <summary>
		/// This method cannot be used, combined with <see cref="TSplineKirchhoffLoveShellElementMaterial"/> as it refers to one-dimensional loads.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="PressureBoundaryCondition"/> was applied to.</param>
		/// <param name="pressure">The <see cref="PressureBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure) => throw new NotImplementedException();

		/// <summary>
		/// Calculates knots physical coordinates for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <returns>A <see cref="double"/> array calculating the coordinates of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		public double[,] CalculatePointsForPostProcessing(TSplineKirchhoffLoveShellElementMaterial element)
		{
			var localCoordinates = new double[4, 2]
			{
				{-1, -1},
				{-1, 1},
				{1, -1},
				{1, 1}
			};

			var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] { -1, 1 });
			var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { -1, 1 });
			var elementControlPoints = element.ControlPoints.ToArray();
			var tsplines = new ShapeTSplines2DFromBezierExtraction(element, elementControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);

			var knotDisplacements = new double[4, 3];
			var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			for (int j = 0; j < localCoordinates.GetLength(0); j++)
			{
				for (int i = 0; i < elementControlPoints.Length; i++)
				{
					knotDisplacements[paraviewKnotRenumbering[j], 0] += tsplines.TSplineValues[i, j] * elementControlPoints[i].X;
					knotDisplacements[paraviewKnotRenumbering[j], 1] += tsplines.TSplineValues[i, j] * elementControlPoints[i].Y;
					knotDisplacements[paraviewKnotRenumbering[j], 2] += tsplines.TSplineValues[i, j] * elementControlPoints[i].Z;
				}
			}

			return knotDisplacements;
		}

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			var elementControlPoints = shellElement.ControlPoints.ToArray();
			var tsplines = new ShapeTSplines2DFromBezierExtraction(shellElement, shellElement.ControlPoints.ToArray());

			for (var j = 0; j < _materialsAtThicknessGp.Keys.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(elementControlPoints, tsplines, j);

				var hessianMatrix = CalculateHessian(elementControlPoints, tsplines, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
				};

				var norm = surfaceBasisVector3.Sum(t => t * t);

				var J1 = Math.Sqrt(norm);

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1, surfaceBasisVector2, elementControlPoints);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, elementControlPoints);

				var membraneStrain = new double[Bmembrane.GetLength(0)];
				var bendingStrain = new double[Bmembrane.GetLength(0)];
				for (var i = 0; i < Bmembrane.GetLength(1); i++)
				{
					for (var k = 0; k < Bmembrane.GetLength(0); k++)
					{
						membraneStrain[i] += Bmembrane[k, i] * localDisplacements[k];
						bendingStrain[i] += Bbending[k, i] * localDisplacements[k];
					}
				}

				foreach (var keyValuePair in _materialsAtThicknessGp[_materialsAtThicknessGp.Keys.ToList()[j]])
				{
					var thicknessPoint = keyValuePair.Key;
					var material = keyValuePair.Value;
					var gpStrain = new double[bendingStrain.Length];
					var z = -thicknessPoint.Zeta;
					for (var i = 0; i < bendingStrain.Length; i++)
					{
						gpStrain[i] += membraneStrain[i] + bendingStrain[i] * z;
					}
					material.UpdateMaterial(gpStrain);
				}
			}
			return new Tuple<double[], double[]>(new double[0], new double[0]);
		}

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
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of a <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="ControlPoint"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			var nurbsElement = (TSplineKirchhoffLoveShellElementMaterial)element;
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
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</returns>
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
		/// <param name="element">An element of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element)
		{
			var shellElement = (TSplineKirchhoffLoveShellElementMaterial)element;
			var gaussPoints = _materialsAtThicknessGp.Keys.ToArray();
			var elementControlPoints = shellElement.ControlPoints.ToArray();
			ShapeTSplines2DFromBezierExtraction tsplines =
				new ShapeTSplines2DFromBezierExtraction(shellElement, elementControlPoints);

			var bRows = 3;
			var bCols = elementControlPoints.Length * 3;
			var stiffnessMatrix = new double[bCols, bCols];
			var BmTranspose = new double[bCols, bRows];
			var BbTranspose = new double[bCols, bRows];

			var BmTransposeMultStiffness = new double[bCols, bRows];
			var BbTransposeMultStiffness = new double[bCols, bRows];
			var BmbTransposeMultStiffness = new double[bCols, bRows];
			var BbmTransposeMultStiffness = new double[bCols, bRows];

			for (int j = 0; j < gaussPoints.Length; j++)
			{
				var jacobianMatrix = CalculateJacobian(elementControlPoints, tsplines, j);

				var hessianMatrix = CalculateHessian(elementControlPoints, tsplines, j);
				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
				};

				foreach (var integrationPointMaterial in _materialsAtThicknessGp[gaussPoints[j]].Values)
				{
					integrationPointMaterial.TangentVectorV1 = surfaceBasisVector1;
					integrationPointMaterial.TangentVectorV2 = surfaceBasisVector2;
					integrationPointMaterial.NormalVectorV3 = surfaceBasisVector3;
				}

				double norm = 0;
				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					norm += surfaceBasisVector3[i] * surfaceBasisVector3[i];
				var J1 = Math.Sqrt(norm);

				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(tsplines, j, surfaceBasisVector1,
					surfaceBasisVector2, elementControlPoints);
				var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, tsplines, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12, elementControlPoints);

				var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
					IntegratedConstitutiveOverThickness(gaussPoints, j);

				double wFactor = J1 * gaussPoints[j].WeightFactor;
				double tempb = 0;
				double tempm = 0;
				Array.Clear(BmTranspose, 0, bRows * bCols);
				Array.Clear(BbTranspose, 0, bRows * bCols);
				for (int i = 0; i < bRows; i++)
				{
					for (int k = 0; k < bCols; k++)
					{
						BmTranspose[k, i] = Bmembrane[i, k] * wFactor;
						BbTranspose[k, i] = Bbending[i, k] * wFactor;
					}
				}

				double tempcm = 0;
				double tempcb = 0;
				double tempcc = 0;
				Array.Clear(BmTransposeMultStiffness, 0, bRows * bCols);
				Array.Clear(BbTransposeMultStiffness, 0, bRows * bCols);
				Array.Clear(BmbTransposeMultStiffness, 0, bRows * bCols);
				Array.Clear(BbmTransposeMultStiffness, 0, bRows * bCols);
				for (int i = 0; i < bCols; i++)
				{
					for (int k = 0; k < bRows; k++)
					{
						tempm = BmTranspose[i, k];
						tempb = BbTranspose[i, k];
						for (int m = 0; m < bRows; m++)
						{
							tempcm = MembraneConstitutiveMatrix[k, m];
							tempcb = BendingConstitutiveMatrix[k, m];
							tempcc = CouplingConstitutiveMatrix[k, m];

							BmTransposeMultStiffness[i, m] += tempm * tempcm;
							BbTransposeMultStiffness[i, m] += tempb * tempcb;
							BmbTransposeMultStiffness[i, m] += tempm * tempcc;
							BbmTransposeMultStiffness[i, m] += tempb * tempcc;
						}
					}
				}

				double tempmb = 0;
				double tempbm = 0;
				double mem = 0;
				double ben = 0;
				for (int i = 0; i < bCols; i++)
				{
					for (int k = 0; k < bRows; k++)
					{
						tempm = BmTransposeMultStiffness[i, k];
						tempb = BbTransposeMultStiffness[i, k];
						tempmb = BmbTransposeMultStiffness[i, k];
						tempbm = BbmTransposeMultStiffness[i, k];

						for (int m = 0; m < bCols; m++)
						{
							mem = Bmembrane[k, m];
							ben = Bbending[k, m];
							stiffnessMatrix[i, m] += tempm * mem + tempb * ben + tempmb * ben + tempbm * mem;
						}
					}
				}
			}
			return Matrix.CreateFromArray(stiffnessMatrix);
		}

		private static double[,] CalculateHessian(ControlPoint[] elementControlPoints, ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			double[,] hessianMatrix = new double[3, 3];
			for (int k = 0; k < elementControlPoints.Length; k++)
			{
				hessianMatrix[0, 0] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * elementControlPoints[k].X;
				hessianMatrix[0, 1] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * elementControlPoints[k].Y;
				hessianMatrix[0, 2] += tsplines.TSplineSecondDerivativesValueKsi[k, j] * elementControlPoints[k].Z;
				hessianMatrix[1, 0] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * elementControlPoints[k].X;
				hessianMatrix[1, 1] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * elementControlPoints[k].Y;
				hessianMatrix[1, 2] += tsplines.TSplineSecondDerivativesValueHeta[k, j] * elementControlPoints[k].Z;
				hessianMatrix[2, 0] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * elementControlPoints[k].X;
				hessianMatrix[2, 1] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * elementControlPoints[k].Y;
				hessianMatrix[2, 2] += tsplines.TSplineSecondDerivativesValueKsiHeta[k, j] * elementControlPoints[k].Z;
			}

			return hessianMatrix;
		}

		private static double[,] CalculateJacobian(ControlPoint[] elementControlPoints, ShapeTSplines2DFromBezierExtraction tsplines, int j)
		{
			var jacobianMatrix = new double[2, 3];
			for (var k = 0; k < elementControlPoints.Length; k++)
			{
				jacobianMatrix[0, 0] += tsplines.TSplineDerivativeValuesKsi[k, j] * elementControlPoints[k].X;
				jacobianMatrix[0, 1] += tsplines.TSplineDerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[0, 2] += tsplines.TSplineDerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[1, 0] += tsplines.TSplineDerivativeValuesHeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[1, 1] += tsplines.TSplineDerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[1, 2] += tsplines.TSplineDerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		private static double[] CalculateSurfaceBasisVector1(double[,] Matrix, int row)
		{
			var surfaceBasisVector1 = new double[3];
			surfaceBasisVector1[0] = Matrix[row, 0];
			surfaceBasisVector1[1] = Matrix[row, 1];
			surfaceBasisVector1[2] = Matrix[row, 2];
			return surfaceBasisVector1;
		}

		private double[,] CalculateBendingDeformationMatrix(double[] surfaceBasisVector3, ShapeTSplines2DFromBezierExtraction tsplines, int j,
									double[] surfaceBasisVector2, double[] surfaceBasisVectorDerivative1, double[] surfaceBasisVector1, double J1,
			double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12, ControlPoint[] elementControlPoints)
		{
			var Bbending = new double[3, elementControlPoints.Length * 3];
			var s1 = Vector.CreateFromArray(surfaceBasisVector1);
			var s2 = Vector.CreateFromArray(surfaceBasisVector2);
			var s3 = Vector.CreateFromArray(surfaceBasisVector3);
			var s11 = Vector.CreateFromArray(surfaceBasisVectorDerivative1);
			var s22 = Vector.CreateFromArray(surfaceBasisVectorDerivative2);
			var s12 = Vector.CreateFromArray(surfaceBasisVectorDerivative12);
			for (int column = 0; column < elementControlPoints.Length * 3; column += 3)
			{
				#region BI1

				var BI1 = s3.CrossProduct(s3);
				BI1.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				var auxVector = s2.CrossProduct(s3);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI1.AddIntoThis(auxVector);
				BI1.ScaleIntoThis(s3.DotProduct(s11));
				auxVector = s1.CrossProduct(s11);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI1.AddIntoThis(auxVector);
				BI1.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueKsi[column / 3, j]);
				BI1.AddIntoThis(auxVector);

				#endregion BI1

				#region BI2

				IVector BI2 = s3.CrossProduct(s3);
				BI2.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				auxVector = s2.CrossProduct(s3);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				BI2.ScaleIntoThis(s3.DotProduct(s22));
				auxVector = s1.CrossProduct(s22);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				auxVector = s22.CrossProduct(s2);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				BI2.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueHeta[column / 3, j]);
				BI2.AddIntoThis(auxVector);

				#endregion BI2

				#region BI3

				Vector BI3 = s3.CrossProduct(s3);
				BI3.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				auxVector = s2.CrossProduct(s3);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				BI3.ScaleIntoThis(s3.DotProduct(s12));
				auxVector = s1.CrossProduct(s12);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesHeta[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				auxVector = s22.CrossProduct(s2);
				auxVector.ScaleIntoThis(tsplines.TSplineDerivativeValuesKsi[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				BI3.ScaleIntoThis(1 / J1);
				auxVector[0] = surfaceBasisVector3[0];
				auxVector[1] = surfaceBasisVector3[1];
				auxVector[2] = surfaceBasisVector3[2];
				auxVector.ScaleIntoThis(-tsplines.TSplineSecondDerivativesValueKsiHeta[column / 3, j]);
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

		private double[,] CalculateMembraneDeformationMatrix(ShapeTSplines2DFromBezierExtraction tsplines, int j, double[] surfaceBasisVector1,
			double[] surfaceBasisVector2, ControlPoint[] elementControlPoints)
		{
			var dRIa = new double[3, elementControlPoints.Length * 3];
			for (var i = 0; i < elementControlPoints.Length; i++)
			{
				for (var m = 0; m < 3; m++)
				{
					dRIa[m, i] = tsplines.TSplineDerivativeValuesHeta[i, j] * surfaceBasisVector1[m] +
								 tsplines.TSplineDerivativeValuesKsi[i, j] * surfaceBasisVector2[m];
				}
			}

			var Bmembrane = new double[3, elementControlPoints.Length * 3];
			for (var column = 0; column < elementControlPoints.Length * 3; column += 3)
			{
				Bmembrane[0, column] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[0];
				Bmembrane[0, column + 1] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[1];
				Bmembrane[0, column + 2] = tsplines.TSplineDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[2];

				Bmembrane[1, column] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[0];
				Bmembrane[1, column + 1] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[1];
				Bmembrane[1, column + 2] = tsplines.TSplineDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[2];

				Bmembrane[2, column] = dRIa[0, column / 3];
				Bmembrane[2, column + 1] = dRIa[1, column / 3];
				Bmembrane[2, column + 2] = dRIa[2, column / 3];
			}

			return Bmembrane;
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(TSplineKirchhoffLoveShellElementMaterial element)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			var medianSurfaceGp = gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta, new List<Knot>
			{
				new Knot() {ID = 0, Ksi = -1, Heta = -1, Zeta = 0},
				new Knot() {ID = 1, Ksi = -1, Heta = 1, Zeta = 0},
				new Knot() {ID = 2, Ksi = 1, Heta = -1, Zeta = 0},
				new Knot() {ID = 3, Ksi = 1, Heta = 1, Zeta = 0},
			});

			foreach (var point in medianSurfaceGp)
			{
				_thicknessIntegrationPoints.Add(point, gauss.CalculateElementGaussPoints(ThicknessIntegrationDegree,
					new List<Knot>
					{
						new Knot() {ID = 0, Ksi = -Thickness / 2, Heta = point.Heta },
						new Knot() {ID = 1, Ksi = Thickness / 2, Heta = point.Heta},
					}).Select(gp => new GaussLegendrePoint3D(point.Ksi, point.Heta, gp.Ksi, gp.WeightFactor)).ToList());
			}

			return medianSurfaceGp;
		}

		private (double[,] MembraneConstitutiveMatrix, double[,] BendingConstitutiveMatrix, double[,] CouplingConstitutiveMatrix) IntegratedConstitutiveOverThickness(IList<GaussLegendrePoint3D> gaussPoints, int j)
		{
			var MembraneConstitutiveMatrix = new double[3, 3];
			var BendingConstitutiveMatrix = new double[3, 3];
			var CouplingConstitutiveMatrix = new double[3, 3];

			foreach (var keyValuePair in _materialsAtThicknessGp[gaussPoints[j]])
			{
				var thicknessPoint = keyValuePair.Key;
				var material = keyValuePair.Value;
				var constitutiveMatrixM = material.ConstitutiveMatrix;
				double tempc = 0;
				double w = thicknessPoint.WeightFactor;
				double z = thicknessPoint.Zeta;
				for (int i = 0; i < 3; i++)
				{
					for (int k = 0; k < 3; k++)
					{
						tempc = constitutiveMatrixM[i, k];
						MembraneConstitutiveMatrix[i, k] += tempc * w;
						CouplingConstitutiveMatrix[i, k] += tempc * w * z;
						BendingConstitutiveMatrix[i, k] += tempc * w * z * z;
					}
				}
			}

			return (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix);
		}

		private (double[] MembraneForces, double[] BendingMoments) IntegratedStressesOverThickness(
			IList<GaussLegendrePoint3D> gaussPoints, int j)
		{
			var MembraneForces = new double[3];
			var BendingMoments = new double[3];

			foreach (var keyValuePair in _materialsAtThicknessGp[gaussPoints[j]])
			{
				var thicknessPoint = keyValuePair.Key;
				var material = keyValuePair.Value;
				var w = thicknessPoint.WeightFactor;
				var z = thicknessPoint.Zeta;
				for (int i = 0; i < 3; i++)
				{
					MembraneForces[i] += material.Stresses[i] * w * Thickness / 2;
					BendingMoments[i] += material.Stresses[i] * w * z * z * Thickness / 2;
				}
			}

			return (MembraneForces, BendingMoments);
		}
	}
}
