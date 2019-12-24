using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;
using Element = ISAAR.MSolve.IGA.Entities.Element;

namespace ISAAR.MSolve.IGA.Elements
{
    public class NurbsKirchhoffLoveShellElementSectionNL : Element, IStructuralIsogeometricElement, ISurfaceLoadedElement
	{
		protected static readonly IDofType[] ControlPointDofTypes = { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private IDofType[][] dofTypes;
		
		internal Dictionary<GaussLegendrePoint3D, IShellSectionMaterial>
			materialsAtMidsurfaceGP = new Dictionary<GaussLegendrePoint3D, IShellSectionMaterial>();

		private bool isInitialized;
		internal double[] _solution;

		public NurbsKirchhoffLoveShellElementSectionNL(IShellSectionMaterial shellMaterial, 
            IList<Knot> elementKnots, IList<ControlPoint> elementControlPoints, Nurbs2D nurbs,
            Patch patch, double thickness, int degreeKsi, int degreeHeta)
		{
			Contract.Requires(shellMaterial != null);
			this.Patch = patch;
			this.Thickness = thickness;
            _nurbs = nurbs;
            _degreeKsi = degreeKsi;
            _degreeHeta = degreeHeta;
			foreach (var knot in elementKnots)
			{
				if (!KnotsDictionary.ContainsKey(knot.ID))
					this.KnotsDictionary.Add(knot.ID, knot);
			}

			_solution = new double[3 * elementControlPoints.Count];

			foreach (var controlPoint in elementControlPoints)
			{
				if (!ControlPointsDictionary.ContainsKey(controlPoint.ID))
					ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
			}

			_midsurfaceGaussPoints=CreateElementGaussPoints(this);
			foreach (var medianSurfaceGP in _midsurfaceGaussPoints)
			{
				materialsAtMidsurfaceGP.Add(medianSurfaceGP, shellMaterial.Clone());
			}
		}

		public CellType CellType { get; } = CellType.Unknown;

		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public bool MaterialModified => false;

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => throw new NotImplementedException();

		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
		{
            throw new NotImplementedException();
			//var nurbsElement = (NurbsKirchhoffLoveShellElementNL)element;
			//var knotParametricCoordinatesKsi = Vector.CreateFromArray(Knots.Select(k => k.Ksi).ToArray());
			//var knotParametricCoordinatesHeta = Vector.CreateFromArray(Knots.Select(k => k.Heta).ToArray());

			//var knotDisplacements = new double[4, 3];
			//var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			//for (var j = 0; j < knotDisplacements.GetLength(0); j++)
			//{
			//	for (int i = 0; i < element.ControlPoints.Count(); i++)
			//	{
			//		knotDisplacements[paraviewKnotRenumbering[j], 0] +=
			//			nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
			//		knotDisplacements[paraviewKnotRenumbering[j], 1] +=
			//			nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
			//		knotDisplacements[paraviewKnotRenumbering[j], 2] +=
			//			nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
			//	}
			//}

			//return knotDisplacements;
		}

		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			var shellElement = (NurbsKirchhoffLoveShellElementSectionNL)element;
			var controlPoints = shellElement.ControlPoints.ToArray();
			var elementNodalForces = new double[shellElement.ControlPointsDictionary.Count * 3];

            var elementMembraneForces = new double[shellElement.ControlPointsDictionary.Count * 3];
            var elementBendingForces = new double[shellElement.ControlPointsDictionary.Count * 3];

            _solution = localDisplacements;
			//var newControlPoints = CurrentControlPoint(controlPoints);
			//var newControlPoints = controlPoints;
			var newControlPoints = CurrentControlPoint(controlPoints);

			var gaussPoints = _midsurfaceGaussPoints.ToArray();

			for (int j = 0; j < gaussPoints.Length; j++)
			{
				var jacobianMatrix = CalculateJacobian(newControlPoints, _nurbs, j);

				var hessianMatrix = CalculateHessian(newControlPoints, _nurbs, j);

				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
				};

				double norm = surfaceBasisVector3.Sum(t => t * t);

				var J1 = Math.Sqrt(norm);

				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(newControlPoints, _nurbs, j, surfaceBasisVector1,
					surfaceBasisVector2);
				var Bbending = CalculateBendingDeformationMatrix(newControlPoints, surfaceBasisVector3, _nurbs, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12);
				var material = materialsAtMidsurfaceGP[gaussPoints[j]];
				var membraneForces = material.MembraneForces;
				var bendingMoments = material.Moments;
				
				var wfactor = InitialJ1[j] * gaussPoints[j].WeightFactor;
				for (int i = 0; i < Bmembrane.GetLength(1); i++)
				{
					for (int k = 0; k < Bmembrane.GetLength(0); k++)
                    {
                        elementMembraneForces[i] += Bmembrane[k, i] * membraneForces[k] * wfactor;
                        elementBendingForces[i] += Bbending[k, i] * bendingMoments[k] * wfactor;
                        elementNodalForces[i] += (+Bmembrane[k, i] * membraneForces[k] * wfactor +
												 +Bbending[k, i] * bendingMoments[k] * wfactor);
					}
				}
			}

			return elementNodalForces;
		}

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann) => throw new NotImplementedException();

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann) => throw new NotImplementedException();

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure) => throw new NotImplementedException();

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure) => throw new NotImplementedException();

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			var shellElement = (NurbsKirchhoffLoveShellElementSectionNL)element;
			var elementControlPoints = shellElement.ControlPoints.ToArray();

			_solution = localDisplacements;

			var newControlPoints = CurrentControlPoint(elementControlPoints);
			//var newControlPoints = elementControlPoints;

			var midsurfaceGP = _midsurfaceGaussPoints.ToArray();
			for (var j = 0; j < midsurfaceGP.Length; j++)
			{
				var jacobianMatrix = CalculateJacobian(newControlPoints, _nurbs, j);

				var hessianMatrix = CalculateHessian(newControlPoints, _nurbs, j);

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

				var unitVector3 = new double[]
				{
					surfaceBasisVector3[0] / J1, surfaceBasisVector3[1] / J1, surfaceBasisVector3[2] / J1
				};

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var A11 = initialSurfaceBasisVectors1[j][0] * initialSurfaceBasisVectors1[j][0] +
						  initialSurfaceBasisVectors1[j][1] * initialSurfaceBasisVectors1[j][1] +
						  initialSurfaceBasisVectors1[j][2] * initialSurfaceBasisVectors1[j][2];

				var A22 = initialSurfaceBasisVectors2[j][0] * initialSurfaceBasisVectors2[j][0] +
						 initialSurfaceBasisVectors2[j][1] * initialSurfaceBasisVectors2[j][1] +
						 initialSurfaceBasisVectors2[j][2] * initialSurfaceBasisVectors2[j][2];

				var A12 = initialSurfaceBasisVectors1[j][0] * initialSurfaceBasisVectors2[j][0] +
						  initialSurfaceBasisVectors1[j][1] * initialSurfaceBasisVectors2[j][1] +
						  initialSurfaceBasisVectors1[j][2] * initialSurfaceBasisVectors2[j][2];

				var a11 = surfaceBasisVector1[0] * surfaceBasisVector1[0] +
						  surfaceBasisVector1[1] * surfaceBasisVector1[1] +
						  surfaceBasisVector1[2] * surfaceBasisVector1[2];

				var a22 = surfaceBasisVector2[0] * surfaceBasisVector2[0] +
						  surfaceBasisVector2[1] * surfaceBasisVector2[1] +
						  surfaceBasisVector2[2] * surfaceBasisVector2[2];

				var a12 = surfaceBasisVector1[0] * surfaceBasisVector2[0] +
						  surfaceBasisVector1[1] * surfaceBasisVector2[1] +
						  surfaceBasisVector1[2] * surfaceBasisVector2[2];

				var membraneStrain = new double[] { 0.5 * (a11 - A11), 0.5 * (a22 - A22), a12 - A12 };

				var B11 = initialSurfaceBasisVectorDerivative1[j][0] * initialUnitSurfaceBasisVectors3[j][0] +
						  initialSurfaceBasisVectorDerivative1[j][1] * initialUnitSurfaceBasisVectors3[j][1] +
						  initialSurfaceBasisVectorDerivative1[j][2] * initialUnitSurfaceBasisVectors3[j][2];

				var B22 = initialSurfaceBasisVectorDerivative2[j][0] * initialUnitSurfaceBasisVectors3[j][0] +
						  initialSurfaceBasisVectorDerivative2[j][1] * initialUnitSurfaceBasisVectors3[j][1] +
						  initialSurfaceBasisVectorDerivative2[j][2] * initialUnitSurfaceBasisVectors3[j][2];

				var B12 = initialSurfaceBasisVectorDerivative12[j][0] * initialUnitSurfaceBasisVectors3[j][0] +
						  initialSurfaceBasisVectorDerivative12[j][1] * initialUnitSurfaceBasisVectors3[j][1] +
						  initialSurfaceBasisVectorDerivative12[j][2] * initialUnitSurfaceBasisVectors3[j][2];

				var b11 = surfaceBasisVectorDerivative1[0] * unitVector3[0] +
						  surfaceBasisVectorDerivative1[1] * unitVector3[1] +
						  surfaceBasisVectorDerivative1[2] * unitVector3[2];

				var b22 = surfaceBasisVectorDerivative2[0] * unitVector3[0] +
						  surfaceBasisVectorDerivative2[1] * unitVector3[1] +
						  surfaceBasisVectorDerivative2[2] * unitVector3[2];

				var b12 = surfaceBasisVectorDerivative12[0] * unitVector3[0] +
						 surfaceBasisVectorDerivative12[1] * unitVector3[1] +
						 surfaceBasisVectorDerivative12[2] * unitVector3[2];

				var bendingStrain = new double[] { (b11 - B11), (b22 - B22), (2 * b12 - 2 * B12) };

				materialsAtMidsurfaceGP[_midsurfaceGaussPoints[j]].UpdateMaterial(membraneStrain, bendingStrain);
			}

			return new Tuple<double[], double[]>(new double[0], new double[0]);
		}

		public Dictionary<int, double> CalculateSurfaceDistributedLoad(Element element, IDofType loadedDof, double loadMagnitude)
		{
			var shellElement = (NurbsKirchhoffLoveShellElementSectionNL)element;
			var elementControlPoints = shellElement.ControlPoints.ToArray();
			var gaussPoints = CreateElementGaussPoints(shellElement);
			var distributedLoad = new Dictionary<int, double>();

			for (var j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(elementControlPoints, _nurbs, j);
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
                                                  _nurbs.Values[i, j] * gaussPoints[j].WeightFactor;
					}
					else
					{
						distributedLoad.Add(dofId, loadMagnitude * _nurbs.Values[i, j] * J1 * gaussPoints[j].WeightFactor);
					}
				}
			}

			return distributedLoad;
		}

		public Dictionary<int, double> CalculateSurfacePressure(Element element, double pressureMagnitude)
		{
			var shellElement = (NurbsKirchhoffLoveShellElementSectionNL)element;
			var elementControlPoints = shellElement.ControlPoints.ToArray();
			var gaussPoints = CreateElementGaussPoints(shellElement);
			var pressureLoad = new Dictionary<int, double>();

			for (var j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(elementControlPoints, _nurbs, j);
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
                                                   _nurbs.Values[i, j] * gaussPoints[j].WeightFactor;
						}
						else
						{
							pressureLoad.Add(dofId, pressureMagnitude * surfaceBasisVector3[k] * _nurbs.Values[i, j] * gaussPoints[j].WeightFactor);
						}
					}
				}
			}

			return pressureLoad;
		}

		public void ClearMaterialState()
		{
		}

		public void ClearMaterialStresses() => throw new NotImplementedException();

		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			dofTypes = new IDofType[element.Nodes.Count][];
			for (var i = 0; i < element.Nodes.Count; i++)
			{
				dofTypes[i] = ControlPointDofTypes;
			}

			return dofTypes;
		}


		public IMatrix MassMatrix(IElement element) => throw new NotImplementedException();

		public void ResetMaterialModified() => throw new NotImplementedException();

		public void SaveMaterialState()
		{
			foreach (var material in materialsAtMidsurfaceGP.Values)
			{
				material.SaveState();
			}
		}

		public IMatrix StiffnessMatrix(IElement element)
		{
			var shellElement = (NurbsKirchhoffLoveShellElementSectionNL)element;
			var gaussPoints = _midsurfaceGaussPoints;

			var controlPoints = shellElement.ControlPoints.ToArray();
            
			if (!isInitialized)
			{
				CalculateInitialConfigurationData(controlPoints, _nurbs, gaussPoints);
				isInitialized = true;
			}

			var elementControlPoints = CurrentControlPoint(controlPoints);

            var KmemTotal = Matrix.CreateZero(elementControlPoints.Length * 3, elementControlPoints.Length * 3);
            var KmemTotalL = Matrix.CreateZero(elementControlPoints.Length * 3, elementControlPoints.Length * 3);
            var KmemTotalNL = Matrix.CreateZero(elementControlPoints.Length * 3, elementControlPoints.Length * 3);

            var KbenTotal = Matrix.CreateZero(elementControlPoints.Length * 3, elementControlPoints.Length * 3);
            var KbenTotalL = Matrix.CreateZero(elementControlPoints.Length * 3, elementControlPoints.Length * 3);
            var KbenTotalNL = Matrix.CreateZero(elementControlPoints.Length * 3, elementControlPoints.Length * 3);


            var bRows = 3;
			var bCols = elementControlPoints.Length * 3;
			var stiffnessMatrix = new double[bCols, bCols];
			var BmTranspose = new double[bCols, bRows];
			var BbTranspose = new double[bCols, bRows];

			var BmTransposeMultStiffness = new double[bCols, bRows];
			var BbTransposeMultStiffness = new double[bCols, bRows];
			var BmbTransposeMultStiffness = new double[bCols, bRows];
			var BbmTransposeMultStiffness = new double[bCols, bRows];

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(elementControlPoints, _nurbs, j);

				var hessianMatrix = CalculateHessian(elementControlPoints, _nurbs, j);
				var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

				var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

				var surfaceBasisVector3 = new[]
				{
					surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
					surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
					surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
				};

				double norm = 0;
				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					norm += surfaceBasisVector3[i] * surfaceBasisVector3[i];
				var J1 = Math.Sqrt(norm);

				for (int i = 0; i < surfaceBasisVector3.Length; i++)
					surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

				var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				var Bmembrane = CalculateMembraneDeformationMatrix(elementControlPoints, _nurbs, j, surfaceBasisVector1,
					surfaceBasisVector2);
				var Bbending = CalculateBendingDeformationMatrix(elementControlPoints, surfaceBasisVector3, _nurbs, j, surfaceBasisVector2,
					surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
					surfaceBasisVectorDerivative12);

				var material = materialsAtMidsurfaceGP[gaussPoints[j]];
				var MembraneConstitutiveMatrix = material.MembraneConstitutiveMatrix;
				var BendingConstitutiveMatrix = material.BendingConstitutiveMatrix;
				var CouplingConstitutiveMatrix = material.CouplingConstitutiveMatrix;

				double wFactor = InitialJ1[j] * gaussPoints[j].WeightFactor;
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

                var MembraneForces = materialsAtMidsurfaceGP[gaussPoints[j]].MembraneForces;
                var BendingMoments = materialsAtMidsurfaceGP[gaussPoints[j]].Moments;

                var KmembraneNL = CalculateKmembraneNL(elementControlPoints, MembraneForces, _nurbs, j);

                var KbendingNL = CalculateKbendingNL(elementControlPoints, BendingMoments, _nurbs,
                    Vector.CreateFromArray(surfaceBasisVector1), Vector.CreateFromArray(surfaceBasisVector2),
                    Vector.CreateFromArray(surfaceBasisVector3),
                    Vector.CreateFromArray(surfaceBasisVectorDerivative1),
                    Vector.CreateFromArray(surfaceBasisVectorDerivative2),
                    Vector.CreateFromArray(surfaceBasisVectorDerivative12), J1, j);

                KmemTotal.AddIntoThis(((KmembraneNL +
                                        Matrix.CreateFromArray(Bmembrane).Transpose() *
                                        MembraneConstitutiveMatrix *
                                        Matrix.CreateFromArray(Bmembrane)).Scale(wFactor)));
                KmemTotalL.AddIntoThis((Matrix.CreateFromArray(Bmembrane).Transpose() *
                                        MembraneConstitutiveMatrix *
                                        Matrix.CreateFromArray(Bmembrane)).Scale(wFactor));
                KmemTotalNL.AddIntoThis(KmembraneNL.Scale(wFactor));


                KbenTotal.AddIntoThis(((KbendingNL +
                                        Matrix.CreateFromArray(Bbending).Transpose() *
                                        BendingConstitutiveMatrix *
                                        Matrix.CreateFromArray(Bbending)).Scale(wFactor)));
                KbenTotalL.AddIntoThis((Matrix.CreateFromArray(Bbending).Transpose() *
                                        BendingConstitutiveMatrix *
                                        Matrix.CreateFromArray(Bbending)).Scale(wFactor));
                KbenTotalNL.AddIntoThis(KbendingNL.Scale(wFactor));


                for (var i = 0; i < stiffnessMatrix.GetLength(0); i++)
                {
                    for (var k = 0; k < stiffnessMatrix.GetLength(1); k++)
                    {
                        stiffnessMatrix[i, k] += KmembraneNL[i, k] * wFactor;
                        stiffnessMatrix[i, k] += KbendingNL[i, k] * wFactor;
                    }
                }
            }

			return Matrix.CreateFromArray(stiffnessMatrix);
		}

		internal ControlPoint[] CurrentControlPoint(ControlPoint[] controlPoints)
		{
			var cp = new ControlPoint[controlPoints.Length];

			for (int i = 0; i < controlPoints.Length; i++)
			{
				cp[i] = new ControlPoint()
				{
					X = controlPoints[i].X + _solution[i * 3],
					Y = controlPoints[i].Y + _solution[i * 3 + 1],
					Z = controlPoints[i].Z + _solution[i * 3 + 2],
					Ksi = controlPoints[i].Ksi,
					Heta = controlPoints[i].Heta,
					Zeta = controlPoints[i].Zeta,
					WeightFactor = controlPoints[i].WeightFactor
				};
			}

			return cp;
		}

		private ControlPoint[] DControlPoint(ControlPoint[] controlPoints)
		{
			var cp = new ControlPoint[controlPoints.Length];

			for (int i = 0; i < controlPoints.Length; i++)
			{
				cp[i] = new ControlPoint()
				{
					X = _solution[i * 3],
					Y = _solution[i * 3 + 1],
					Z = _solution[i * 3 + 2],
					Ksi = controlPoints[i].Ksi,
					Heta = controlPoints[i].Heta,
					Zeta = controlPoints[i].Zeta,
					WeightFactor = controlPoints[i].WeightFactor
				};
			}

			return cp;
		}

		internal double[,] CalculateHessian(ControlPoint[] controlPoints, Nurbs2D nurbs, int j)
		{
			var hessianMatrix = new double[3, 3];
			for (var k = 0; k < controlPoints.Length; k++)
			{
				hessianMatrix[0, 0] +=
					nurbs.SecondDerivativeValuesKsi[k, j] * controlPoints[k].X;
				hessianMatrix[0, 1] +=
					nurbs.SecondDerivativeValuesKsi[k, j] * controlPoints[k].Y;
				hessianMatrix[0, 2] +=
					nurbs.SecondDerivativeValuesKsi[k, j] * controlPoints[k].Z;
				hessianMatrix[1, 0] +=
					nurbs.SecondDerivativeValuesHeta[k, j] * controlPoints[k].X;
				hessianMatrix[1, 1] +=
					nurbs.SecondDerivativeValuesHeta[k, j] * controlPoints[k].Y;
				hessianMatrix[1, 2] +=
					nurbs.SecondDerivativeValuesHeta[k, j] * controlPoints[k].Z;
				hessianMatrix[2, 0] +=
					nurbs.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].X;
				hessianMatrix[2, 1] +=
					nurbs.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].Y;
				hessianMatrix[2, 2] +=
					nurbs.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].Z;
			}

			return hessianMatrix;
		}

		internal double[,] CalculateJacobian(ControlPoint[] controlPoints, Nurbs2D nurbs, int j)
		{
			var jacobianMatrix = new double[2, 3];
			for (var k = 0; k < controlPoints.Length; k++)
			{
				jacobianMatrix[0, 0] += nurbs.DerivativeValuesKsi[k, j] * controlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.DerivativeValuesKsi[k, j] * controlPoints[k].Y;
				jacobianMatrix[0, 2] += nurbs.DerivativeValuesKsi[k, j] * controlPoints[k].Z;
				jacobianMatrix[1, 0] += nurbs.DerivativeValuesHeta[k, j] * controlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.DerivativeValuesHeta[k, j] * controlPoints[k].Y;
				jacobianMatrix[1, 2] += nurbs.DerivativeValuesHeta[k, j] * controlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		internal double[] CalculateSurfaceBasisVector1(double[,] Matrix, int row)
		{
			var surfaceBasisVector1 = new double[3];
			surfaceBasisVector1[0] = Matrix[row, 0];
			surfaceBasisVector1[1] = Matrix[row, 1];
			surfaceBasisVector1[2] = Matrix[row, 2];
			return surfaceBasisVector1;
		}

		

		internal double[,] CalculateBendingDeformationMatrix(ControlPoint[] controlPoints, double[] surfaceBasisVector3,
			Nurbs2D nurbs, int j, double[] surfaceBasisVector2, double[] surfaceBasisVectorDerivative1, double[] surfaceBasisVector1,
			double J1, double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12)
		{
			var Bbending = new double[3, controlPoints.Length * 3];
			var s1 = Vector.CreateFromArray(surfaceBasisVector1);
			var s2 = Vector.CreateFromArray(surfaceBasisVector2);
			var s3 = Vector.CreateFromArray(surfaceBasisVector3);
			var s11 = Vector.CreateFromArray(surfaceBasisVectorDerivative1);
			var s22 = Vector.CreateFromArray(surfaceBasisVectorDerivative2);
			var s12 = Vector.CreateFromArray(surfaceBasisVectorDerivative12);
			for (int column = 0; column < controlPoints.Length * 3; column += 3)
			{
				#region BI1

				var BI1 = s3.CrossProduct(s1);
				BI1.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				var auxVector = s2.CrossProduct(s3);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI1.AddIntoThis(auxVector);

				BI1.ScaleIntoThis(s3.DotProduct(s11));
				auxVector = s1.CrossProduct(s11);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				BI1.AddIntoThis(auxVector);

				auxVector = s11.CrossProduct(s2);
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

				IVector BI2 = s3.CrossProduct(s1);
				BI2.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				auxVector = s2.CrossProduct(s3);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				BI2.ScaleIntoThis(s3.DotProduct(s22));
				auxVector = s1.CrossProduct(s22);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				BI2.AddIntoThis(auxVector);
				auxVector = s22.CrossProduct(s2);
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

				var BI3 = s3.CrossProduct(s1);
				BI3.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				auxVector = s2.CrossProduct(s3);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesKsi[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				BI3.ScaleIntoThis(s3.DotProduct(s12));
				auxVector = s1.CrossProduct(s12);
				auxVector.ScaleIntoThis(nurbs.DerivativeValuesHeta[column / 3, j]);
				BI3.AddIntoThis(auxVector);
				auxVector = s22.CrossProduct(s2);
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

		private double[] CalculateCrossProduct(double[] vector1, double[] vector2)
		{
			return new[]
			{
				vector1[1] * vector2[2] - vector1[2] * vector2[1],
				vector1[2] * vector2[0] - vector1[0] * vector2[2],
				vector1[0] * vector2[1] - vector1[1] * vector2[0]
			};
		}

		private double[][] initialSurfaceBasisVectors1;
		private double[][] initialSurfaceBasisVectors2;
		private double[][] initialUnitSurfaceBasisVectors3;

		private double[][] initialSurfaceBasisVectorDerivative1;
		private double[][] initialSurfaceBasisVectorDerivative2;
		private double[][] initialSurfaceBasisVectorDerivative12;

		private double[] InitialJ1;
		private IList<GaussLegendrePoint3D> _midsurfaceGaussPoints;

		internal void CalculateInitialConfigurationData(ControlPoint[] controlPoints,
			Nurbs2D nurbs, IList<GaussLegendrePoint3D> gaussPoints)
		{
			var numberOfGP = gaussPoints.Count;
			InitialJ1=new double[numberOfGP];
			initialSurfaceBasisVectors1 = new double[numberOfGP][];
			initialSurfaceBasisVectors2 = new double[numberOfGP][];
			initialUnitSurfaceBasisVectors3 = new double[numberOfGP][];
			initialSurfaceBasisVectorDerivative1 = new double[numberOfGP][];
			initialSurfaceBasisVectorDerivative2 = new double[numberOfGP][];
			initialSurfaceBasisVectorDerivative12 = new double[numberOfGP][];

			for (int j = 0; j < gaussPoints.Count; j++)
			{
				var jacobianMatrix = CalculateJacobian(controlPoints, nurbs, j);

				var hessianMatrix = CalculateHessian(controlPoints, nurbs, j);
				initialSurfaceBasisVectors1[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
				initialSurfaceBasisVectors2[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
				var s3= CalculateCrossProduct(initialSurfaceBasisVectors1[j], initialSurfaceBasisVectors2[j]);
				var norm = s3.Sum(t => t * t);
				InitialJ1[j] = Math.Sqrt(norm);
				 var vector3= CalculateCrossProduct(initialSurfaceBasisVectors1[j], initialSurfaceBasisVectors2[j]);
				 initialUnitSurfaceBasisVectors3[j]= new double[]
				 {
					 vector3[0]/InitialJ1[j],
					 vector3[1]/InitialJ1[j],
					 vector3[2]/InitialJ1[j],
				 };

				initialSurfaceBasisVectorDerivative1[j] = CalculateSurfaceBasisVector1(hessianMatrix, 0);
				initialSurfaceBasisVectorDerivative2[j] = CalculateSurfaceBasisVector1(hessianMatrix, 1);
				initialSurfaceBasisVectorDerivative12[j] = CalculateSurfaceBasisVector1(hessianMatrix, 2);

				materialsAtMidsurfaceGP[gaussPoints[j]].TangentVectorV1= initialSurfaceBasisVectors1[j];
				materialsAtMidsurfaceGP[gaussPoints[j]].TangentVectorV2 = initialSurfaceBasisVectors2[j];
				materialsAtMidsurfaceGP[gaussPoints[j]].NormalVectorV3 = initialUnitSurfaceBasisVectors3[j];
			}
		}

		internal Matrix CalculateKbendingNL(ControlPoint[] controlPoints,
		   double[] bendingMoments, Nurbs2D nurbs, Vector surfaceBasisVector1,
		   Vector surfaceBasisVector2, Vector surfaceBasisVector3, Vector surfaceBasisVectorDerivative1, Vector surfaceBasisVectorDerivative2,
		   Vector surfaceBasisVectorDerivative12, double J1, int j)
		{
			var KbendingNL =
				Matrix.CreateZero(controlPoints.Length * 3, controlPoints.Length * 3);

            for (int i = 0; i < controlPoints.Length; i++)
			{
				var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.DerivativeValuesKsi[i, j]);
				var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.DerivativeValuesHeta[i, j]);

				var a11r = Matrix3by3.CreateIdentity().Scale(nurbs.SecondDerivativeValuesKsi[i, j]);
				var a22r = Matrix3by3.CreateIdentity().Scale(nurbs.SecondDerivativeValuesHeta[i, j]);
				var a12r = Matrix3by3.CreateIdentity().Scale(nurbs.SecondDerivativeValuesKsiHeta[i, j]);
				for (int k = 0; k < controlPoints.Length; k++)
				{
                    var a1s = Matrix3by3.CreateIdentity().Scale(nurbs.DerivativeValuesKsi[k, j]);
                    var a2s = Matrix3by3.CreateIdentity().Scale(nurbs.DerivativeValuesHeta[k, j]);
                    
                    var a11s = Matrix3by3.CreateIdentity().Scale(nurbs.SecondDerivativeValuesKsi[k, j]);
					var a22s = Matrix3by3.CreateIdentity().Scale(nurbs.SecondDerivativeValuesHeta[k, j]);
					var a12s = Matrix3by3.CreateIdentity().Scale(nurbs.SecondDerivativeValuesKsiHeta[k, j]);

                    var a3r = CalculateA3r(nurbs, i, j, surfaceBasisVector2, surfaceBasisVector1, a1r, a2r, surfaceBasisVector3, J1);
                    var a3s = CalculateA3r(nurbs, k, j, surfaceBasisVector2, surfaceBasisVector1, a1s, a2s, surfaceBasisVector3, J1);


                    #region B

                    var term1_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            var temp = a1r.GetColumn(m).CrossProduct(a2s.GetColumn(n)) +
                                       a1s.GetColumn(n).CrossProduct(a2r.GetColumn(m));
                            temp.ScaleIntoThis(J1);
                            term1_532[m, n] = temp;
                        }
                    }

                    var term2_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        // 5.24 Kiendl Thesis
                        var a3r_dashed = a1r.GetColumn(m).CrossProduct(surfaceBasisVector2) +
                                         surfaceBasisVector1.CrossProduct(a2r.GetColumn(m));
                        for (int n = 0; n < 3; n++)
                        {
                            //TODO: a3s_dashed, a3r_dashed calculated out of the loop for all cp
                            var a3s_dashed = a1s.GetColumn(n).CrossProduct(surfaceBasisVector2) +
                                             surfaceBasisVector1.CrossProduct(a2s.GetColumn(n));
                            // 5.25 Kiendl Thesis
                            var term_525 = surfaceBasisVector3 * a3s_dashed;
                            term2_532[m, n] = a3r_dashed.Scale(-term_525 / J1 / J1);
                        }
                    }

                    var term3_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        var a3r_dashed = a1r.GetColumn(m).CrossProduct(surfaceBasisVector2) +
                                         surfaceBasisVector1.CrossProduct(a2r.GetColumn(m));
                        for (int n = 0; n < 3; n++)
                        {
                            var a3s_dashed = a1s.GetColumn(n).CrossProduct(surfaceBasisVector2) +
                                             surfaceBasisVector1.CrossProduct(a2s.GetColumn(n));
                            var term_525 = surfaceBasisVector3 * a3r_dashed;
                            term3_532[m, n] = a3s_dashed.Scale(-term_525 / J1 / J1);
                        }
                    }

                    var term4_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        var a3r_dashed = a1r.GetColumn(m).CrossProduct(surfaceBasisVector2) +
                                         surfaceBasisVector1.CrossProduct(a2r.GetColumn(m));
                        for (int n = 0; n < 3; n++)
                        {
                            var a3s_dashed = a1s.GetColumn(n).CrossProduct(surfaceBasisVector2) +
                                             surfaceBasisVector1.CrossProduct(a2s.GetColumn(n));
                            // term 5_31
                            var a3_rs = ((term1_532[m, n] * J1) *
                                        (surfaceBasisVector3 * J1)
                                        + a3r_dashed * a3s_dashed) / J1 -
                                        ((a3r_dashed * surfaceBasisVector3 * J1) *
                                        (a3s_dashed * surfaceBasisVector3 * J1)) / J1 / J1 / J1;

                            term4_532[m, n] = surfaceBasisVector3.Scale(-a3_rs / J1);
                        }
                    }

                    var term5_532 = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        var a3r_dashed = a1r.GetColumn(m).CrossProduct(surfaceBasisVector2) +
                                         surfaceBasisVector1.CrossProduct(a2r.GetColumn(m));
                        var term_525_r = surfaceBasisVector3 * a3r_dashed;
                        for (int n = 0; n < 3; n++)
                        {
                            var a3s_dashed = a1s.GetColumn(n).CrossProduct(surfaceBasisVector2) +
                                             surfaceBasisVector1.CrossProduct(a2s.GetColumn(n));
                            var term_525_s = surfaceBasisVector3 * a3s_dashed;

                            term5_532[m, n] = surfaceBasisVector3.
                                Scale(2 / J1 / J1 * term_525_r * term_525_s);
                        }
                    }

                    var a3rs = new Vector[3, 3];
                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            a3rs[m, n] = term1_532[m, n] + term2_532[m, n] + term3_532[m, n] + term4_532[m, n] +
                                         term5_532[m, n];
                        }
                    }

                    #endregion B
                    
                    var ktemp = Matrix3by3.CreateZero();
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            var Bab_rs = new double[3]
                            {(a11r.GetColumn(l) * a3s.GetColumn(m) + a11s.GetColumn(m) * a3r.GetColumn(l)),
                                (a22r.GetColumn(l) * a3s.GetColumn(m) + a22s.GetColumn(m) * a3r.GetColumn(l)),
                                (a12r.GetColumn(l) * a3s.GetColumn(m) + a12s.GetColumn(m) * a3r.GetColumn(l)) * 2};

                            Bab_rs[0] += surfaceBasisVectorDerivative1 * a3rs[l, m];
                            Bab_rs[1] += surfaceBasisVectorDerivative2 * a3rs[l, m];
                            Bab_rs[2] += surfaceBasisVectorDerivative12 * a3rs[l, m];


                            ktemp[l, m] += Vector.CreateFromArray(bendingMoments) * Vector.CreateFromArray(Bab_rs);
                        }
                    }

                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            KbendingNL[i * 3 + l, k * 3 + m] -= ktemp[l, m];
                        }
                    }
                }
			}

			return KbendingNL;
		}


        private Matrix3by3 CalculateA3r(Nurbs2D nurbs, int i, int j,
            Vector surfaceBasisVector2, Vector surfaceBasisVector1, Matrix3by3 a1r, Matrix3by3 a2r, Vector surfaceBasisVector3, double J1)
        {
            Matrix3by3 da3_tilde_dr = Matrix3by3.CreateZero(); //r1, r2 kai r3 sthles

            for (int i1 = 0; i1 < 3; i1++)
            {
                var col = (a1r.GetColumn(i1).CrossProduct(surfaceBasisVector2) + surfaceBasisVector1.CrossProduct(a2r.GetColumn(i1)));
                for (int i2 = 0; i2 < 3; i2++)
                {
                    da3_tilde_dr[i2, i1] = col[i2];
                }

            }

            double[] dnorma3_dr = new double[3]; //r1, r2 kai r3 sthles
            for (int i1 = 0; i1 < 3; i1++)
            {
                dnorma3_dr[i1] = surfaceBasisVector3.DotProduct(da3_tilde_dr.GetColumn(i1));
            }

            Matrix3by3 da3_unit_dr = Matrix3by3.CreateZero();
            for (int i1 = 0; i1 < 3; i1++)
            {
                var col = (1 / J1) * (da3_tilde_dr.GetColumn(i1) - surfaceBasisVector3.Scale(dnorma3_dr[i1]));
                for (int i2 = 0; i2 < 3; i2++)
                {
                    da3_unit_dr[i2, i1] = col[i2];
                }
            }

            return da3_unit_dr;
        }




        internal Matrix CalculateKmembraneNL(ControlPoint[] controlPoints, double[] membraneForces, Nurbs2D nurbs, int j)
		{
			var kmembraneNl =
				Matrix.CreateZero(controlPoints.Length * 3, controlPoints.Length * 3);

			for (var i = 0; i < controlPoints.Length; i++)
			{
				var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.DerivativeValuesKsi[i, j]);
				var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.DerivativeValuesHeta[i, j]);
				for (int k = 0; k < controlPoints.Length; k++)
				{
					var a1s = Matrix3by3.CreateIdentity().Scale(nurbs.DerivativeValuesKsi[k, j]);
					var a2s = Matrix3by3.CreateIdentity().Scale(nurbs.DerivativeValuesHeta[k, j]);

					var klocal = membraneForces[0] * a1r * a1s + membraneForces[1] * a2r * a2s +
								 membraneForces[2] * (a1r * a2s + a1s * a2r);

					for (int l = 0; l < 3; l++)
					{
						for (int m = 0; m < 3; m++)
						{
							kmembraneNl[i * 3 + l, k * 3 + m] += klocal[l, m];
						}
					}
				}
			}

			return kmembraneNl;
		}

		private double[,] CreateDiagonal3by3WithValue(double value)
		{
			var matrix = new double[3, 3];
			matrix[0, 0] = value;
			matrix[1, 1] = value;
			matrix[2, 2] = value;
			return matrix;
		}

		
		internal double[,] CalculateMembraneDeformationMatrix(ControlPoint[] controlPoints, Nurbs2D nurbs, int j,
			double[] surfaceBasisVector1,
			double[] surfaceBasisVector2)
		{
			var bmembrane = new double[3, controlPoints.Length * 3];
			for (int column = 0; column < controlPoints.Length * 3; column += 3)
			{
				bmembrane[0, column] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[0];
				bmembrane[0, column + 1] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[1];
				bmembrane[0, column + 2] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[2];

				bmembrane[1, column] = nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[0];
				bmembrane[1, column + 1] = nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[1];
				bmembrane[1, column + 2] = nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[2];

				bmembrane[2, column] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector2[0]+
                                       nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector1[0] ;
                bmembrane[2, column + 1] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector2[1] +
                                           nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector1[1];
				bmembrane[2, column + 2] = nurbs.DerivativeValuesKsi[column / 3, j] * surfaceBasisVector2[2] +
                                           nurbs.DerivativeValuesHeta[column / 3, j] * surfaceBasisVector1[2];
			}

			return bmembrane;
		}

		private double[,] CopyConstitutiveMatrix(double[,] f)
		{
			var g = new double[f.GetLength(0), f.GetLength(1)];
			Array.Copy(f, 0, g, 0, f.Length);
			return g;
		}

		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(NurbsKirchhoffLoveShellElementSectionNL shellElement)
		{
			var gauss = new GaussQuadrature();
			var medianSurfaceGP = gauss.CalculateElementGaussPoints(_degreeKsi, _degreeHeta, shellElement.Knots.ToList());
			return medianSurfaceGP;
		}

		private const int ThicknessIntegrationDegree = 2;

		public double Thickness { get; set; }

        internal readonly Nurbs2D _nurbs;
        private readonly int _degreeKsi;
        private int _degreeHeta;
    }
}
