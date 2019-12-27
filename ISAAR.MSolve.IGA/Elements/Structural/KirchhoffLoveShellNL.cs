using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using System.Text;
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
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Elements.Structural
{
    public class KirchhoffLoveShellNL : Element, IStructuralIsogeometricElement, ISurfaceLoadedElement
    {
        protected static readonly IDofType[] ControlPointDofTypes =
            {StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ};

        private IDofType[][] dofTypes;

        private Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>> thicknessIntegrationPoints =
            new Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>>();

        internal Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>
            materialsAtThicknessGP =
                new Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>();

        private bool isInitialized;
        internal double[] _solution;

        public KirchhoffLoveShellNL(IShellMaterial shellMaterial, IList<Knot> elementKnots,
            IShapeFunction2D shapeFunctions, IList<ControlPoint> elementControlPoints, Patch patch, double thickness,
            int degreeKsi, int degreeHeta)
        {
            Contract.Requires(shellMaterial != null);
            this.Patch = patch;
            this.Thickness = thickness;
            _degreeKsi = degreeKsi;
            _degreeHeta = degreeHeta;
            foreach (var knot in elementKnots)
            {
                if (!KnotsDictionary.ContainsKey(knot.ID))
                    this.KnotsDictionary.Add(knot.ID, knot);
            }

            _shapeFunctions = shapeFunctions;
            _solution = new double[3 * elementControlPoints.Count];
            
            CreateElementGaussPoints(this);
            foreach (var medianSurfaceGP in thicknessIntegrationPoints.Keys)
            {
                materialsAtThicknessGP.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IShellMaterial>());
                foreach (var point in thicknessIntegrationPoints[medianSurfaceGP])
                {
                    materialsAtThicknessGP[medianSurfaceGP].Add(point, shellMaterial.Clone());
                }
            }

            _controlPoints = elementControlPoints.ToArray();
        }

        public CellType CellType { get; } = CellType.Unknown;

        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public bool MaterialModified => false;

        public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) =>
            throw new NotImplementedException();

        public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
        {
            throw new NotImplementedException();
            //var nurbsElement = (NurbsKirchhoffLoveShellElementNL) element;
            //var knotParametricCoordinatesKsi = Vector.CreateFromArray(Knots.Select(k => k.Ksi).ToArray());
            //var knotParametricCoordinatesHeta = Vector.CreateFromArray(Knots.Select(k => k.Heta).ToArray());

            //var shapeFunctions = new Nurbs2D(nurbsElement, nurbsElement.ControlPoints.ToArray(), knotParametricCoordinatesKsi,
            //    knotParametricCoordinatesHeta);

            //var knotDisplacements = new double[4, 3];
            //var paraviewKnotRenumbering = new int[] {0, 3, 1, 2};
            //for (var j = 0; j < knotDisplacements.GetLength(0); j++)
            //{
            //    for (int i = 0; i < element.ControlPoints.Count(); i++)
            //    {
            //        knotDisplacements[paraviewKnotRenumbering[j], 0] +=
            //            shapeFunctions.NurbsValues[i, j] * localDisplacements[i, 0];
            //        knotDisplacements[paraviewKnotRenumbering[j], 1] +=
            //            shapeFunctions.NurbsValues[i, j] * localDisplacements[i, 1];
            //        knotDisplacements[paraviewKnotRenumbering[j], 2] +=
            //            shapeFunctions.NurbsValues[i, j] * localDisplacements[i, 2];
            //    }
            //}

            //return knotDisplacements;
        }

        public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            var shellElement = (Elements.NurbsKirchhoffLoveShellElementNL)element;
            var elementNodalForces = new double[shellElement.ControlPointsDictionary.Count * 3];

            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(_controlPoints);
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();

            var Bmembrane = new double[3, _controlPoints.Length * 3];
            var Bbending = new double[3, _controlPoints.Length * 3];
            var numberOfControlPoints = _controlPoints.Length;
            var MembraneForces = new Elements.Forces();
            var BendingMoments = new Elements.Forces();

            for (int j = 0; j < gaussPoints.Length; j++)
            {
                CalculateJacobian(newControlPoints, _shapeFunctions, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(newControlPoints, _shapeFunctions, j);

                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

                var surfaceBasisVector3 = new[]
                {
                    surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                    surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                    surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
                };

                var J1 = Math.Sqrt(surfaceBasisVector3[0] * surfaceBasisVector3[0] +
                                   surfaceBasisVector3[1] * surfaceBasisVector3[1] +
                                   surfaceBasisVector3[2] * surfaceBasisVector3[2]);

                surfaceBasisVector3[0] /= J1;
                surfaceBasisVector3[1] /= J1;
                surfaceBasisVector3[2] /= J1;


                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                CalculateMembraneDeformationMatrix(numberOfControlPoints, _shapeFunctions, j, surfaceBasisVector1,
                    surfaceBasisVector2, Bmembrane);
                CalculateBendingDeformationMatrix(numberOfControlPoints, surfaceBasisVector3, _shapeFunctions, j,
                    surfaceBasisVector2,
                    surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12, Bbending);

                IntegratedStressesOverThickness(gaussPoints[j], ref MembraneForces, ref BendingMoments);

                var wfactor = InitialJ1[j] * gaussPoints[j].WeightFactor;

                for (int i = 0; i < Bmembrane.GetLength(1); i++)
                {
                    elementNodalForces[i] +=
                        (Bmembrane[0, i] * MembraneForces.v0 * wfactor + Bbending[0, i] * BendingMoments.v0 * wfactor) +
                        (Bmembrane[1, i] * MembraneForces.v1 * wfactor + Bbending[1, i] * BendingMoments.v1 * wfactor) +
                        (Bmembrane[2, i] * MembraneForces.v2 * wfactor + Bbending[2, i] * BendingMoments.v2 * wfactor);
                }
            }

            return elementNodalForces;
        }

        public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
            NeumannBoundaryCondition neumann) => throw new NotImplementedException();

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
            NeumannBoundaryCondition neumann) => throw new NotImplementedException();

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
            PressureBoundaryCondition pressure) => throw new NotImplementedException();

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
            PressureBoundaryCondition pressure) => throw new NotImplementedException();

        double[,] jacobianMatrix = new double[2, 3];
        public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
            double[] localdDisplacements)
        {
            var shellElement = (Elements.NurbsKirchhoffLoveShellElementNL)element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();

            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(elementControlPoints);
            var midsurfaceGP = materialsAtThicknessGP.Keys.ToArray();

            for (var j = 0; j < midsurfaceGP.Length; j++)
            {
                CalculateJacobian(newControlPoints, _shapeFunctions, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(newControlPoints, _shapeFunctions, j);

                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

                var surfaceBasisVector3 = new[]
                {
                    surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                    surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                    surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
                };

                var J1 = Math.Sqrt(surfaceBasisVector3[0] * surfaceBasisVector3[0] +
                                   surfaceBasisVector3[1] * surfaceBasisVector3[1] +
                                   surfaceBasisVector3[2] * surfaceBasisVector3[2]);

                surfaceBasisVector3[0] /= J1;
                surfaceBasisVector3[1] /= J1;
                surfaceBasisVector3[2] /= J1;

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

                var b11 = surfaceBasisVectorDerivative1[0] * surfaceBasisVector3[0] +
                          surfaceBasisVectorDerivative1[1] * surfaceBasisVector3[1] +
                          surfaceBasisVectorDerivative1[2] * surfaceBasisVector3[2];

                var b22 = surfaceBasisVectorDerivative2[0] * surfaceBasisVector3[0] +
                          surfaceBasisVectorDerivative2[1] * surfaceBasisVector3[1] +
                          surfaceBasisVectorDerivative2[2] * surfaceBasisVector3[2];

                var b12 = surfaceBasisVectorDerivative12[0] * surfaceBasisVector3[0] +
                          surfaceBasisVectorDerivative12[1] * surfaceBasisVector3[1] +
                          surfaceBasisVectorDerivative12[2] * surfaceBasisVector3[2];

                var bendingStrain = new double[] { b11 - B11, b22 - B22, 2 * b12 - 2 * B12 };

                //var bendingStrain = new double[] { -(b11 - B11), -(b22 - B22), -(2 * b12 - 2 * B12) };

                foreach (var keyValuePair in materialsAtThicknessGP[midsurfaceGP[j]])
                {
                    var thicknessPoint = keyValuePair.Key;
                    var material = keyValuePair.Value;
                    var gpStrain = new double[bendingStrain.Length];
                    var z = thicknessPoint.Zeta;
                    for (var i = 0; i < bendingStrain.Length; i++)
                    {
                        gpStrain[i] += membraneStrain[i] + bendingStrain[i] * z;
                    }

                    material.UpdateMaterial(gpStrain);
                }
            }

            return new Tuple<double[], double[]>(new double[0], new double[0]);
        }

        public Dictionary<int, double> CalculateSurfaceDistributedLoad(Element element, IDofType loadedDof,
            double loadMagnitude)
        {
            var shellElement = (KirchhoffLoveShellNL)element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();
            var gaussPoints = CreateElementGaussPoints(shellElement);
            var distributedLoad = new Dictionary<int, double>();

            for (var j = 0; j < gaussPoints.Count; j++)
            {
                CalculateJacobian(elementControlPoints, _shapeFunctions, j, jacobianMatrix);
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
                                                  _shapeFunctions.Values[i, j] * gaussPoints[j].WeightFactor;
                    }
                    else
                    {
                        distributedLoad.Add(dofId,
                            loadMagnitude * _shapeFunctions.Values[i, j] * J1 * gaussPoints[j].WeightFactor);
                    }
                }
            }

            return distributedLoad;
        }

        public Dictionary<int, double> CalculateSurfacePressure(Element element, double pressureMagnitude)
        {
            var shellElement = (KirchhoffLoveShellNL)element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();
            var gaussPoints = CreateElementGaussPoints(shellElement);
            var pressureLoad = new Dictionary<int, double>();

            for (var j = 0; j < gaussPoints.Count; j++)
            {
                CalculateJacobian(elementControlPoints, _shapeFunctions, j, jacobianMatrix);
                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
                var surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
                var J1 = surfaceBasisVector3.Norm2();
                surfaceBasisVector3.ScaleIntoThis(1 / J1);

                for (int i = 0; i < elementControlPoints.Length; i++)
                {
                    for (int k = 0; k < ControlPointDofTypes.Length; k++)
                    {
                        int dofId = element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[i],
                            ControlPointDofTypes[k]];

                        if (pressureLoad.ContainsKey(dofId))
                        {
                            pressureLoad[dofId] += pressureMagnitude * surfaceBasisVector3[k] *
                                                   _shapeFunctions.Values[i, j] * gaussPoints[j].WeightFactor;
                        }
                        else
                        {
                            pressureLoad.Add(dofId,
                                pressureMagnitude * surfaceBasisVector3[k] * _shapeFunctions.Values[i, j] *
                                gaussPoints[j].WeightFactor);
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

        internal (double[,] MembraneConstitutiveMatrix, double[,] BendingConstitutiveMatrix, double[,]
            CouplingConstitutiveMatrix) IntegratedConstitutiveOverThickness(GaussLegendrePoint3D midSurfaceGaussPoint)
        {
            var MembraneConstitutiveMatrix = new double[3, 3];
            var BendingConstitutiveMatrix = new double[3, 3];
            var CouplingConstitutiveMatrix = new double[3, 3];

            foreach (var keyValuePair in materialsAtThicknessGP[midSurfaceGaussPoint])
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

        internal void IntegratedStressesOverThickness(
            GaussLegendrePoint3D midSurfaceGaussPoint, ref Elements.Forces MembraneForces, ref Elements.Forces BendingMoments)
        {
            MembraneForces = new Elements.Forces();
            BendingMoments = new Elements.Forces();
            var thicknessPoints = thicknessIntegrationPoints[midSurfaceGaussPoint];

            for (int i = 0; i < thicknessPoints.Count; i++)
            {
                var thicknessPoint = thicknessPoints[i];
                var material = materialsAtThicknessGP[midSurfaceGaussPoint][thicknessPoints[i]];
                var w = thicknessPoint.WeightFactor;
                var z = thicknessPoint.Zeta;
                MembraneForces.v0 += material.Stresses[0] * w;
                MembraneForces.v1 += material.Stresses[1] * w;
                MembraneForces.v2 += material.Stresses[2] * w;

                BendingMoments.v0 -= material.Stresses[0] * w * z;
                BendingMoments.v1 -= material.Stresses[1] * w * z;
                BendingMoments.v2 -= material.Stresses[2] * w * z;
            }
        }

        public IMatrix MassMatrix(IElement element) => throw new NotImplementedException();

        public void ResetMaterialModified() => throw new NotImplementedException();

        public void SaveMaterialState()
        {
            foreach (var gp in materialsAtThicknessGP.Keys)
            {
                foreach (var material in materialsAtThicknessGP[gp].Values)
                {
                    material.SaveState();
                }
            }
        }

        internal IShapeFunction2D _shapeFunctions;

        public IMatrix StiffnessMatrix(IElement element)
        {
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();

            if (!isInitialized)
            {
                CalculateInitialConfigurationData(_controlPoints, _shapeFunctions, gaussPoints);
                isInitialized = true;
            }

            var elementControlPoints = CurrentControlPoint(_controlPoints);

            var bRows = 3;
            var bCols = elementControlPoints.Length * 3;
            var stiffnessMatrix = new double[bCols, bCols];
            var KmembraneNL = new double[bCols, bCols];
            var KbendingNL = new double[bCols, bCols];
            var Bmembrane = new double[bRows, bCols];
            var Bbending = new double[bRows, bCols];

            var BmTranspose = new double[bCols, bRows];
            var BbTranspose = new double[bCols, bRows];

            var BmTransposeMultStiffness = new double[bCols, bRows];
            var BbTransposeMultStiffness = new double[bCols, bRows];
            var BmbTransposeMultStiffness = new double[bCols, bRows];
            var BbmTransposeMultStiffness = new double[bCols, bRows];
            var MembraneForces = new Elements.Forces();
            var BendingMoments = new Elements.Forces();

            for (int j = 0; j < gaussPoints.Length; j++)
            {
                CalculateJacobian(elementControlPoints, _shapeFunctions, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(elementControlPoints, _shapeFunctions, j);
                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

                var surfaceBasisVector3 = new[]
                {
                    surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                    surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                    surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
                };

                var J1 = Math.Sqrt(surfaceBasisVector3[0] * surfaceBasisVector3[0] +
                                   surfaceBasisVector3[1] * surfaceBasisVector3[1] +
                                   surfaceBasisVector3[2] * surfaceBasisVector3[2]);

                surfaceBasisVector3[0] /= J1;
                surfaceBasisVector3[1] /= J1;
                surfaceBasisVector3[2] /= J1;

                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                double wFactor = InitialJ1[j] * gaussPoints[j].WeightFactor;

                CalculateLinearStiffness(elementControlPoints, _shapeFunctions, j, surfaceBasisVector1, surfaceBasisVector2,
                    Bmembrane, surfaceBasisVector3, surfaceBasisVectorDerivative1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12, Bbending, gaussPoints, BmTranspose, bRows, bCols, BbTranspose,
                    wFactor, BmTransposeMultStiffness, BbTransposeMultStiffness, BmbTransposeMultStiffness,
                    BbmTransposeMultStiffness, stiffnessMatrix);
                CalculateNonLinearStiffness(gaussPoints, j, KmembraneNL, bCols, KbendingNL, elementControlPoints, _shapeFunctions,
                    surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, surfaceBasisVectorDerivative1,
                    surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, J1, stiffnessMatrix, wFactor,
                    ref MembraneForces, ref BendingMoments);
            }

            return Matrix.CreateFromArray(stiffnessMatrix);
        }

        private void CalculateNonLinearStiffness(GaussLegendrePoint3D[] gaussPoints, int j, double[,] KmembraneNL, int bCols,
            double[,] KbendingNL, ControlPoint[] elementControlPoints, IShapeFunction2D shapeFunctions, double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3, double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12, double J1,
            double[,] stiffnessMatrix, double wFactor, ref Elements.Forces MembraneForces, ref Elements.Forces BendingMoments)
        {
            IntegratedStressesOverThickness(gaussPoints[j], ref MembraneForces, ref BendingMoments);

            Array.Clear(KmembraneNL, 0, bCols * bCols);
            Array.Clear(KbendingNL, 0, bCols * bCols);

            CalculateKmembraneNL(elementControlPoints, ref MembraneForces, shapeFunctions, j, KmembraneNL);
            CalculateKbendingNL(elementControlPoints, ref BendingMoments, shapeFunctions,
                surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3,
                surfaceBasisVectorDerivative1,
                surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12, J1, j, KbendingNL);

            for (var i = 0; i < stiffnessMatrix.GetLength(0); i++)
            {
                for (var k = 0; k < stiffnessMatrix.GetLength(1); k++)
                {
                    stiffnessMatrix[i, k] += (KmembraneNL[i, k] + KbendingNL[i, k]) * wFactor;
                }
            }
        }

        private void CalculateLinearStiffness(ControlPoint[] elementControlPoints, IShapeFunction2D shapeFunctions, int j,
            double[] surfaceBasisVector1, double[] surfaceBasisVector2, double[,] Bmembrane, double[] surfaceBasisVector3,
            double[] surfaceBasisVectorDerivative1, double J1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double[,] Bbending, GaussLegendrePoint3D[] gaussPoints,
            double[,] BmTranspose, int bRows, int bCols, double[,] BbTranspose, double wFactor,
            double[,] BmTransposeMultStiffness, double[,] BbTransposeMultStiffness, double[,] BmbTransposeMultStiffness,
            double[,] BbmTransposeMultStiffness, double[,] stiffnessMatrix)
        {
            CalculateMembraneDeformationMatrix(elementControlPoints.Length, shapeFunctions, j, surfaceBasisVector1,
                surfaceBasisVector2, Bmembrane);
            CalculateBendingDeformationMatrix(elementControlPoints.Length, surfaceBasisVector3, shapeFunctions, j,
                surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12, Bbending);

            var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
                IntegratedConstitutiveOverThickness(gaussPoints[j]);


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


        internal ControlPoint[] CurrentControlPoint(ControlPoint[] controlPoints)
        {
            var cp = new ControlPoint[controlPoints.Length];

            for (int i = 0; i < controlPoints.Length; i++)
            {
                cp[i] = new ControlPoint()
                {
                    ID = controlPoints[i].ID,
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

        internal double[,] CalculateHessian(ControlPoint[] controlPoints, IShapeFunction2D shapeFunctions, int j)
        {
            var hessianMatrix = new double[3, 3];
            for (var k = 0; k < controlPoints.Length; k++)
            {
                hessianMatrix[0, 0] +=
                    shapeFunctions.SecondDerivativeValuesKsi[k, j] * controlPoints[k].X;
                hessianMatrix[0, 1] +=
                    shapeFunctions.SecondDerivativeValuesKsi[k, j] * controlPoints[k].Y;
                hessianMatrix[0, 2] +=
                    shapeFunctions.SecondDerivativeValuesKsi[k, j] * controlPoints[k].Z;
                hessianMatrix[1, 0] +=
                    shapeFunctions.SecondDerivativeValuesHeta[k, j] * controlPoints[k].X;
                hessianMatrix[1, 1] +=
                    shapeFunctions.SecondDerivativeValuesHeta[k, j] * controlPoints[k].Y;
                hessianMatrix[1, 2] +=
                    shapeFunctions.SecondDerivativeValuesHeta[k, j] * controlPoints[k].Z;
                hessianMatrix[2, 0] +=
                    shapeFunctions.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].X;
                hessianMatrix[2, 1] +=
                    shapeFunctions.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].Y;
                hessianMatrix[2, 2] +=
                    shapeFunctions.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].Z;
            }

            return hessianMatrix;
        }

        internal void CalculateJacobian(ControlPoint[] controlPoints, IShapeFunction2D shapeFunctions, int j, double[,] jacobianOut)
        {
            jacobianOut[0, 0] = jacobianOut[0, 1] = jacobianOut[0, 2] =
                jacobianOut[1, 0] = jacobianOut[1, 1] = jacobianOut[1, 2] = 0.0;
            for (var k = 0; k < controlPoints.Length; k++)
            {
                jacobianOut[0, 0] += shapeFunctions.DerivativeValuesKsi[k, j] * controlPoints[k].X;
                jacobianOut[0, 1] += shapeFunctions.DerivativeValuesKsi[k, j] * controlPoints[k].Y;
                jacobianOut[0, 2] += shapeFunctions.DerivativeValuesKsi[k, j] * controlPoints[k].Z;
                jacobianOut[1, 0] += shapeFunctions.DerivativeValuesHeta[k, j] * controlPoints[k].X;
                jacobianOut[1, 1] += shapeFunctions.DerivativeValuesHeta[k, j] * controlPoints[k].Y;
                jacobianOut[1, 2] += shapeFunctions.DerivativeValuesHeta[k, j] * controlPoints[k].Z;
            }
        }

        internal double[] CalculateSurfaceBasisVector1(double[,] Matrix, int row)
        {
            var surfaceBasisVector1 = new double[3];
            surfaceBasisVector1[0] = Matrix[row, 0];
            surfaceBasisVector1[1] = Matrix[row, 1];
            surfaceBasisVector1[2] = Matrix[row, 2];
            return surfaceBasisVector1;
        }

        internal void CalculateBendingDeformationMatrix(int controlPointsCount, double[] surfaceBasisVector3,
            IShapeFunction2D shapeFunctions, int j, double[] surfaceBasisVector2, double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVector1, double J1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double[,] BbendingOut)
        {
            var s10 = surfaceBasisVector1[0];
            var s11 = surfaceBasisVector1[1];
            var s12 = surfaceBasisVector1[2];

            var s20 = surfaceBasisVector2[0];
            var s21 = surfaceBasisVector2[1];
            var s22 = surfaceBasisVector2[2];

            var s30 = surfaceBasisVector3[0];
            var s31 = surfaceBasisVector3[1];
            var s32 = surfaceBasisVector3[2];

            var s11_0 = surfaceBasisVectorDerivative1[0];
            var s11_1 = surfaceBasisVectorDerivative1[1];
            var s11_2 = surfaceBasisVectorDerivative1[2];

            var s22_0 = surfaceBasisVectorDerivative2[0];
            var s22_1 = surfaceBasisVectorDerivative2[1];
            var s22_2 = surfaceBasisVectorDerivative2[2];

            var s12_0 = surfaceBasisVectorDerivative12[0];
            var s12_1 = surfaceBasisVectorDerivative12[1];
            var s12_2 = surfaceBasisVectorDerivative12[2];

            for (int column = 0; column < controlPointsCount * 3; column += 3)
            {
                var dksi = shapeFunctions.DerivativeValuesKsi[column / 3, j];
                var dheta = shapeFunctions.DerivativeValuesHeta[column / 3, j];
                var d2Ksi = shapeFunctions.SecondDerivativeValuesKsi[column / 3, j];
                var d2Heta = shapeFunctions.SecondDerivativeValuesHeta[column / 3, j];
                var d2KsiHeta = shapeFunctions.SecondDerivativeValuesKsiHeta[column / 3, j];

                BbendingOut[0, column] = -d2Ksi * s30 - ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s11 * s11_2 - s12 * s11_1) + dksi * (s21 * s11_2 - s22 * s11_1)) / J1;
                BbendingOut[0, column + 1] = ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s10 * s11_2 - s12 * s11_0) + dksi * (s20 * s11_2 - s22 * s11_0)) / J1 - d2Ksi * s31;
                BbendingOut[0, column + 2] = -d2Ksi * s32 - ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s10 * s11_1 - s11 * s11_0) + dksi * (s20 * s11_1 - s21 * s11_0)) / J1;

                BbendingOut[1, column] = -d2Heta * s30 - ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s11 * s22_2 - s12 * s22_1) + dksi * (s21 * s22_2 - s22 * s22_1)) / J1;
                BbendingOut[1, column + 1] = ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s10 * s22_2 - s12 * s22_0) + dksi * (s20 * s22_2 - s22 * s22_0)) / J1 - d2Heta * s31;
                BbendingOut[1, column + 2] = -d2Heta * s32 - ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s10 * s22_1 - s11 * s22_0) + dksi * (s20 * s22_1 - s21 * s22_0)) / J1;

                BbendingOut[2, column] = -2 * d2KsiHeta * s30 - (2 * ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s11 * s12_2 - s12 * s12_1) + dksi * (s21 * s12_2 - s22 * s12_1))) / J1;
                BbendingOut[2, column + 1] = (2 * ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s10 * s12_2 - s12 * s12_0) + dksi * (s20 * s12_2 - s22 * s12_0))) / J1 - 2 * d2KsiHeta * s31;
                BbendingOut[2, column + 2] = -2 * d2KsiHeta * s32 - (2 * ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s10 * s12_1 - s11 * s12_0) + dksi * (s20 * s12_1 - s21 * s12_0))) / J1;
            }
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

        internal void CalculateInitialConfigurationData(ControlPoint[] controlPoints,
            IShapeFunction2D shapeFunctions, IList<GaussLegendrePoint3D> gaussPoints)
        {
            var numberOfGP = gaussPoints.Count;
            InitialJ1 = new double[numberOfGP];
            initialSurfaceBasisVectors1 = new double[numberOfGP][];
            initialSurfaceBasisVectors2 = new double[numberOfGP][];
            initialUnitSurfaceBasisVectors3 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative1 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative2 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative12 = new double[numberOfGP][];

            for (int j = 0; j < gaussPoints.Count; j++)
            {
                CalculateJacobian(controlPoints, shapeFunctions, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(controlPoints, shapeFunctions, j);
                initialSurfaceBasisVectors1[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
                initialSurfaceBasisVectors2[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
                var s3 = CalculateCrossProduct(initialSurfaceBasisVectors1[j], initialSurfaceBasisVectors2[j]);
                var norm = s3.Sum(t => t * t);
                InitialJ1[j] = Math.Sqrt(norm);
                var vector3 = CalculateCrossProduct(initialSurfaceBasisVectors1[j], initialSurfaceBasisVectors2[j]);
                initialUnitSurfaceBasisVectors3[j] = new double[]
                {
                    vector3[0] / InitialJ1[j],
                    vector3[1] / InitialJ1[j],
                    vector3[2] / InitialJ1[j],
                };

                initialSurfaceBasisVectorDerivative1[j] = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                initialSurfaceBasisVectorDerivative2[j] = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                initialSurfaceBasisVectorDerivative12[j] = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                foreach (var integrationPointMaterial in materialsAtThicknessGP[gaussPoints[j]].Values)
                {
                    integrationPointMaterial.TangentVectorV1 = initialSurfaceBasisVectors1[j];
                    integrationPointMaterial.TangentVectorV2 = initialSurfaceBasisVectors2[j];
                    integrationPointMaterial.NormalVectorV3 = initialUnitSurfaceBasisVectors3[j];
                }
            }
        }

        private Elements.a3rs a3rs = new Elements.a3rs();
        private Elements.Bab_rs Bab_rs = new Elements.Bab_rs();
        private Elements.a3r a3r = new Elements.a3r();
        private Elements.a3r a3s = new Elements.a3r();
        private ControlPoint[] _controlPoints;

        internal void CalculateKbendingNL(ControlPoint[] controlPoints,
            ref Elements.Forces bendingMoments, IShapeFunction2D shapeFunctions, double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double[] surfaceBasisVectorDerivative1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double J1, int j, double[,] KbendingNLOut)
        {
            for (int i = 0; i < controlPoints.Length; i++)
            {
                var dksi_r = shapeFunctions.DerivativeValuesKsi[i, j];
                var dheta_r = shapeFunctions.DerivativeValuesHeta[i, j];

                var d2Ksi_dr2 = shapeFunctions.SecondDerivativeValuesKsi[i, j];
                var d2Heta_dr2 = shapeFunctions.SecondDerivativeValuesHeta[i, j];
                var d2KsiHeta_dr2 = shapeFunctions.SecondDerivativeValuesKsiHeta[i, j];

                for (int k = 0; k < controlPoints.Length; k++)
                {
                    var d2Ksi_ds2 = shapeFunctions.SecondDerivativeValuesKsi[k, j];
                    var d2Heta_ds2 = shapeFunctions.SecondDerivativeValuesHeta[k, j];
                    var d2KsiHeta_ds2 = shapeFunctions.SecondDerivativeValuesKsiHeta[k, j];

                    var dksi_s = shapeFunctions.DerivativeValuesKsi[k, j];
                    var dheta_s = shapeFunctions.DerivativeValuesHeta[k, j];

                    CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r,
                        dheta_r, J1, ref a3r);
                    CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_s,
                        dheta_s, J1, ref a3s);
                    a3rs = new Elements.a3rs();//Clear struct values
                    Calculate_a3rs(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, J1, dksi_r,
                        dheta_r, dksi_s, dheta_s, ref a3rs);

                    CalculateBab_rs(surfaceBasisVectorDerivative1, surfaceBasisVectorDerivative2,
                        surfaceBasisVectorDerivative12, d2Ksi_dr2, ref a3s, d2Ksi_ds2, ref a3r, ref a3rs, d2Heta_dr2,
                        d2Heta_ds2, d2KsiHeta_dr2, d2KsiHeta_ds2, ref Bab_rs);


                    KbendingNLOut[i * 3 + 0, k * 3 + 0] -= (Bab_rs.Bab_rs00_0 * bendingMoments.v0 + Bab_rs.Bab_rs00_1 * bendingMoments.v1 + Bab_rs.Bab_rs00_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 0, k * 3 + 1] -= (Bab_rs.Bab_rs01_0 * bendingMoments.v0 + Bab_rs.Bab_rs01_1 * bendingMoments.v1 + Bab_rs.Bab_rs01_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 0, k * 3 + 2] -= (Bab_rs.Bab_rs02_0 * bendingMoments.v0 + Bab_rs.Bab_rs02_1 * bendingMoments.v1 + Bab_rs.Bab_rs02_2 * bendingMoments.v2);

                    KbendingNLOut[i * 3 + 1, k * 3 + 0] -= (Bab_rs.Bab_rs10_0 * bendingMoments.v0 + Bab_rs.Bab_rs10_1 * bendingMoments.v1 + Bab_rs.Bab_rs10_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 1, k * 3 + 1] -= (Bab_rs.Bab_rs11_0 * bendingMoments.v0 + Bab_rs.Bab_rs11_1 * bendingMoments.v1 + Bab_rs.Bab_rs11_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 1, k * 3 + 2] -= (Bab_rs.Bab_rs12_0 * bendingMoments.v0 + Bab_rs.Bab_rs12_1 * bendingMoments.v1 + Bab_rs.Bab_rs12_2 * bendingMoments.v2);

                    KbendingNLOut[i * 3 + 2, k * 3 + 0] -= (Bab_rs.Bab_rs20_0 * bendingMoments.v0 + Bab_rs.Bab_rs20_1 * bendingMoments.v1 + Bab_rs.Bab_rs20_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 2, k * 3 + 1] -= (Bab_rs.Bab_rs21_0 * bendingMoments.v0 + Bab_rs.Bab_rs21_1 * bendingMoments.v1 + Bab_rs.Bab_rs21_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 2, k * 3 + 2] -= (Bab_rs.Bab_rs22_0 * bendingMoments.v0 + Bab_rs.Bab_rs22_1 * bendingMoments.v1 + Bab_rs.Bab_rs22_2 * bendingMoments.v2);

                }
            }
        }

        private static void CalculateBab_rs(double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12, double d2Ksi_dr2, ref Elements.a3r a3s,
            double d2Ksi_ds2, ref Elements.a3r a3r, ref Elements.a3rs a3rs, double d2Heta_dr2, double d2Heta_ds2, double d2KsiHeta_dr2,
            double d2KsiHeta_ds2, ref Elements.Bab_rs Bab_rsOut)
        {
            var s11_0 = surfaceBasisVectorDerivative1[0];
            var s11_1 = surfaceBasisVectorDerivative1[1];
            var s11_2 = surfaceBasisVectorDerivative1[2];

            var s22_0 = surfaceBasisVectorDerivative2[0];
            var s22_1 = surfaceBasisVectorDerivative2[1];
            var s22_2 = surfaceBasisVectorDerivative2[2];

            var s12_0 = surfaceBasisVectorDerivative12[0];
            var s12_1 = surfaceBasisVectorDerivative12[1];
            var s12_2 = surfaceBasisVectorDerivative12[2];

            Bab_rsOut.Bab_rs00_0 = d2Ksi_dr2 * a3s.a3r00 + d2Ksi_ds2 * a3r.a3r00 +
                                 s11_0 * a3rs.a3rs00_0 +
                                 s11_1 * a3rs.a3rs00_1 +
                                 s11_2 * a3rs.a3rs00_2;

            Bab_rsOut.Bab_rs00_1 = d2Heta_dr2 * a3s.a3r00 + d2Heta_ds2 * a3r.a3r00 +
                                 s22_0 * a3rs.a3rs00_0 +
                                 s22_1 * a3rs.a3rs00_1 +
                                 s22_2 * a3rs.a3rs00_2;

            Bab_rsOut.Bab_rs00_2 = (d2KsiHeta_dr2 * a3s.a3r00 + d2KsiHeta_ds2 * a3r.a3r00) * 2 +
                                 s12_0 * a3rs.a3rs00_0 +
                                 s12_1 * a3rs.a3rs00_1 +
                                 s12_2 * a3rs.a3rs00_2;

            Bab_rsOut.Bab_rs01_0 = d2Ksi_dr2 * a3s.a3r01 + d2Ksi_ds2 * a3r.a3r10 +
                                 s11_0 * a3rs.a3rs01_0 +
                                 s11_1 * a3rs.a3rs01_1 +
                                 s11_2 * a3rs.a3rs01_2;
            Bab_rsOut.Bab_rs01_1 = d2Heta_dr2 * a3s.a3r01 + d2Heta_ds2 * a3r.a3r10 +
                                 s22_0 * a3rs.a3rs01_0 +
                                 s22_1 * a3rs.a3rs01_1 +
                                 s22_2 * a3rs.a3rs01_2;
            Bab_rsOut.Bab_rs01_2 = (d2KsiHeta_dr2 * a3s.a3r01 + d2KsiHeta_ds2 * a3r.a3r10) * 2 +
                                s12_0 * a3rs.a3rs01_0 +
                                s12_1 * a3rs.a3rs01_1 +
                                s12_2 * a3rs.a3rs01_2;

            Bab_rsOut.Bab_rs02_0 = d2Ksi_dr2 * a3s.a3r02 + d2Ksi_ds2 * a3r.a3r20 +
                                 s11_0 * a3rs.a3rs02_0 +
                                 s11_1 * a3rs.a3rs02_1 +
                                 s11_2 * a3rs.a3rs02_2;
            Bab_rsOut.Bab_rs02_1 = d2Heta_dr2 * a3s.a3r02 + d2Heta_ds2 * a3r.a3r20 +
                                 s22_0 * a3rs.a3rs02_0 +
                                 s22_1 * a3rs.a3rs02_1 +
                                 s22_2 * a3rs.a3rs02_2;
            Bab_rsOut.Bab_rs02_2 = (d2KsiHeta_dr2 * a3s.a3r02 + d2KsiHeta_ds2 * a3r.a3r20) * 2 +
                                 s12_0 * a3rs.a3rs02_0 +
                                 s12_1 * a3rs.a3rs02_1 +
                                 s12_2 * a3rs.a3rs02_2;

            Bab_rsOut.Bab_rs10_0 = d2Ksi_dr2 * a3s.a3r10 + d2Ksi_ds2 * a3r.a3r01 +
                                 s11_0 * a3rs.a3rs10_0 +
                                 s11_1 * a3rs.a3rs10_1 +
                                 s11_2 * a3rs.a3rs10_2;
            Bab_rsOut.Bab_rs10_1 = d2Heta_dr2 * a3s.a3r10 + d2Heta_ds2 * a3r.a3r01 +
                                 s22_0 * a3rs.a3rs10_0 +
                                 s22_1 * a3rs.a3rs10_1 +
                                 s22_2 * a3rs.a3rs10_2;
            Bab_rsOut.Bab_rs10_2 = (d2KsiHeta_dr2 * a3s.a3r10 + d2KsiHeta_ds2 * a3r.a3r01) * 2 +
                                 s12_0 * a3rs.a3rs10_0 +
                                 s12_1 * a3rs.a3rs10_1 +
                                 s12_2 * a3rs.a3rs10_2;

            Bab_rsOut.Bab_rs11_0 = d2Ksi_dr2 * a3s.a3r11 + d2Ksi_ds2 * a3r.a3r11 +
                                 s11_0 * a3rs.a3rs11_0 +
                                 s11_1 * a3rs.a3rs11_1 +
                                 s11_2 * a3rs.a3rs11_2;
            Bab_rsOut.Bab_rs11_1 = d2Heta_dr2 * a3s.a3r11 + d2Heta_ds2 * a3r.a3r11 +
                                 s22_0 * a3rs.a3rs11_0 +
                                 s22_1 * a3rs.a3rs11_1 +
                                 s22_2 * a3rs.a3rs11_2;
            Bab_rsOut.Bab_rs11_2 = (d2KsiHeta_dr2 * a3s.a3r11 + d2KsiHeta_ds2 * a3r.a3r11) * 2 +
                                 s12_0 * a3rs.a3rs11_0 +
                                 s12_1 * a3rs.a3rs11_1 +
                                 s12_2 * a3rs.a3rs11_2;

            Bab_rsOut.Bab_rs12_0 = d2Ksi_dr2 * a3s.a3r12 + d2Ksi_ds2 * a3r.a3r21 +
                                 s11_0 * a3rs.a3rs12_0 +
                                 s11_1 * a3rs.a3rs12_1 +
                                 s11_2 * a3rs.a3rs12_2;
            Bab_rsOut.Bab_rs12_1 = d2Heta_dr2 * a3s.a3r12 + d2Heta_ds2 * a3r.a3r21 +
                                 s22_0 * a3rs.a3rs12_0 +
                                 s22_1 * a3rs.a3rs12_1 +
                                 s22_2 * a3rs.a3rs12_2;
            Bab_rsOut.Bab_rs12_2 = (d2KsiHeta_dr2 * a3s.a3r12 + d2KsiHeta_ds2 * a3r.a3r21) * 2 +
                                 s12_0 * a3rs.a3rs12_0 +
                                 s12_1 * a3rs.a3rs12_1 +
                                 s12_2 * a3rs.a3rs12_2;

            Bab_rsOut.Bab_rs20_0 = d2Ksi_dr2 * a3s.a3r20 + d2Ksi_ds2 * a3r.a3r02 +
                                 s11_0 * a3rs.a3rs20_0 +
                                 s11_1 * a3rs.a3rs20_1 +
                                 s11_2 * a3rs.a3rs20_2;
            Bab_rsOut.Bab_rs20_1 = d2Heta_dr2 * a3s.a3r20 + d2Heta_ds2 * a3r.a3r02 +
                                 s22_0 * a3rs.a3rs20_0 +
                                 s22_1 * a3rs.a3rs20_1 +
                                 s22_2 * a3rs.a3rs20_2;
            Bab_rsOut.Bab_rs20_2 = (d2KsiHeta_dr2 * a3s.a3r20 + d2KsiHeta_ds2 * a3r.a3r02) * 2 +
                                 s12_0 * a3rs.a3rs20_0 +
                                 s12_1 * a3rs.a3rs20_1 +
                                 s12_2 * a3rs.a3rs20_2;

            Bab_rsOut.Bab_rs21_0 = d2Ksi_dr2 * a3s.a3r21 + d2Ksi_ds2 * a3r.a3r12 +
                                 s11_0 * a3rs.a3rs21_0 +
                                 s11_1 * a3rs.a3rs21_1 +
                                 s11_2 * a3rs.a3rs21_2;
            Bab_rsOut.Bab_rs21_1 = d2Heta_dr2 * a3s.a3r21 + d2Heta_ds2 * a3r.a3r12 +
                                 s22_0 * a3rs.a3rs21_0 +
                                 s22_1 * a3rs.a3rs21_1 +
                                 s22_2 * a3rs.a3rs21_2;
            Bab_rsOut.Bab_rs21_2 = (d2KsiHeta_dr2 * a3s.a3r21 + d2KsiHeta_ds2 * a3r.a3r12) * 2 +
                                 s12_0 * a3rs.a3rs21_0 +
                                 s12_1 * a3rs.a3rs21_1 +
                                 s12_2 * a3rs.a3rs21_2;

            Bab_rsOut.Bab_rs22_0 = d2Ksi_dr2 * a3s.a3r22 + d2Ksi_ds2 * a3r.a3r22 +
                                 s11_0 * a3rs.a3rs22_0 +
                                 s11_1 * a3rs.a3rs22_1 +
                                 s11_2 * a3rs.a3rs22_2;
            Bab_rsOut.Bab_rs22_1 = d2Heta_dr2 * a3s.a3r22 + d2Heta_ds2 * a3r.a3r22 +
                                 s22_0 * a3rs.a3rs22_0 +
                                 s22_1 * a3rs.a3rs22_1 +
                                 s22_2 * a3rs.a3rs22_2;
            Bab_rsOut.Bab_rs22_2 = (d2KsiHeta_dr2 * a3s.a3r22 + d2KsiHeta_ds2 * a3r.a3r22) * 2 +
                                 s12_0 * a3rs.a3rs22_0 +
                                 s12_1 * a3rs.a3rs22_1 +
                                 s12_2 * a3rs.a3rs22_2;
        }

        private static void Calculate_a3rs(double[] surfaceBasisVector1, double[] surfaceBasisVector2,
            double[] surfaceBasisVector3, double J1, double dksi_r, double dheta_r, double dksi_s, double dheta_s, ref Elements.a3rs a3rsOut)
        {
            #region Initializations
            var s10 = surfaceBasisVector1[0];
            var s11 = surfaceBasisVector1[1];
            var s12 = surfaceBasisVector1[2];

            var s20 = surfaceBasisVector2[0];
            var s21 = surfaceBasisVector2[1];
            var s22 = surfaceBasisVector2[2];

            var s30 = surfaceBasisVector3[0];
            var s31 = surfaceBasisVector3[1];
            var s32 = surfaceBasisVector3[2];

            var aux1Term1 = (dheta_s * dksi_r - dheta_r * dksi_s) * J1;
            var aux2Term1 = dheta_r * dksi_s - dheta_s * dksi_r;
            var aux3Term1 = aux2Term1 * J1;

            var aux1Term2 = dheta_s * s11 - dksi_s * s21;
            var aux2Term2 = dheta_s * s12 - dksi_s * s22;
            var aux3Term2 = dheta_r * s12 - dksi_r * s22;
            var aux4Term2 = dheta_r * s11 - dksi_r * s21;
            var aux5Term2 = dheta_s * s10 - dksi_s * s20;
            var aux6Term2 = dheta_r * s10 - dksi_r * s20;
            var aux7Term2 = s32 * aux1Term2 - s31 * aux2Term2;
            var aux8Term2 = s32 * aux5Term2 - s30 * aux2Term2;
            var aux9Term2 = s31 * aux5Term2 - s30 * aux1Term2;
            var J1squared = J1 * J1;

            var aux1Term3 = s32 * aux4Term2 - s31 * aux3Term2;
            var aux2Term3 = s32 * aux6Term2 - s30 * aux3Term2;
            var aux3Term3 = s31 * aux6Term2 - s30 * aux4Term2;

            var aux1 = aux4Term2 * aux1Term2;
            var aux2 = aux3Term2 * aux2Term2;
            var aux3 = aux1Term3 * aux7Term2;
            var aux4 = aux7Term2 * aux3Term2;
            var aux5 = aux1Term3 * aux2Term2;
            var aux6 = aux7Term2 * aux4Term2;
            var aux7 = aux1Term3 * aux1Term2;
            var aux8 = aux1 + aux2 - aux3;
            var aux9 = aux1Term3 * aux8Term2;
            var aux10 = aux4Term2 * aux5Term2;
            var aux11 = J1 * s32 * J1 * aux2Term1;
            var aux12 = aux8Term2 * aux4Term2;
            var aux13 = aux8Term2 * aux3Term2;
            var aux14 = aux1Term3 * aux5Term2;
            var aux15 = (aux10 + aux11 - aux9);
            var aux16 = aux9Term2 * aux1Term3;
            var aux17 = aux3Term2 * aux5Term2;
            var aux18 = J1 * s31 * J1 * aux2Term1;
            var aux19 = aux9Term2 * aux3Term2;
            var aux20 = aux9Term2 * aux4Term2;
            var aux21 = (aux17 - aux18 + aux16);
            var aux22 = aux2Term3 * aux7Term2;
            var aux23 = aux6Term2 * aux1Term2;
            var aux24 = aux2Term3 * aux2Term2;
            var aux25 = aux7Term2 * aux6Term2;
            var aux26 = aux2Term3 * aux1Term2;
            var aux27 = aux23 - aux11 - aux22;
            var aux28 = aux2Term3 * aux8Term2;
            var aux29 = aux6Term2 * aux5Term2;
            var aux30 = aux8Term2 * aux6Term2;
            var aux31 = aux2Term3 * aux5Term2;
            var aux32 = (aux29 + aux2 - aux28);
            var aux33 = aux2Term3 * aux9Term2;
            var aux34 = J1 * s30 * J1 * aux2Term1;
            var aux35 = aux3Term2 * aux1Term2;
            var aux36 = aux9Term2 * aux6Term2;
            var aux37 = (aux35 + aux34 - aux33);
            var aux38 = aux3Term3 * aux7Term2;
            var aux39 = aux6Term2 * aux2Term2;
            var aux40 = aux3Term3 * aux2Term2;
            var aux41 = aux3Term3 * aux1Term2;
            var aux42 = (aux39 + aux18 + aux38);
            var aux43 = aux3Term3 * aux8Term2;
            var aux44 = aux4Term2 * aux2Term2;
            var aux45 = aux3Term3 * aux5Term2;
            var aux46 = (aux44 - aux34 - aux43);
            var aux47 = aux3Term3 * aux9Term2;
            var aux48 = (aux29 + aux1 - aux47);
            #endregion

            #region Term1



            a3rsOut.a3rs00_0 = (2 * s30 * aux3 - s30 * aux8) / J1squared;
            a3rsOut.a3rs00_1 = (aux4 + aux5 + 2 * s31 * aux3 - s31 * aux8) / J1squared;
            a3rsOut.a3rs00_2 = -(aux6 + aux7 - 2 * s32 * aux3 + s32 * aux8) / J1squared;

            a3rsOut.a3rs01_0 = -(aux5 + 2 * s30 * aux9 - s30 * aux15) / J1squared;
            a3rsOut.a3rs01_1 = -(aux13 + 2 * s31 * aux9 - s31 * aux15) / J1squared;
            a3rsOut.a3rs01_2 = aux1Term1 + (aux12 + aux14 - 2 * s32 * aux9 + s32 * aux15) / J1squared;

            a3rsOut.a3rs02_0 = (aux7 + 2 * s30 * aux16 + s30 * aux21) / J1squared;
            a3rsOut.a3rs02_1 = aux3Term1 + (aux19 - aux14 + 2 * s31 * aux16 + s31 * aux21) / J1squared;
            a3rsOut.a3rs02_2 = -(aux20 - 2 * s32 * aux16 - s32 * aux21) / J1squared;

            a3rsOut.a3rs10_0 = -(aux4 + 2 * s30 * aux22 - s30 * aux27) / J1squared;
            a3rsOut.a3rs10_1 = -(aux24 + 2 * s31 * aux22 - s31 * aux27) / J1squared;
            a3rsOut.a3rs10_2 = aux3Term1 + (aux25 + aux26 - 2 * s32 * aux22 + s32 * aux27) / J1squared;

            a3rsOut.a3rs11_0 = (aux13 + aux24 + 2 * s30 * aux28 - s30 * aux32) / J1squared;
            a3rsOut.a3rs11_1 = (2 * s31 * aux28 - s31 * aux32) / J1squared;
            a3rsOut.a3rs11_2 = -(aux30 + aux31 - 2 * s32 * aux28 + s32 * aux32) / J1squared;

            a3rsOut.a3rs12_0 = aux1Term1 - (aux19 + aux26 + 2 * s30 * aux33 - s30 * aux37) / J1squared;
            a3rsOut.a3rs12_1 = (aux31 - 2 * s31 * aux33 + s31 * aux37) / J1squared;
            a3rsOut.a3rs12_2 = (aux36 - 2 * s32 * aux33 + s32 * aux37) / J1squared;

            a3rsOut.a3rs20_0 = (aux6 + 2 * s30 * aux38 + s30 * aux42) / J1squared;
            a3rsOut.a3rs20_1 = aux1Term1 - (aux25 - aux40 - 2 * s31 * aux38 - s31 * aux42) / J1squared;
            a3rsOut.a3rs20_2 = -(aux41 - 2 * s32 * aux38 - s32 * aux42) / J1squared;

            a3rsOut.a3rs21_0 = aux3Term1 - (aux12 + aux40 + 2 * s30 * aux43 - s30 * aux46) / J1squared;
            a3rsOut.a3rs21_1 = (aux30 - 2 * s31 * aux43 + s31 * aux46) / J1squared;
            a3rsOut.a3rs21_2 = (aux45 - 2 * s32 * aux43 + s32 * aux46) / J1squared;

            a3rsOut.a3rs22_0 = (aux20 + aux41 + 2 * s30 * aux47 - s30 * aux48) / J1squared;
            a3rsOut.a3rs22_1 = -(aux36 + aux45 - 2 * s31 * aux47 + s31 * aux48) / J1squared;
            a3rsOut.a3rs22_2 = (2 * s32 * aux47 - s32 * aux48) / J1squared;

            #endregion
        }


        private void CalculateA3r(double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double dksi_r, double dheta_r, double J1, ref Elements.a3r da3_unit_dr_out)
        {
            var s30 = surfaceBasisVector3[0];
            var s31 = surfaceBasisVector3[1];
            var s32 = surfaceBasisVector3[2];

            var da3_tilde_dr10 = dheta_r * surfaceBasisVector1[2] - dksi_r * surfaceBasisVector2[2];
            var da3_tilde_dr20 = dksi_r * surfaceBasisVector2[1] - dheta_r * surfaceBasisVector1[1];

            var da3_tilde_dr01 = dksi_r * surfaceBasisVector2[2] - dheta_r * surfaceBasisVector1[2];
            var da3_tilde_dr21 = dheta_r * surfaceBasisVector1[0] - dksi_r * surfaceBasisVector2[0];

            var da3_tilde_dr02 = dheta_r * surfaceBasisVector1[1] - dksi_r * surfaceBasisVector2[1];
            var da3_tilde_dr12 = dksi_r * surfaceBasisVector2[0] - dheta_r * surfaceBasisVector1[0];


            var dnorma3_dr0 = s31 * da3_tilde_dr10 +
                              s32 * da3_tilde_dr20;

            var dnorma3_dr1 = s30 * da3_tilde_dr01 +
                              s32 * da3_tilde_dr21;

            var dnorma3_dr2 = s30 * da3_tilde_dr02 +
                              s31 * da3_tilde_dr12;

            da3_unit_dr_out.a3r00 = -s30 * dnorma3_dr0;
            da3_unit_dr_out.a3r10 = da3_tilde_dr10 - s31 * dnorma3_dr0;
            da3_unit_dr_out.a3r20 = da3_tilde_dr20 - s32 * dnorma3_dr0;

            da3_unit_dr_out.a3r01 = da3_tilde_dr01 - s30 * dnorma3_dr1;
            da3_unit_dr_out.a3r11 = -s31 * dnorma3_dr1;
            da3_unit_dr_out.a3r21 = da3_tilde_dr21 - s32 * dnorma3_dr1;

            da3_unit_dr_out.a3r02 = da3_tilde_dr02 - s30 * dnorma3_dr2;
            da3_unit_dr_out.a3r12 = da3_tilde_dr12 - s31 * dnorma3_dr2;
            da3_unit_dr_out.a3r22 = -s32 * dnorma3_dr2;
        }

        internal void CalculateKmembraneNL(ControlPoint[] controlPoints, ref Elements.Forces membraneForces, IShapeFunction2D shapeFunctions,
            int j, double[,] KmembraneNLOut)
        {
            for (var i = 0; i < controlPoints.Length; i++)
            {
                var dksi_r = shapeFunctions.DerivativeValuesKsi[i, j];
                var dheta_r = shapeFunctions.DerivativeValuesHeta[i, j];

                for (int k = 0; k < controlPoints.Length; k++)
                {
                    var dksi_s = shapeFunctions.DerivativeValuesKsi[k, j];
                    var dheta_s = shapeFunctions.DerivativeValuesHeta[k, j];


                    var aux = membraneForces.v0 * dksi_r * dksi_s +
                              membraneForces.v1 * dheta_r * dheta_s +
                              membraneForces.v2 * (dksi_r * dheta_s + dksi_s * dheta_r);

                    KmembraneNLOut[i * 3, k * 3] += aux;
                    KmembraneNLOut[i * 3 + 1, k * 3 + 1] += aux;
                    KmembraneNLOut[i * 3 + 2, k * 3 + 2] += aux;
                }
            }
        }

        internal void CalculateMembraneDeformationMatrix(int controlPointsCount, IShapeFunction2D shapeFunctions, int j,
            double[] surfaceBasisVector1, double[] surfaceBasisVector2, double[,] BmembraneOut)
        {
            var s1_0 = surfaceBasisVector1[0];
            var s1_1 = surfaceBasisVector1[1];
            var s1_2 = surfaceBasisVector1[2];

            var s2_0 = surfaceBasisVector2[0];
            var s2_1 = surfaceBasisVector2[1];
            var s2_2 = surfaceBasisVector2[2];

            for (int column = 0; column < controlPointsCount * 3; column += 3)
            {
                var dKsi = shapeFunctions.DerivativeValuesKsi[column / 3, j];
                var dHeta = shapeFunctions.DerivativeValuesHeta[column / 3, j];

                BmembraneOut[0, column] = dKsi * s1_0;
                BmembraneOut[0, column + 1] = dKsi * s1_1;
                BmembraneOut[0, column + 2] = dKsi * s1_2;

                BmembraneOut[1, column] = dHeta * s2_0;
                BmembraneOut[1, column + 1] = dHeta * s2_1;
                BmembraneOut[1, column + 2] = dHeta * s2_2;

                BmembraneOut[2, column] = dHeta * s1_0 + dKsi * s2_0;
                BmembraneOut[2, column + 1] = dHeta * s1_1 + dKsi * s2_1;
                BmembraneOut[2, column + 2] = dHeta * s1_2 + dKsi * s2_2;
            }
        }

        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(KirchhoffLoveShellNL shellElement)
        {
            var gauss = new GaussQuadrature();
            var medianSurfaceGP = gauss.CalculateElementGaussPoints(_degreeKsi,
                _degreeHeta, shellElement.Knots.ToList());
            foreach (var point in medianSurfaceGP)
            {
                var gp = gauss.CalculateElementGaussPoints(ThicknessIntegrationDegree,
                    new List<Knot>
                    {
                        new Knot() {ID = 0, Ksi = -shellElement.Thickness / 2, Heta = point.Heta},
                        new Knot() {ID = 1, Ksi = shellElement.Thickness / 2, Heta = point.Heta},
                    }).ToList();

                thicknessIntegrationPoints.Add(point,
                    gp.Select(g => new GaussLegendrePoint3D(point.Ksi, point.Heta, g.Ksi, g.WeightFactor))
                        .ToList());
            }

            return medianSurfaceGP;
        }

        private const int ThicknessIntegrationDegree = 2;

        public double Thickness { get; set; }

        private readonly int _degreeKsi;
        private readonly int _degreeHeta;
    }

    public struct a3r
    {
        public double a3r00;
        public double a3r01;
        public double a3r02;

        public double a3r10;
        public double a3r11;
        public double a3r12;

        public double a3r20;
        public double a3r21;
        public double a3r22;
    }

    public struct a3rs
    {
        public double a3rs00_0;
        public double a3rs00_1;
        public double a3rs00_2;

        public double a3rs01_0;
        public double a3rs01_1;
        public double a3rs01_2;

        public double a3rs02_0;
        public double a3rs02_1;
        public double a3rs02_2;

        public double a3rs10_0;
        public double a3rs10_1;
        public double a3rs10_2;

        public double a3rs11_0;
        public double a3rs11_1;
        public double a3rs11_2;

        public double a3rs12_0;
        public double a3rs12_1;
        public double a3rs12_2;

        public double a3rs20_0;
        public double a3rs20_1;
        public double a3rs20_2;

        public double a3rs21_0;
        public double a3rs21_1;
        public double a3rs21_2;

        public double a3rs22_0;
        public double a3rs22_1;
        public double a3rs22_2;
    }

    public struct Bab_rs
    {
        public double Bab_rs00_0;
        public double Bab_rs00_1;
        public double Bab_rs00_2;

        public double Bab_rs01_0;
        public double Bab_rs01_1;
        public double Bab_rs01_2;

        public double Bab_rs02_0;
        public double Bab_rs02_1;
        public double Bab_rs02_2;

        public double Bab_rs10_0;
        public double Bab_rs10_1;
        public double Bab_rs10_2;

        public double Bab_rs11_0;
        public double Bab_rs11_1;
        public double Bab_rs11_2;

        public double Bab_rs12_0;
        public double Bab_rs12_1;
        public double Bab_rs12_2;

        public double Bab_rs20_0;
        public double Bab_rs20_1;
        public double Bab_rs20_2;

        public double Bab_rs21_0;
        public double Bab_rs21_1;
        public double Bab_rs21_2;

        public double Bab_rs22_0;
        public double Bab_rs22_1;
        public double Bab_rs22_2;
    }

    public struct Forces
    {
        public double v0;
        public double v1;
        public double v2;
    }
}
