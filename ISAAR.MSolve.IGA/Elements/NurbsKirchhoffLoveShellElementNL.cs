using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using System.Runtime.CompilerServices;
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

[assembly: InternalsVisibleTo("ISAAR.MSolve.IGA.Tests")]

namespace ISAAR.MSolve.IGA.Elements
{
    public class NurbsKirchhoffLoveShellElementNL : Element, IStructuralIsogeometricElement, ISurfaceLoadedElement
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

        public NurbsKirchhoffLoveShellElementNL(IShellMaterial shellMaterial, IList<Knot> elementKnots,
            IList<ControlPoint> elementControlPoints, Patch patch, double thickness)
        {
            Contract.Requires(shellMaterial != null);
            this.Patch = patch;
            this.Thickness = thickness;
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

            CreateElementGaussPoints(this);
            foreach (var medianSurfaceGP in thicknessIntegrationPoints.Keys)
            {
                materialsAtThicknessGP.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IShellMaterial>());
                foreach (var point in thicknessIntegrationPoints[medianSurfaceGP])
                {
                    materialsAtThicknessGP[medianSurfaceGP].Add(point, shellMaterial.Clone());
                }
            }
        }

        public CellType CellType { get; } = CellType.Unknown;

        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public bool MaterialModified => false;

        public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) =>
            throw new NotImplementedException();

        public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
        {
            var nurbsElement = (NurbsKirchhoffLoveShellElementNL) element;
            var knotParametricCoordinatesKsi = Vector.CreateFromArray(Knots.Select(k => k.Ksi).ToArray());
            var knotParametricCoordinatesHeta = Vector.CreateFromArray(Knots.Select(k => k.Heta).ToArray());

            var nurbs = new Nurbs2D(nurbsElement, nurbsElement.ControlPoints.ToArray(), knotParametricCoordinatesKsi,
                knotParametricCoordinatesHeta);

            var knotDisplacements = new double[4, 3];
            var paraviewKnotRenumbering = new int[] {0, 3, 1, 2};
            for (var j = 0; j < knotDisplacements.GetLength(0); j++)
            {
                for (int i = 0; i < element.ControlPoints.Count(); i++)
                {
                    knotDisplacements[paraviewKnotRenumbering[j], 0] +=
                        nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
                    knotDisplacements[paraviewKnotRenumbering[j], 1] +=
                        nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
                    knotDisplacements[paraviewKnotRenumbering[j], 2] +=
                        nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
                }
            }

            return knotDisplacements;
        }

        public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNL) element;
            var controlPoints = shellElement.ControlPoints.ToArray();
            var elementNodalForces = new double[shellElement.ControlPointsDictionary.Count * 3];

            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(controlPoints);
            var nurbs = CalculateShapeFunctions(shellElement, shellElement.ControlPoints);
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();

            for (int j = 0; j < gaussPoints.Length; j++)
            {
                var jacobianMatrix = CalculateJacobian(newControlPoints, nurbs, j);

                var hessianMatrix = CalculateHessian(newControlPoints, nurbs, j);

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

                var Bmembrane = CalculateMembraneDeformationMatrix(newControlPoints, nurbs, j, surfaceBasisVector1,
                    surfaceBasisVector2);
                var Bbending = CalculateBendingDeformationMatrix(newControlPoints, surfaceBasisVector3, nurbs, j,
                    surfaceBasisVector2,
                    surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12);

                var (membraneForces, bendingMoments) =
                    IntegratedStressesOverThickness(gaussPoints[j]);

                var wfactor = InitialJ1[j] * gaussPoints[j].WeightFactor;


                for (int i = 0; i < Bmembrane.GetLength(1); i++)
                {
                    for (int k = 0; k < Bmembrane.GetLength(0); k++)
                    {
                        elementNodalForces[i] += (Bmembrane[k, i] * membraneForces[k] * wfactor +
                                                  Bbending[k, i] * bendingMoments[k] * wfactor);
                    }
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

        public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
            double[] localdDisplacements)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNL) element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();
            var nurbs = CalculateShapeFunctions(shellElement, shellElement.ControlPoints);

            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(elementControlPoints);


            //var newControlPoints = elementControlPoints;

            var midsurfaceGP = materialsAtThicknessGP.Keys.ToArray();
            for (var j = 0; j < midsurfaceGP.Length; j++)
            {
                var jacobianMatrix = CalculateJacobian(newControlPoints, nurbs, j);

                var hessianMatrix = CalculateHessian(newControlPoints, nurbs, j);

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

                var membraneStrain = new double[] {0.5 * (a11 - A11), 0.5 * (a22 - A22), a12 - A12};

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

                var bendingStrain = new double[] {b11 - B11, b22 - B22, 2 * b12 - 2 * B12};

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
            var shellElement = (NurbsKirchhoffLoveShellElementNL) element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();
            var gaussPoints = CreateElementGaussPoints(shellElement);
            var distributedLoad = new Dictionary<int, double>();
            var nurbs = new Nurbs2D(shellElement, elementControlPoints);

            for (var j = 0; j < gaussPoints.Count; j++)
            {
                var jacobianMatrix = CalculateJacobian(elementControlPoints, nurbs, j);
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
                                                  nurbs.NurbsValues[i, j] * gaussPoints[j].WeightFactor;
                    }
                    else
                    {
                        distributedLoad.Add(dofId,
                            loadMagnitude * nurbs.NurbsValues[i, j] * J1 * gaussPoints[j].WeightFactor);
                    }
                }
            }

            return distributedLoad;
        }

        public Dictionary<int, double> CalculateSurfacePressure(Element element, double pressureMagnitude)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNL) element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();
            var gaussPoints = CreateElementGaussPoints(shellElement);
            var pressureLoad = new Dictionary<int, double>();
            var nurbs = new Nurbs2D(shellElement, elementControlPoints);

            for (var j = 0; j < gaussPoints.Count; j++)
            {
                var jacobianMatrix = CalculateJacobian(elementControlPoints, nurbs, j);
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
                                                   nurbs.NurbsValues[i, j] * gaussPoints[j].WeightFactor;
                        }
                        else
                        {
                            pressureLoad.Add(dofId,
                                pressureMagnitude * surfaceBasisVector3[k] * nurbs.NurbsValues[i, j] *
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

        internal (double[] MembraneForces, double[] BendingMoments) IntegratedStressesOverThickness(
            GaussLegendrePoint3D midSurfaceGaussPoint)
        {
            var MembraneForces = new double[3];
            var BendingMoments = new double[3];
            var thicknessPoints = thicknessIntegrationPoints[midSurfaceGaussPoint];

            for (int i = 0; i < thicknessPoints.Count; i++)
            {
                var thicknessPoint = thicknessPoints[i];
                var material = materialsAtThicknessGP[midSurfaceGaussPoint][thicknessPoints[i]];
                var w = thicknessPoint.WeightFactor;
                var z = thicknessPoint.Zeta;
                for (int j = 0; j < 3; j++)
                {
                    MembraneForces[j] += material.Stresses[j] * w;
                    BendingMoments[j] -= material.Stresses[j] * w * z;
                }
            }

            return (MembraneForces, BendingMoments);
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

        private Nurbs2D _nurbs;

        public IMatrix StiffnessMatrix(IElement element)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNL) element;
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();

            var controlPoints = shellElement.ControlPoints.ToArray();
            var nurbs = CalculateShapeFunctions(shellElement, shellElement.ControlPoints);

            if (!isInitialized)
            {
                CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
                isInitialized = true;
            }

            var elementControlPoints = CurrentControlPoint(controlPoints);

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
                var jacobianMatrix = CalculateJacobian(elementControlPoints, nurbs, j);

                var hessianMatrix = CalculateHessian(elementControlPoints, nurbs, j);
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

                var Bmembrane = CalculateMembraneDeformationMatrix(elementControlPoints, nurbs, j, surfaceBasisVector1,
                    surfaceBasisVector2);
                var Bbending = CalculateBendingDeformationMatrix(elementControlPoints, surfaceBasisVector3, nurbs, j,
                    surfaceBasisVector2,
                    surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12);

                var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
                    IntegratedConstitutiveOverThickness(gaussPoints[j]);

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

                var (MembraneForces, BendingMoments) = IntegratedStressesOverThickness(gaussPoints[j]);

                var KmembraneNL = CalculateKmembraneNL(elementControlPoints, MembraneForces, nurbs, j);
                var KbendingNL = CalculateKbendingNL(elementControlPoints, BendingMoments, nurbs,
                    surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3,
                    surfaceBasisVectorDerivative1,
                    surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12, J1, j);

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

        private Nurbs2D CalculateShapeFunctions(NurbsKirchhoffLoveShellElementNL shellElement,
            IEnumerable<ControlPoint> controlPoints)
        {
            return _nurbs ?? (_nurbs = new Nurbs2D(shellElement, controlPoints.ToArray()));
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

        internal double[,] CalculateHessian(ControlPoint[] controlPoints, Nurbs2D nurbs, int j)
        {
            var hessianMatrix = new double[3, 3];
            for (var k = 0; k < controlPoints.Length; k++)
            {
                hessianMatrix[0, 0] +=
                    nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].X;
                hessianMatrix[0, 1] +=
                    nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].Y;
                hessianMatrix[0, 2] +=
                    nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].Z;
                hessianMatrix[1, 0] +=
                    nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].X;
                hessianMatrix[1, 1] +=
                    nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].Y;
                hessianMatrix[1, 2] +=
                    nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].Z;
                hessianMatrix[2, 0] +=
                    nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].X;
                hessianMatrix[2, 1] +=
                    nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].Y;
                hessianMatrix[2, 2] +=
                    nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].Z;
            }

            return hessianMatrix;
        }

        internal double[,] CalculateJacobian(ControlPoint[] controlPoints, Nurbs2D nurbs, int j)
        {
            var jacobianMatrix = new double[2, 3];
            for (var k = 0; k < controlPoints.Length; k++)
            {
                jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * controlPoints[k].X;
                jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * controlPoints[k].Y;
                jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * controlPoints[k].Z;
                jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * controlPoints[k].X;
                jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * controlPoints[k].Y;
                jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * controlPoints[k].Z;
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
            Nurbs2D nurbs, int j, double[] surfaceBasisVector2, double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVector1,
            double J1, double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12)
        {
            var Bbending = new double[3, controlPoints.Length * 3];
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

            //var s1 = Vector.CreateFromArray(surfaceBasisVector1);
            //var s2 = Vector.CreateFromArray(surfaceBasisVector2);
            //var s3 = Vector.CreateFromArray(surfaceBasisVector3);
            //var s11 = Vector.CreateFromArray(surfaceBasisVectorDerivative1);
            //var s22 = Vector.CreateFromArray(surfaceBasisVectorDerivative2);
            //var s12 = Vector.CreateFromArray(surfaceBasisVectorDerivative12);
            for (int column = 0; column < controlPoints.Length * 3; column += 3)
            {
                var dksi = nurbs.NurbsDerivativeValuesKsi[column / 3, j];
                var dheta = nurbs.NurbsDerivativeValuesHeta[column / 3, j];
                var d2Ksi = nurbs.NurbsSecondDerivativeValueKsi[column / 3, j];
                var d2Heta = nurbs.NurbsSecondDerivativeValueHeta[column / 3, j];
                var d2KsiHeta = nurbs.NurbsSecondDerivativeValueKsiHeta[column / 3, j];

                //#region BI1

                //var BI1 = s3.CrossProduct(s1);
                //BI1.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);

                //var auxVector = s2.CrossProduct(s3);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                //BI1.AddIntoThis(auxVector);

                //BI1.ScaleIntoThis(s3.DotProduct(s11));
                //auxVector = s1.CrossProduct(s11);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                //BI1.AddIntoThis(auxVector);

                //auxVector = s11.CrossProduct(s2);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                //BI1.AddIntoThis(auxVector);

                //BI1.ScaleIntoThis(1 / J1);
                //auxVector[0] = surfaceBasisVector3[0];
                //auxVector[1] = surfaceBasisVector3[1];
                //auxVector[2] = surfaceBasisVector3[2];
                //auxVector.ScaleIntoThis(-nurbs.NurbsSecondDerivativeValueKsi[column / 3, j]);
                //BI1.AddIntoThis(auxVector);

                //#endregion BI1

                //#region BI2

                //IVector BI2 = s3.CrossProduct(s1);
                //BI2.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                //auxVector = s2.CrossProduct(s3);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                //BI2.AddIntoThis(auxVector);
                //BI2.ScaleIntoThis(s3.DotProduct(s22));
                //auxVector = s1.CrossProduct(s22);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                //BI2.AddIntoThis(auxVector);
                //auxVector = s22.CrossProduct(s2);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                //BI2.AddIntoThis(auxVector);
                //BI2.ScaleIntoThis(1 / J1);
                //auxVector[0] = surfaceBasisVector3[0];
                //auxVector[1] = surfaceBasisVector3[1];
                //auxVector[2] = surfaceBasisVector3[2];
                //auxVector.ScaleIntoThis(-nurbs.NurbsSecondDerivativeValueHeta[column / 3, j]);
                //BI2.AddIntoThis(auxVector);

                //#endregion BI2

                //#region BI3

                //var BI3 = s3.CrossProduct(s1);
                //BI3.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                //auxVector = s2.CrossProduct(s3);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                //BI3.AddIntoThis(auxVector);
                //BI3.ScaleIntoThis(s3.DotProduct(s12));
                //auxVector = s1.CrossProduct(s12);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                //BI3.AddIntoThis(auxVector);
                //auxVector = s22.CrossProduct(s2);
                //auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                //BI3.AddIntoThis(auxVector);
                //BI3.ScaleIntoThis(1 / J1);
                //auxVector[0] = surfaceBasisVector3[0];
                //auxVector[1] = surfaceBasisVector3[1];
                //auxVector[2] = surfaceBasisVector3[2];
                //auxVector.ScaleIntoThis(-nurbs.NurbsSecondDerivativeValueKsiHeta[column / 3, j]);
                //BI3.AddIntoThis(auxVector);

                //#endregion BI3


                //Bbending[0, column] = BI1[0];
                //Bbending[0, column + 1] = BI1[1];
                //Bbending[0, column + 2] = BI1[2];

                //Bbending[1, column] = BI2[0];
                //Bbending[1, column + 1] = BI2[1];
                //Bbending[1, column + 2] = BI2[2];

                //Bbending[2, column] = 2 * BI3[0];
                //Bbending[2, column + 1] = 2 * BI3[1];
                //Bbending[2, column + 2] = 2 * BI3[2];



                Bbending[0, column] = -d2Ksi * s30 - ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s11 * s11_2 - s12 * s11_1) + dksi * (s21 * s11_2 - s22 * s11_1)) / J1;
                Bbending[0, column + 1] = ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s10 * s11_2 - s12 * s11_0) + dksi * (s20 * s11_2 - s22 * s11_0)) / J1 - d2Ksi * s31;
                Bbending[0, column + 2] = -d2Ksi * s32 - ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s10 * s11_1 - s11 * s11_0) + dksi * (s20 * s11_1 - s21 * s11_0)) / J1;

                Bbending[1, column] = -d2Heta * s30 - ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s11 * s22_2 - s12 * s22_1) + dksi * (s21 * s22_2 - s22 * s22_1)) / J1;
                Bbending[1, column + 1] = ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s10 * s22_2 - s12 * s22_0) + dksi * (s20 * s22_2 - s22 * s22_0)) / J1 - d2Heta * s31;
                Bbending[1, column + 2] = -d2Heta * s32 - ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s10 * s22_1 - s11 * s22_0) + dksi * (s20 * s22_1 - s21 * s22_0)) / J1;

                Bbending[2, column] = -2 * d2KsiHeta * s30 - (2 * ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s11 * s12_2 - s12 * s12_1) + dksi * (s21 * s12_2 - s22 * s12_1))) / J1;
                Bbending[2, column + 1] = (2 * ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s10 * s12_2 - s12 * s12_0) + dksi * (s20 * s12_2 - s22 * s12_0))) / J1 - 2 * d2KsiHeta * s31;
                Bbending[2, column + 2] = -2 * d2KsiHeta * s32 - (2 * ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s10 * s12_1 - s11 * s12_0) + dksi * (s20 * s12_1 - s21 * s12_0))) / J1;
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

        internal void CalculateInitialConfigurationData(ControlPoint[] controlPoints,
            Nurbs2D nurbs, IList<GaussLegendrePoint3D> gaussPoints)
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
                var jacobianMatrix = CalculateJacobian(controlPoints, nurbs, j);

                var hessianMatrix = CalculateHessian(controlPoints, nurbs, j);
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

        internal double[,] CalculateKbendingNL(ControlPoint[] controlPoints,
            double[] bendingMoments, Nurbs2D nurbs, double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double[] surfaceBasisVectorDerivative1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double J1, int j)
        {
            var KbendingNL = new double[controlPoints.Length * 3, controlPoints.Length * 3];

            for (int i = 0; i < controlPoints.Length; i++)
            {
                var dksi_r = nurbs.NurbsDerivativeValuesKsi[i, j];
                var dheta_r = nurbs.NurbsDerivativeValuesHeta[i, j];

                var d2Ksi_dr2 = nurbs.NurbsSecondDerivativeValueKsi[i, j];
                var d2Heta_dr2 = nurbs.NurbsSecondDerivativeValueHeta[i, j];
                var d2KsiHeta_dr2 = nurbs.NurbsSecondDerivativeValueKsiHeta[i, j];

                for (int k = 0; k < controlPoints.Length; k++)
                {
                    var d2Ksi_ds2 = nurbs.NurbsSecondDerivativeValueKsi[k, j];
                    var d2Heta_ds2 = nurbs.NurbsSecondDerivativeValueHeta[k, j];
                    var d2KsiHeta_ds2 = nurbs.NurbsSecondDerivativeValueKsiHeta[k, j];

                    var dksi_s = nurbs.NurbsDerivativeValuesKsi[k, j];
                    var dheta_s = nurbs.NurbsDerivativeValuesHeta[k, j];

                    var a3r = CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r,
                        dheta_r, J1);
                    var a3s = CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_s,
                        dheta_s, J1);

                    var a3rs = Calculate_a3rs(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, J1, dksi_r,
                        dheta_r, dksi_s, dheta_s);
                    
                    var Bab_rs = CalculateBab_rs(surfaceBasisVectorDerivative1, surfaceBasisVectorDerivative2, 
                        surfaceBasisVectorDerivative12, d2Ksi_dr2, a3s, d2Ksi_ds2, a3r, a3rs, d2Heta_dr2, 
                        d2Heta_ds2, d2KsiHeta_dr2, d2KsiHeta_ds2);

                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            KbendingNL[i * 3 + l, k * 3 + m] -= (Bab_rs[l, m][0] * bendingMoments[0] +
                                                                 Bab_rs[l, m][1] * bendingMoments[1] +
                                                                 Bab_rs[l, m][2] * bendingMoments[2]);
                        }
                    }
                }
            }

            return KbendingNL;
        }

        private static double[,][] CalculateBab_rs(double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12, double d2Ksi_dr2, double[,] a3s,
            double d2Ksi_ds2, double[,] a3r, double[,][] a3rs, double d2Heta_dr2, double d2Heta_ds2, double d2KsiHeta_dr2,
            double d2KsiHeta_ds2)
        {
            var Bab_rs = new double[3, 3][];

            Bab_rs[0, 0] = new double[3]
            {
                d2Ksi_dr2 * a3s[0, 0] + d2Ksi_ds2 * a3r[0, 0] +
                surfaceBasisVectorDerivative1[0] * a3rs[0, 0][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[0, 0][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[0, 0][2],
                d2Heta_dr2 * a3s[0, 0] + d2Heta_ds2 * a3r[0, 0] +
                surfaceBasisVectorDerivative2[0] * a3rs[0, 0][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[0, 0][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[0, 0][2],
                (d2KsiHeta_dr2 * a3s[0, 0] + d2KsiHeta_ds2 * a3r[0, 0]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[0, 0][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[0, 0][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[0, 0][2]
            };
            Bab_rs[0, 1] = new double[3]
            {
                d2Ksi_dr2 * a3s[0, 1] + d2Ksi_ds2 * a3r[1, 0] +
                surfaceBasisVectorDerivative1[0] * a3rs[0, 1][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[0, 1][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[0, 1][2],
                d2Heta_dr2 * a3s[0, 1] + d2Heta_ds2 * a3r[1, 0] +
                surfaceBasisVectorDerivative2[0] * a3rs[0, 1][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[0, 1][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[0, 1][2],
                (d2KsiHeta_dr2 * a3s[0, 1] + d2KsiHeta_ds2 * a3r[1, 0]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[0, 1][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[0, 1][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[0, 1][2]
            };
            Bab_rs[0, 2] = new double[3]
            {
                d2Ksi_dr2 * a3s[0, 2] + d2Ksi_ds2 * a3r[2, 0] +
                surfaceBasisVectorDerivative1[0] * a3rs[0, 2][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[0, 2][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[0, 2][2],
                d2Heta_dr2 * a3s[0, 2] + d2Heta_ds2 * a3r[2, 0] +
                surfaceBasisVectorDerivative2[0] * a3rs[0, 2][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[0, 2][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[0, 2][2],
                (d2KsiHeta_dr2 * a3s[0, 2] + d2KsiHeta_ds2 * a3r[2, 0]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[0, 2][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[0, 2][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[0, 2][2]
            };

            Bab_rs[1, 0] = new double[3]
            {
                d2Ksi_dr2 * a3s[1, 0] + d2Ksi_ds2 * a3r[0, 1] +
                surfaceBasisVectorDerivative1[0] * a3rs[1, 0][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[1, 0][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[1, 0][2],
                d2Heta_dr2 * a3s[1, 0] + d2Heta_ds2 * a3r[0, 1] +
                surfaceBasisVectorDerivative2[0] * a3rs[1, 0][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[1, 0][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[1, 0][2],
                (d2KsiHeta_dr2 * a3s[1, 0] + d2KsiHeta_ds2 * a3r[0, 1]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[1, 0][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[1, 0][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[1, 0][2]
            };
            Bab_rs[1, 1] = new double[3]
            {
                d2Ksi_dr2 * a3s[1, 1] + d2Ksi_ds2 * a3r[1, 1] +
                surfaceBasisVectorDerivative1[0] * a3rs[1, 1][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[1, 1][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[1, 1][2],
                d2Heta_dr2 * a3s[1, 1] + d2Heta_ds2 * a3r[1, 1] +
                surfaceBasisVectorDerivative2[0] * a3rs[1, 1][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[1, 1][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[1, 1][2],
                (d2KsiHeta_dr2 * a3s[1, 1] + d2KsiHeta_ds2 * a3r[1, 1]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[1, 1][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[1, 1][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[1, 1][2]
            };
            Bab_rs[1, 2] = new double[3]
            {
                d2Ksi_dr2 * a3s[1, 2] + d2Ksi_ds2 * a3r[2, 1] +
                surfaceBasisVectorDerivative1[0] * a3rs[1, 2][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[1, 2][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[1, 2][2],
                d2Heta_dr2 * a3s[1, 2] + d2Heta_ds2 * a3r[2, 1] +
                surfaceBasisVectorDerivative2[0] * a3rs[1, 2][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[1, 2][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[1, 2][2],
                (d2KsiHeta_dr2 * a3s[1, 2] + d2KsiHeta_ds2 * a3r[2, 1]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[1, 2][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[1, 2][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[1, 2][2]
            };

            Bab_rs[2, 0] = new double[3]
            {
                d2Ksi_dr2 * a3s[2, 0] + d2Ksi_ds2 * a3r[0, 2] +
                surfaceBasisVectorDerivative1[0] * a3rs[2, 0][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[2, 0][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[2, 0][2],
                d2Heta_dr2 * a3s[2, 0] + d2Heta_ds2 * a3r[0, 2] +
                surfaceBasisVectorDerivative2[0] * a3rs[2, 0][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[2, 0][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[2, 0][2],
                (d2KsiHeta_dr2 * a3s[2, 0] + d2KsiHeta_ds2 * a3r[0, 2]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[2, 0][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[2, 0][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[2, 0][2]
            };
            Bab_rs[2, 1] = new double[3]
            {
                d2Ksi_dr2 * a3s[2, 1] + d2Ksi_ds2 * a3r[1, 2] +
                surfaceBasisVectorDerivative1[0] * a3rs[2, 1][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[2, 1][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[2, 1][2],
                d2Heta_dr2 * a3s[2, 1] + d2Heta_ds2 * a3r[1, 2] +
                surfaceBasisVectorDerivative2[0] * a3rs[2, 1][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[2, 1][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[2, 1][2],
                (d2KsiHeta_dr2 * a3s[2, 1] + d2KsiHeta_ds2 * a3r[1, 2]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[2, 1][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[2, 1][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[2, 1][2]
            };
            Bab_rs[2, 2] = new double[3]
            {
                d2Ksi_dr2 * a3s[2, 2] + d2Ksi_ds2 * a3r[2, 2] +
                surfaceBasisVectorDerivative1[0] * a3rs[2, 2][0] +
                surfaceBasisVectorDerivative1[1] * a3rs[2, 2][1] +
                surfaceBasisVectorDerivative1[2] * a3rs[2, 2][2],
                d2Heta_dr2 * a3s[2, 2] + d2Heta_ds2 * a3r[2, 2] +
                surfaceBasisVectorDerivative2[0] * a3rs[2, 2][0] +
                surfaceBasisVectorDerivative2[1] * a3rs[2, 2][1] +
                surfaceBasisVectorDerivative2[2] * a3rs[2, 2][2],
                (d2KsiHeta_dr2 * a3s[2, 2] + d2KsiHeta_ds2 * a3r[2, 2]) * 2 +
                surfaceBasisVectorDerivative12[0] * a3rs[2, 2][0] +
                surfaceBasisVectorDerivative12[1] * a3rs[2, 2][1] +
                surfaceBasisVectorDerivative12[2] * a3rs[2, 2][2]
            };
            return Bab_rs;
        }

        private static double[,][] Calculate_a3rs(double[] surfaceBasisVector1, double[] surfaceBasisVector2,
            double[] surfaceBasisVector3, double J1, double dksi_r, double dheta_r, double dksi_s, double dheta_s)
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

            var a3rs = new double[3, 3][];
            a3rs[0, 0] = new double[3];
            a3rs[0, 1] = new double[3];
            a3rs[0, 2] = new double[3];

            a3rs[1, 0] = new double[3];
            a3rs[1, 1] = new double[3];
            a3rs[1, 2] = new double[3];

            a3rs[2, 0] = new double[3];
            a3rs[2, 1] = new double[3];
            a3rs[2, 2] = new double[3];
            #endregion
            
            #region Term1
            var aux1Term1 = (dheta_s * dksi_r - dheta_r * dksi_s) * J1;
            var aux2Term1 = dheta_r * dksi_s - dheta_s * dksi_r;
            var aux3Term1 = aux2Term1 * J1;
            a3rs[0, 1][2] = aux1Term1;
            a3rs[0, 2][1] = aux3Term1;

            a3rs[1, 0][2] = aux3Term1;
            a3rs[1, 2][0] = aux1Term1;

            a3rs[2, 0][1] = aux1Term1;
            a3rs[2, 1][0] = aux3Term1;
            #endregion

            #region Term2
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


            a3rs[0, 0][1] += (aux7Term2 * aux3Term2) / J1squared;
            a3rs[0, 0][2] += -(aux7Term2 * aux4Term2) / J1squared;

            a3rs[0, 1][1] += -(aux8Term2 * aux3Term2) / J1squared;
            a3rs[0, 1][2] += (aux8Term2 * aux4Term2) / J1squared;

            a3rs[0, 2][1] += (aux9Term2 * aux3Term2) / J1squared;
            a3rs[0, 2][2] += -(aux9Term2 * aux4Term2) / J1squared;

            a3rs[1, 0][0] += -(aux7Term2 * aux3Term2) / J1squared;
            a3rs[1, 0][2] += (aux7Term2 * aux6Term2) / J1squared;

            a3rs[1, 1][0] += (aux8Term2 * aux3Term2) / J1squared;
            a3rs[1, 1][2] += -(aux8Term2 * aux6Term2) / J1squared;

            a3rs[1, 2][0] += -(aux9Term2 * aux3Term2) / J1squared;
            a3rs[1, 2][2] += (aux9Term2 * aux6Term2) / J1squared;

            a3rs[2, 0][0] += (aux7Term2 * aux4Term2) / J1squared;
            a3rs[2, 0][1] += -(aux7Term2 * aux6Term2) / J1squared;

            a3rs[2, 1][0] += -(aux8Term2 * aux4Term2) / J1squared;
            a3rs[2, 1][1] += (aux8Term2 * aux6Term2) / J1squared;

            a3rs[2, 2][0] += (aux9Term2 * aux4Term2) / J1squared;
            a3rs[2, 2][1] += -(aux9Term2 * aux6Term2) / J1squared;
            #endregion

            #region Term3
            var aux1Term3 = s32 * aux4Term2 - s31 * aux3Term2;
            var aux2Term3 = s32 * aux6Term2 - s30 * aux3Term2;
            var aux3Term3 = s31 * aux6Term2 - s30 * aux4Term2;

            a3rs[0, 0][1] += (aux1Term3 * aux2Term2) / J1squared;
            a3rs[0, 0][2] += -(aux1Term3 * aux1Term2) / J1squared;

            a3rs[0, 1][0] += -(aux1Term3 * aux2Term2) / J1squared;
            a3rs[0, 1][2] += (aux1Term3 * aux5Term2) / J1squared;

            a3rs[0, 2][0] += (aux1Term3 * aux1Term2) / J1squared;
            a3rs[0, 2][1] += -(aux1Term3 * aux5Term2) / J1squared;

            a3rs[1, 0][1] += -(aux2Term3 * aux2Term2) / J1squared;
            a3rs[1, 0][2] += (aux2Term3 * aux1Term2) / J1squared;

            a3rs[1, 1][0] += (aux2Term3 * aux2Term2) / J1squared;
            a3rs[1, 1][2] += -(aux2Term3 * aux5Term2) / J1squared;

            a3rs[1, 2][0] += -(aux2Term3 * aux1Term2) / J1squared;
            a3rs[1, 2][1] += (aux2Term3 * aux5Term2) / J1squared;

            a3rs[2, 0][1] += (aux3Term3 * aux2Term2) / J1squared;
            a3rs[2, 0][2] += -(aux3Term3 * aux1Term2) / J1squared;

            a3rs[2, 1][0] += -(aux3Term3 * aux2Term2) / J1squared;
            a3rs[2, 1][2] += (aux3Term3 * aux5Term2) / J1squared;

            a3rs[2, 2][0] += (aux3Term3 * aux1Term2) / J1squared;
            a3rs[2, 2][1] += -(aux3Term3 * aux5Term2) / J1squared;
            #endregion

            #region Term4
            a3rs[0, 0][0] += -(s30 * ((aux4Term2 * aux1Term2 + aux3Term2 * aux2Term2) / J1 -
                                      (aux1Term3 * aux7Term2) / J1)) / J1;
            a3rs[0, 0][1] += -(s31 * ((aux4Term2 * aux1Term2 + aux3Term2 * aux2Term2) / J1 -
                                      (aux1Term3 * aux7Term2) / J1)) / J1;
            a3rs[0, 0][2] += -(s32 * ((aux4Term2 * aux1Term2 + aux3Term2 * aux2Term2) / J1 -
                                      (aux1Term3 * aux7Term2) / J1)) / J1;

            a3rs[0, 1][0] += (s30 * ((aux4Term2 * aux5Term2 + J1 * s32 * J1 * aux2Term1) / J1 -
                                     (aux1Term3 * aux8Term2) / J1)) / J1;
            a3rs[0, 1][1] += (s31 * ((aux4Term2 * aux5Term2 + J1 * s32 * J1 * aux2Term1) / J1 -
                                     (aux1Term3 * aux8Term2) / J1)) / J1;
            a3rs[0, 1][2] += (s32 * ((aux4Term2 * aux5Term2 + J1 * s32 * J1 * aux2Term1) / J1 -
                                     (aux1Term3 * aux8Term2) / J1)) / J1;

            a3rs[0, 2][0] += (s30 * ((aux3Term2 * aux5Term2 - J1 * s31 * J1 * aux2Term1) / J1 +
                                     (aux9Term2 * aux1Term3) / J1)) / J1;
            a3rs[0, 2][1] += (s31 * ((aux3Term2 * aux5Term2 - J1 * s31 * J1 * aux2Term1) / J1 +
                                     (aux9Term2 * aux1Term3) / J1)) / J1;
            a3rs[0, 2][2] += (s32 * ((aux3Term2 * aux5Term2 - J1 * s31 * J1 * aux2Term1) / J1 +
                                     (aux9Term2 * aux1Term3) / J1)) / J1;

            a3rs[1, 0][0] += (s30 * ((aux6Term2 * aux1Term2 - J1 * s32 * J1 * aux2Term1) / J1 -
                                     (aux2Term3 * aux7Term2) / J1)) / J1;
            a3rs[1, 0][1] += (s31 * ((aux6Term2 * aux1Term2 - J1 * s32 * J1 * aux2Term1) / J1 -
                                     (aux2Term3 * aux7Term2) / J1)) / J1;
            a3rs[1, 0][2] += (s32 * ((aux6Term2 * aux1Term2 - J1 * s32 * J1 * aux2Term1) / J1 -
                                     (aux2Term3 * aux7Term2) / J1)) / J1;

            a3rs[1, 1][0] += -(s30 * ((aux6Term2 * aux5Term2 + aux3Term2 * aux2Term2) / J1 -
                                      (aux2Term3 * aux8Term2) / J1)) / J1;
            a3rs[1, 1][1] += -(s31 * ((aux6Term2 * aux5Term2 + aux3Term2 * aux2Term2) / J1 -
                                      (aux2Term3 * aux8Term2) / J1)) / J1;
            a3rs[1, 1][2] += -(s32 * ((aux6Term2 * aux5Term2 + aux3Term2 * aux2Term2) / J1 -
                                      (aux2Term3 * aux8Term2) / J1)) / J1;

            a3rs[1, 2][0] += (s30 * ((aux3Term2 * aux1Term2 + J1 * s30 * J1 * aux2Term1) / J1 -
                                     (aux2Term3 * aux9Term2) / J1)) / J1;
            a3rs[1, 2][1] += (s31 * ((aux3Term2 * aux1Term2 + J1 * s30 * J1 * aux2Term1) / J1 -
                                     (aux2Term3 * aux9Term2) / J1)) / J1;
            a3rs[1, 2][2] += (s32 * ((aux3Term2 * aux1Term2 + J1 * s30 * J1 * aux2Term1) / J1 -
                                     (aux2Term3 * aux9Term2) / J1)) / J1;

            a3rs[2, 0][0] += (s30 * ((aux6Term2 * aux2Term2 + J1 * s31 * J1 * aux2Term1) / J1 +
                                     (aux3Term3 * aux7Term2) / J1)) / J1;
            a3rs[2, 0][1] += (s31 * ((aux6Term2 * aux2Term2 + J1 * s31 * J1 * aux2Term1) / J1 +
                                     (aux3Term3 * aux7Term2) / J1)) / J1;
            a3rs[2, 0][2] += (s32 * ((aux6Term2 * aux2Term2 + J1 * s31 * J1 * aux2Term1) / J1 +
                                     (aux3Term3 * aux7Term2) / J1)) / J1;

            a3rs[2, 1][0] += (s30 * ((aux4Term2 * aux2Term2 - J1 * s30 * J1 * aux2Term1) / J1 -
                                     (aux3Term3 * aux8Term2) / J1)) / J1;
            a3rs[2, 1][1] += (s31 * ((aux4Term2 * aux2Term2 - J1 * s30 * J1 * aux2Term1) / J1 -
                                     (aux3Term3 * aux8Term2) / J1)) / J1;
            a3rs[2, 1][2] += (s32 * ((aux4Term2 * aux2Term2 - J1 * s30 * J1 * aux2Term1) / J1 -
                                     (aux3Term3 * aux8Term2) / J1)) / J1;

            a3rs[2, 2][0] += -(s30 * ((aux6Term2 * aux5Term2 + aux4Term2 * aux1Term2) / J1 -
                                      (aux3Term3 * aux9Term2) / J1)) / J1;
            a3rs[2, 2][1] += -(s31 * ((aux6Term2 * aux5Term2 + aux4Term2 * aux1Term2) / J1 -
                                      (aux3Term3 * aux9Term2) / J1)) / J1;
            a3rs[2, 2][2] += -(s32 * ((aux6Term2 * aux5Term2 + aux4Term2 * aux1Term2) / J1 -
                                      (aux3Term3 * aux9Term2) / J1)) / J1;
            #endregion

            #region Term5

            a3rs[0, 0][0] += (2 * s30 * aux1Term3 * aux7Term2) / J1squared;
            a3rs[0, 0][1] += (2 * s31 * aux1Term3 * aux7Term2) / J1squared;
            a3rs[0, 0][2] += (2 * s32 * aux1Term3 * aux7Term2) / J1squared;

            a3rs[0, 1][0] += -(2 * s30 * aux1Term3 * aux8Term2) / J1squared;
            a3rs[0, 1][1] += -(2 * s31 * aux1Term3 * aux8Term2) / J1squared;
            a3rs[0, 1][2] += -(2 * s32 * aux1Term3 * aux8Term2) / J1squared;

            a3rs[0, 2][0] += (2 * s30 * aux9Term2 * aux1Term3) / J1squared;
            a3rs[0, 2][1] += (2 * s31 * aux9Term2 * aux1Term3) / J1squared;
            a3rs[0, 2][2] += (2 * s32 * aux9Term2 * aux1Term3) / J1squared;

            a3rs[1, 0][0] += -(2 * s30 * aux2Term3 * aux7Term2) / J1squared;
            a3rs[1, 0][1] += -(2 * s31 * aux2Term3 * aux7Term2) / J1squared;
            a3rs[1, 0][2] += -(2 * s32 * aux2Term3 * aux7Term2) / J1squared;

            a3rs[1, 1][0] += (2 * s30 * aux2Term3 * aux8Term2) / J1squared;
            a3rs[1, 1][1] += (2 * s31 * aux2Term3 * aux8Term2) / J1squared;
            a3rs[1, 1][2] += (2 * s32 * aux2Term3 * aux8Term2) / J1squared;

            a3rs[1, 2][0] += -(2 * s30 * aux2Term3 * aux9Term2) / J1squared;
            a3rs[1, 2][1] += -(2 * s31 * aux2Term3 * aux9Term2) / J1squared;
            a3rs[1, 2][2] += -(2 * s32 * aux2Term3 * aux9Term2) / J1squared;

            a3rs[2, 0][0] += (2 * s30 * aux3Term3 * aux7Term2) / J1squared;
            a3rs[2, 0][1] += (2 * s31 * aux3Term3 * aux7Term2) / J1squared;
            a3rs[2, 0][2] += (2 * s32 * aux3Term3 * aux7Term2) / J1squared;

            a3rs[2, 1][0] += -(2 * s30 * aux3Term3 * aux8Term2) / J1squared;
            a3rs[2, 1][1] += -(2 * s31 * aux3Term3 * aux8Term2) / J1squared;
            a3rs[2, 1][2] += -(2 * s32 * aux3Term3 * aux8Term2) / J1squared;

            a3rs[2, 2][0] += (2 * s30 * aux3Term3 * aux9Term2) / J1squared;
            a3rs[2, 2][1] += (2 * s31 * aux3Term3 * aux9Term2) / J1squared;
            a3rs[2, 2][2] += (2 * s32 * aux3Term3 * aux9Term2) / J1squared;
            #endregion
            return a3rs;
        }


        private double[,] CalculateA3r(double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double dksi_r, double dheta_r, double J1)
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

            var da3_unit_dr = new double[3, 3];
            da3_unit_dr[0, 0] = -s30 * dnorma3_dr0;
            da3_unit_dr[1, 0] = da3_tilde_dr10 - s31 * dnorma3_dr0;
            da3_unit_dr[2, 0] = da3_tilde_dr20 - s32 * dnorma3_dr0;

            da3_unit_dr[0, 1] = da3_tilde_dr01 - s30 * dnorma3_dr1;
            da3_unit_dr[1, 1] = -s31 * dnorma3_dr1;
            da3_unit_dr[2, 1] = da3_tilde_dr21 - s32 * dnorma3_dr1;

            da3_unit_dr[0, 2] = da3_tilde_dr02 - s30 * dnorma3_dr2;
            da3_unit_dr[1, 2] = da3_tilde_dr12 - s31 * dnorma3_dr2;
            da3_unit_dr[2, 2] = -s32 * dnorma3_dr2;

            return da3_unit_dr;
        }

        internal double[,] CalculateKmembraneNL(ControlPoint[] controlPoints, double[] membraneForces, Nurbs2D nurbs,
            int j)
        {
            var kmembraneNl = new double[controlPoints.Length * 3, controlPoints.Length * 3];

            for (var i = 0; i < controlPoints.Length; i++)
            {
                var dksi_r = nurbs.NurbsDerivativeValuesKsi[i, j];
                var dheta_r = nurbs.NurbsDerivativeValuesHeta[i, j];

                for (int k = 0; k < controlPoints.Length; k++)
                {
                    var dksi_s = nurbs.NurbsDerivativeValuesKsi[k, j];
                    var dheta_s = nurbs.NurbsDerivativeValuesHeta[k, j];

                    var aux = membraneForces[0] * dksi_r * dksi_s + membraneForces[1] * dheta_r * dheta_s +
                              membraneForces[2] * (dksi_r * dheta_s + dksi_s * dheta_r);

                    kmembraneNl[i * 3, k * 3] += aux;
                    kmembraneNl[i * 3 + 1, k * 3 + 1] += aux;
                    kmembraneNl[i * 3 + 2, k * 3 + 2] += aux;
                }
            }

            return kmembraneNl;
        }
        
        internal double[,] CalculateMembraneDeformationMatrix(ControlPoint[] controlPoints, Nurbs2D nurbs, int j,
            double[] surfaceBasisVector1,
            double[] surfaceBasisVector2)
        {
            var dRIa = new double[3, controlPoints.Length * 3];
            for (int i = 0; i < controlPoints.Length; i++)
            {
                for (int m = 0; m < 3; m++)
                {
                    dRIa[m, i] = nurbs.NurbsDerivativeValuesHeta[i, j] * surfaceBasisVector1[m] +
                                 nurbs.NurbsDerivativeValuesKsi[i, j] * surfaceBasisVector2[m];
                }
            }

            var bmembrane = new double[3, controlPoints.Length * 3];
            for (int column = 0; column < controlPoints.Length * 3; column += 3)
            {
                bmembrane[0, column] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[0];
                bmembrane[0, column + 1] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[1];
                bmembrane[0, column + 2] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[2];

                bmembrane[1, column] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[0];
                bmembrane[1, column + 1] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[1];
                bmembrane[1, column + 2] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[2];

                bmembrane[2, column] = dRIa[0, column / 3];
                bmembrane[2, column + 1] = dRIa[1, column / 3];
                bmembrane[2, column + 2] = dRIa[2, column / 3];
            }

            return bmembrane;
        }
        
        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(NurbsKirchhoffLoveShellElementNL shellElement)
        {
            var gauss = new GaussQuadrature();
            var medianSurfaceGP = gauss.CalculateElementGaussPoints(shellElement.Patch.DegreeKsi,
                shellElement.Patch.DegreeHeta, shellElement.Knots.ToList());
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
    }
}