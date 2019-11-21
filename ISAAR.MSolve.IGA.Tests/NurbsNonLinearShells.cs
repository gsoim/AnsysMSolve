using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
    public class NurbsNonLinearShells
    {
        private const double Tolerance = 1e-9;

        #region Fields
        private List<ControlPoint> ElementControlPoints()
        {
            return new List<ControlPoint>
            {
                new ControlPoint {ID = 0, X = 0.0, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 1, X = 0.0, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 2, X = 0.0, Y =  1, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 3, X = 0.3125, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 4, X = 0.3125, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 5, X = 0.3125, Y =  1, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 6, X = 0.9375, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 7, X = 0.9375, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
                new ControlPoint {ID = 8, X = 0.9375, Y =  1, Z = 0.0, WeightFactor = 1.0},
            };
        }

        private List<Knot> ElementKnots()
        {
            return new List<Knot>()
            {
                new Knot(){ID=0,Ksi=0.0,Heta=0.0,Zeta =0.0 },
                new Knot(){ID=1,Ksi=0.0,Heta=1.0,Zeta =0.0 },
                new Knot(){ID=2,Ksi=0.625,Heta=0.0,Zeta =0.0 },
                new Knot(){ID=3,Ksi=0.625,Heta=1,Zeta =0.0 }
            };
        }

        private Vector KnotValueVectorKsi()
        {
            return Vector.CreateFromArray(new double[]
            {
                0.0000000000000, 0.0000000000000, 0.0000000000000, 0.6250000000000, 1.2500000000000, 1.8750000000000,
                2.5000000000000, 3.1250000000000, 3.7500000000000, 4.3750000000000, 5.0000000000000, 5.6250000000000,
                6.2500000000000, 6.8750000000000, 7.5000000000000, 8.1250000000000, 8.7500000000000, 9.3750000000000,
                10.0000000000000, 10.0000000000000, 10.0000000000000,
            });
        }

        private Vector KnotValueVectorHeta()
        {
            return Vector.CreateFromArray(new double[]
            {
                0.0, 0.0, 0.0, 1.0, 1.0, 1.0
            });
        }

        private ShellElasticMaterial2D Material => new ShellElasticMaterial2D()
        {
            YoungModulus = 1200000,
            PoissonRatio = 0.0,
            TangentVectorV1 = new double[] { 1.000000000000000000000000000000000000, 0.000000000000000053342746886286800000, 0.000000000000000000000000000000000000 },
            TangentVectorV2 = new double[] { 3.90312782094781000000000000000000E-18, 9.99999999999999000000000000000000E-01, 0.00000000000000000000000000000000E+00, },
            NormalVectorV3 = new double[] { 0, 0, 1 }
        };
        
        private NurbsKirchhoffLoveShellElementNL Element
        {
            get
            {
                var patch = new Patch();
                patch.DegreeKsi = 2;
                patch.DegreeHeta = 2;
                patch.NumberOfControlPointsHeta = 3;
                patch.KnotValueVectorKsi = KnotValueVectorKsi();
                patch.KnotValueVectorHeta = KnotValueVectorHeta();
                var element =
                    new NurbsKirchhoffLoveShellElementNL(Material, ElementKnots(), ElementControlPoints(), patch, 0.1);
                element._solution= localSolution;
                return element;
            }
        }

        private double[] localSolution => new double[]
        {
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            0.000000000000000000000000000000000000000, 0.000000000000000000000000000000000000000,
            -0.000001762030370944800000000000000000000, 0.000000000000000018440454162078200000000,
            0.001664779366416430000000000000000000000, -0.000001762030370985190000000000000000000,
            0.000000000000000001055844403292750000000, 0.001664779366416240000000000000000000000,
            -0.000001762030370843590000000000000000000, -0.000000000000000079088297192718400000000,
            0.001664779366416050000000000000000000000
        };

        private double[] MembraneForces => new double[]
        {
            -0.0327209647710278000000000000000000000000000000000,
            -0.0000000000000000000000000000193609885449360000000,
            0.0000000000001690940892323080000000000000000000000,
        };

        private double[] BendingMoments => new double[]
        {
            -0.42618363401210900000000000000000000000000000000000,
            -0.00000000000000000075669834331405100000000000000000,
            0.00000000000000683985998006887000000000000000000000,
        };

            #endregion

        [Fact]
        public void IsogeometricCantileverShell()
        {
            Model model = new Model();
            var filename = "CantileverShellBenchmark16x1";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(),"InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinear);

            Value verticalDistributedLoad = delegate (double x, double y, double z)
            {
                return new double[] { 0, 0, 4 };
            };
            model.Patches[0].EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));

            for (int i = 0; i < 6; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            // Solvers
            var solverBuilder = new DenseMatrixSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 1000);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var logger = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "CantileverBenchmarkLog16x1.txt");
            childAnalyzer.IncrementalLogs.Add(0, logger);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        [Fact]
        public void JacobianTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = shellElement.CalculateJacobian(elementControlPoints, nurbs, 0);

            var expectedJacobian = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "jacobian");
            
            for (var i = 0; i < jacobianMatrix.GetLength(0); i++)
            {
                for (var j = 0; j < jacobianMatrix.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedJacobian[i, j], jacobianMatrix[i, j],
                        Tolerance));
                }
            }
        }

        [Fact]
        public void HessianTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var expectedHessian = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "hessian");

            for (var i = 0; i < hessian.GetLength(0); i++)
            {
                for (var j = 0; j < hessian.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedHessian[i, j], hessian[i, j],
                        Tolerance));
                }
            }
        }

        [Fact]
        public void MembraneDeformationMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = shellElement.CalculateJacobian(elementControlPoints, nurbs, 0);
            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 1);
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

            var Bmembrane = shellElement.CalculateMembraneDeformationMatrix(elementControlPoints, nurbs, 0, surfaceBasisVector1, surfaceBasisVector2);

            var expectedBmembrane = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Bmembrane");

            for (var i = 0; i < Bmembrane.GetLength(0); i++)
            {
                for (var j = 0; j < Bmembrane.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedBmembrane[i, j], Bmembrane[i, j],
                        Tolerance));
                }
            }
        }
        
        [Fact]
        public void BendingDeformationMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = shellElement.CalculateJacobian(elementControlPoints, nurbs, 0);
            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 1);
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

            var surfaceBasisVectorDerivative1 = shellElement.CalculateSurfaceBasisVector1(hessian, 0);
            var surfaceBasisVectorDerivative2 = shellElement.CalculateSurfaceBasisVector1(hessian, 1);
            var surfaceBasisVectorDerivative12 = shellElement.CalculateSurfaceBasisVector1(hessian, 2);

            var Bbending = shellElement.CalculateBendingDeformationMatrix(elementControlPoints, surfaceBasisVector3, nurbs, 0, surfaceBasisVector2,
                surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12);

            var expectedBbending = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Bbending");

            for (var i = 0; i < Bbending.GetLength(0); i++)
            {
                for (var j = 0; j < Bbending.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedBbending[i, j], Bbending[i, j],
                        Tolerance));
                }
            }
        }
        
        
        [Fact]
        public void StiffnessMatrixMembraneNLTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);

            var KmembraneNL = shellElement.CalculateKmembraneNL(elementControlPoints, MembraneForces, nurbs, 0);

            var expectedKmembraneNL = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "KmembraneNL");

            for (var i = 0; i < KmembraneNL.NumRows; i++)
            {
                for (var j = 0; j < KmembraneNL.NumColumns; j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKmembraneNL[i, j], KmembraneNL[i, j],
                        Tolerance));
                }
            }
        }

        [Fact]
        public void StiffnessMatrixBendingNLTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            shellElement._solution = localSolution;
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = shellElement.CalculateJacobian(elementControlPoints, nurbs, 0);
            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 1);
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

            var surfaceBasisVectorDerivative1 = shellElement.CalculateSurfaceBasisVector1(hessian, 0);
            var surfaceBasisVectorDerivative2 = shellElement.CalculateSurfaceBasisVector1(hessian, 1);
            var surfaceBasisVectorDerivative12 = shellElement.CalculateSurfaceBasisVector1(hessian, 2);


            var KbendingNL = shellElement.CalculateKbendingNL(elementControlPoints, BendingMoments, nurbs,
                Vector.CreateFromArray(surfaceBasisVector1), 
                Vector.CreateFromArray(surfaceBasisVector2),
                Vector.CreateFromArray(surfaceBasisVector3),
                Vector.CreateFromArray(surfaceBasisVectorDerivative1),
                Vector.CreateFromArray(surfaceBasisVectorDerivative2),
                Vector.CreateFromArray(surfaceBasisVectorDerivative12), J1, 0);

            var expectedKbendingNL = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "KbendingNL");


            for (var i = 0; i < KbendingNL.NumRows; i++)
            {
                for (var j = 0; j < KbendingNL.NumColumns; j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKbendingNL[i, j], KbendingNL[i, j],
                        Tolerance));
                }
            }
        }


        [Fact]
        public void ConstitutiveMatrixThicknessIntegration()
        {
            throw new NotImplementedException();
        }

        [Fact]
        public void StressesThicknessIntegration()
        {
            throw new NotImplementedException();
        }

        [Fact]
        public void TotalElementStiffnessMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);
            var Ktotal=shellElement.StiffnessMatrix(shellElement);

            var expectedKtotal = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Ktotal");
            
            for (var i = 0; i < Ktotal.NumRows; i++)
            {
                for (var j = 0; j < Ktotal.NumColumns; j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKtotal[i, j], Ktotal[i, j],
                        Tolerance));
                }
            }
        }
    }
}
