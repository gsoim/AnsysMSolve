using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.IGA.Elements.Boundary;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.MSAnalysis.RveTemplatesPaper;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using JetBrains.dotMemoryUnit;
using MathNet.Numerics.LinearAlgebra;
using Troschuetz.Random;
using Xunit;
using MatlabWriter = MathNet.Numerics.Data.Matlab.MatlabWriter;

[assembly: SuppressXUnitOutputException]

namespace ISAAR.MSolve.IGA.Tests
{
    public class PaperExamples
    {
        [Fact]
        public static void ScordelisLoShell()
        {
            var filename = "ScordelisLoShell";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt")
                .ToString(CultureInfo.InvariantCulture);
            var numberOfRealizations = 1;
            var trandom = new TRandom();

            var youngModulusSolutionPairs = new double[numberOfRealizations, 2];
            for (var realization =0; realization<numberOfRealizations; realization++)
            {
                // Data from https://www.researchgate.net/figure/Mechanical-properties-of-the-nano-hydroxyapatite-polyetheretherketone-nha-PeeK_tbl1_265175039
                var randomInnerE = trandom.Normal(3.4e9, 0.2e9);
                youngModulusSolutionPairs[realization, 0] = randomInnerE;
                var outterMaterial = new ElasticMaterial3DtotalStrain()
                {
                    YoungModulus = 4.3210e8,
                    PoissonRatio = 0.0
                };
                var innerMaterial = new ElasticMaterial3DtotalStrain()
                {
                    YoungModulus = 3.4e9,
                    PoissonRatio = 0.0
                };

                var material3= new ElasticMaterial2D(StressState2D.PlaneStress)
                {
                    YoungModulus = 4.3210e8,
                    PoissonRatio = 0.0
                };
                var homogeneousRveBuilder1 =
                    new CompositeMaterialModeluilderTet2(outterMaterial, innerMaterial, 100, 100, 100);
                //var material = new MicrostructureShell2D(homogeneousRveBuilder1,
                //    microModel => (new SuiteSparseSolver.Builder()).BuildSolver(microModel), false, 1);

                var material4 = new Shell2dRVEMaterialHostConst(1, 1, 1, homogeneousRveBuilder1,
                    constModel => (new SuiteSparseSolver.Builder()).BuildSolver(constModel));

                var modelReader = new IsogeometricShellReader(GeometricalFormulation.NonLinear, filepath, material4);
                var model = modelReader.GenerateModelFromFile();

                model.SurfaceLoads.Add(new SurfaceDistributedLoad(-90, StructuralDof.TranslationY));

                // Rigid diaphragm for AB
                for (var i = 0; i < 19; i++)
                {
                    model.ControlPointsDictionary[i * 19].Constraints
                        .Add(new Constraint() {DOF = StructuralDof.TranslationX});
                    model.ControlPointsDictionary[i * 19].Constraints
                        .Add(new Constraint() {DOF = StructuralDof.TranslationY});
                }

                // Symmetry for CD
                for (var i = 0; i < 19; i++)
                {
                    model.ControlPointsDictionary[i * 19 + 18].Constraints
                        .Add(new Constraint() {DOF = StructuralDof.TranslationZ});

                    model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                        new NodalDof(model.ControlPointsDictionary[i * 19 + 18], StructuralDof.TranslationX),
                        new NodalDof(model.ControlPointsDictionary[i * 19 + 17], StructuralDof.TranslationX)));
                    model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                        new NodalDof(model.ControlPointsDictionary[i * 19 + 18], StructuralDof.TranslationY),
                        new NodalDof(model.ControlPointsDictionary[i * 19 + 17], StructuralDof.TranslationY)));
                }

                // Symmetry for AD
                for (var j = 0; j < 19; j++)
                {
                    model.ControlPointsDictionary[j].Constraints
                        .Add(new Constraint() {DOF = StructuralDof.TranslationX});
                    model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                        new NodalDof(model.ControlPointsDictionary[j], StructuralDof.TranslationY),
                        new NodalDof(model.ControlPointsDictionary[j + 19], StructuralDof.TranslationY)));
                    model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                        new NodalDof(model.ControlPointsDictionary[j], StructuralDof.TranslationZ),
                        new NodalDof(model.ControlPointsDictionary[j + 19], StructuralDof.TranslationZ)));
                }

                // Solvers
                var solverBuilder = new SkylineSolver.Builder();
                ISolver solver = solverBuilder.BuildSolver(model);

                // Structural problem provider
                var provider = new ProblemStructural(model, solver);

                // Linear static analysis
                var childAnalyzer = new LinearAnalyzer(model, solver, provider);
                var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

                // Run the analysis
                parentAnalyzer.Initialize();
                parentAnalyzer.Solve();


                var cp = model.ControlPointsDictionary.Values.Last();
                var dofA = model.GlobalDofOrdering.GlobalFreeDofs[cp, StructuralDof.TranslationY];

                var solution = solver.LinearSystems[0].Solution[dofA];

                youngModulusSolutionPairs[realization, 1] = solution;
            }

            var writer = new Array2DWriter();
            writer.WriteToFile(youngModulusSolutionPairs,
                Path.Combine(Directory.GetCurrentDirectory(), "ScordelisLoMultiscaleResults"));
        }

        [Fact]
        public void SimpleHoodBenchmark()
        {
            var filename = "attempt2";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.iga");

            for (int d = 0; d < 1; d++)
            {
                var numberOfRealizations = 1;
                var trandom = new TRandom();
                var youngModulusSolutionPairs = new double[numberOfRealizations, 2];

                for (int realization = 0; realization < numberOfRealizations; realization++)
                {
                    var randomCnts = trandom.DiscreteUniform(220, 790);
                    youngModulusSolutionPairs[realization, 0] = randomCnts;
                    var outterMaterial = new ElasticMaterial3DtotalStrain()
                    {
                        YoungModulus = 4, //2.79e9,
                        PoissonRatio = 0.4 //0.4
                    };
                    var homogeneousRveBuilder1 = new CntReinforcedElasticNanocomposite(outterMaterial, 200);

                    var material4 = new Shell2dRVEMaterialHostConst(1, 1, 1, homogeneousRveBuilder1,
                        constModel => (new SuiteSparseSolver.Builder()).BuildSolver(constModel));

                    var model = new Model();
                    //var material4= new ShellElasticMaterial2Dtransformationb()
                    //{
                    //    YoungModulus = 100,
                    //    PoissonRatio = 0.3
                    //};
                    var modelReader = new IgaFileReader(model, filepath);
                    modelReader.CreateTSplineShellsModelFromFile(IgaFileReader.TSplineShellType.Thickness, material4);

                    for (int i = 0; i < 100; i++)
                    {
                        var id = model.ControlPoints.ToList()[i].ID;
                        model.ControlPointsDictionary[id].Constraints
                            .Add(new Constraint() {DOF = StructuralDof.TranslationX});
                        model.ControlPointsDictionary[id].Constraints
                            .Add(new Constraint() {DOF = StructuralDof.TranslationY});
                        model.ControlPointsDictionary[id].Constraints
                            .Add(new Constraint() {DOF = StructuralDof.TranslationZ});
                    }

                    for (int i = model.ControlPoints.Count() - 100; i < model.ControlPoints.Count(); i++)
                    {
                        var id = model.ControlPoints.ToList()[i].ID;
                        model.Loads.Add(new Load()
                        {
                            Amount = 100 / 10e4,
                            Node = model.ControlPointsDictionary[id],
                            DOF = StructuralDof.TranslationZ
                        });
                    }

                    var solverBuilder = new SuiteSparseSolver.Builder();
                    ISolver solver = solverBuilder.BuildSolver(model);

                    // Structural problem provider
                    var provider = new ProblemStructural(model, solver);

                    // Linear static analysis
                    var childAnalyzer = new LinearAnalyzer(model, solver, provider);
                    var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

                    // Run the analysis
                    parentAnalyzer.Initialize();
                    var k = solver.LinearSystems[0].Matrix;
                    Matrix<double> kmatlab = MathNet.Numerics.LinearAlgebra.CreateMatrix.Sparse<double>(k.NumRows, k.NumColumns);
                    for (int i = 0; i < k.NumRows; i++)
                    {
                        for (int j = 0; j < k.NumColumns; j++)
                        {
                            kmatlab[i, j] = k[i, j];
                        }
                    }
                    MatlabWriter.Write(Path.Combine(Directory.GetCurrentDirectory(),"KBumper.mat"), kmatlab, "Kff");


                    parentAnalyzer.Solve();


                    //var max=solver.LinearSystems[0].Solution.CopyToArray().Select(Math.Abs).Max();

                    youngModulusSolutionPairs[realization, 1] = solver.LinearSystems[0].Solution[46448];
                }

                var writer = new Array2DWriter();
                writer.WriteToFile(youngModulusSolutionPairs,
                    Path.Combine(Directory.GetCurrentDirectory(), $"BumperMultiscaleResults_{d}"));
            }


            //var paraview= new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution,filename);
            //paraview.CreateParaviewFile();
        }

        [Fact]
        public void TestElasticAndMultiscaleMatricesCnt()
        {
            var numberOfRealizations = 20000;
            var youngModulusSolutionPairs = new double[numberOfRealizations, 2];
            for (int i = 0; i < numberOfRealizations; i++)
            {
                var trandom = new TRandom();
                var randomCnts = trandom.Normal(200, 20);

                var material3 = new ShellElasticMaterial2Dtransformationb()
                {
                    YoungModulus = 4,
                    PoissonRatio = 0.4,
                    TangentVectorV1 = new double[3] {1, 0, 0},
                    TangentVectorV2 = new double[3] {0, 1, 0}
                };

                var outterMaterial = new ElasticMaterial3DtotalStrain()
                {
                    YoungModulus = 4, //2.79e9,
                    PoissonRatio = 0.4 //0.4
                };
                var homogeneousRveBuilder1 = new CntReinforcedElasticNanocomposite(outterMaterial, (int)randomCnts);

                var material4 = new MicrostructureShell2D(homogeneousRveBuilder1,
                    model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1)
                {
                    TangentVectorV1 = new double[3] {1, 0, 0},
                    TangentVectorV2 = new double[3] {0, 1, 0}
                };
                material4.UpdateMaterial(new double[]{0,0,0});


                youngModulusSolutionPairs[i, 0] = (int)randomCnts;
                youngModulusSolutionPairs[i, 0] = material4.ConstitutiveMatrix[0,0];
            }
            var writer = new Array2DWriter();
            writer.WriteToFile(youngModulusSolutionPairs,
                Path.Combine(Directory.GetCurrentDirectory(), "materialCnt200"));
        }


        [Fact]
        public void TestElasticAndMultiscaleMatricesTet()
        {
            var material3 = new ShellElasticMaterial2Dtransformationb()
            {
                YoungModulus = 4.3210,
                PoissonRatio = 0.0,
                TangentVectorV1 = new double[3] {1, 0, 0},
                TangentVectorV2 = new double[3] {0, 1, 0}
            };

            //var trandom = new TRandom();
            //var randomInnerE = trandom.Normal(3.4e9, 0.2e9);
            var outterMaterial = new ElasticMaterial3DtotalStrain()
            {
                YoungModulus = 4.3210,
                PoissonRatio = 0.0
            };
            var innerMaterial = new ElasticMaterial3DtotalStrain()
            {
                YoungModulus = 34,
                PoissonRatio = 0.0
            };

            var homogeneousRveBuilder1 =
                new CompositeMaterialModeluilderTet2(outterMaterial, innerMaterial, 100, 100, 100);

            var material4 = new MicrostructureShell2D(homogeneousRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1)
            {
                TangentVectorV1 = new double[3] {1, 0, 0},
                TangentVectorV2 = new double[3] {0, 1, 0}
            };
            material4.UpdateMaterial(new double[]{0,0,0});
        }
    }
}