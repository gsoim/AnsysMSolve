using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
    public class NurbsNonLinearShells
    {
        [Fact]
        public void IsogeometricCantileverShell()
        {
            Model model = new Model();
            var filename = "CantileverShellBenchmark16x1";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), $"{filename}.txt");
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
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 10000);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var logger = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 10000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "CantileverBenchmarkLog16x1.txt");
            childAnalyzer.IncrementalLogs.Add(0, logger);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}
