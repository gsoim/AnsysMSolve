﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Filters;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Optimization.Structural.Tests.Topology.SIMP
{
    /// <summary>
    /// Tests for <see cref="TopologySimpLinear2D"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class TopologySimpTests
    {
        private const int subdomainID = 0;

        [Fact]
        private static void TestCantileverBeamGeneral()
        {
            //TODO: Also need a general filter

            // Define the fem analysis and the filter
            double thickness = 1.0;
            var material = new ElasticMaterial2D_v2(StressState2D.PlaneStress) { YoungModulus = 1.0, PoissonRatio = 0.3 };
            var dynamicProperties = new DynamicMaterial(1.0, 0.0, 0.0);

            // Model with 1 subdomain
            var model = new Model_v2();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Generate mesh
            int numElementsX = 32, numElementsY = 20;
            double lengthX = numElementsX;
            double lengthY = numElementsY;
            var mesher = new UniformMeshGenerator2D_v2(0, 0, lengthX, lengthY, numElementsX, numElementsY);
            (IReadOnlyList<Node_v2> nodes, IReadOnlyList<CellConnectivity_v2> connectivity) = mesher.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < nodes.Count; ++n) model.NodesDictionary.Add(n, nodes[n]);

            // Add Quad4 elements to the model
            var factory = new ContinuumElement2DFactory(thickness, material, dynamicProperties);
            for (int e = 0; e < connectivity.Count; ++e)
            {
                ContinuumElement2D element = factory.CreateElement(connectivity[e].CellType, connectivity[e].Vertices);
                var elementWrapper = new Element_v2() { ID = e, ElementType = element };
                foreach (Node_v2 node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(e, elementWrapper);
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }

            // Clamp boundary condition at left edge
            double tol = 1E-10; //TODO: this should be chosen w.r.t. the element size along X
            foreach (var node in model.Nodes.Where(node => Math.Abs(node.X) <= tol))
            {
                node.Constraints.Add(new Constraint() { DOF = DOFType.X, Amount = 0.0 });
                node.Constraints.Add(new Constraint() { DOF = DOFType.Y, Amount = 0.0 });
            }

            // Apply concentrated load at the bottom right corner
            double load = 1.0;
            var cornerNode = model.Nodes.Where(
                node => (Math.Abs(node.X - lengthX) <= tol) && (Math.Abs(node.Y - lengthY) <= tol)); //TODO: this means Y is downwards to match the filer. Change this (and the load)
            Assert.True(cornerNode.Count() == 1);
            model.Loads.Add(new Load_v2() { Amount = load, Node = cornerNode.First(), DOF = DOFType.Y });

            // Define the solver
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);

            // Define the fem analysis and the filter
            double filterAreaRadius = 1.2;
            var fem = new LinearFemAnalysis2DGeneral(model, solver);
            var filter = new MeshIndependentSensitivityFilter2DUniform(numElementsX, numElementsY, filterAreaRadius);

            // Run the test
            TestCantileverBeam(fem, filter);
        }

        [Fact]
        private static void TestCantileverBeamHardCoded()
        {
            // Parameters
            int numElementsX = 32, numElementsY = 20;
            double filterAreaRadius = 1.2;

            // Define the fem analysis and the filter
            var material = new ElasticMaterial2D(StressState2D.PlaneStress) { YoungModulus = 1.0, PoissonRatio = 0.3 };
            var fem = new LinearFemAnalysis2DUniformHardcoded(numElementsX, numElementsY, material,
                LinearFemAnalysis2DUniformHardcoded.BoundaryConditions.ShortCantilever);
            var filter = new MeshIndependentSensitivityFilter2DUniform(numElementsX, numElementsY, filterAreaRadius);

            // Run the test
            TestCantileverBeam(fem, filter);
        }

        private static void TestCantileverBeam(ILinearFemAnalysis fem, IDensityFilter filter)
        {
            // Parameters
            double volumeFraction = 0.4, penalty = 3.0;

            // Define the optimization
            var simp = new TopologySimpLinear2D(fem, filter, volumeFraction);
            simp.PenalizationExponent = penalty;
            var logger = new ObjectiveFunctionLogger();
            simp.Logger = logger;

            // Run the optimization
            simp.Initialize();
            (double compliance, Vector densities) = simp.Optimize();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(TopologySimp99LinesTests.ComplianceHistoryForCantilever);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-11));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_shortcantilever_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Vector densitiesExpected = reader.ReadFile(densitiesPath).Reshape(false);
            Assert.True(densitiesExpected.Equals(densities, 1E-10));
        }

        [Fact]
        private static void TestCantileverBeamWith2LoadCasesHardCoded()
        {
            // Parameters
            int numElementsX = 30, numElementsY = 30;
            double filterAreaRadius = 1.2, volumeFraction = 0.4, penalty = 3.0;

            // Define the optimization
            var material = new ElasticMaterial2D(StressState2D.PlaneStress) { YoungModulus = 1.0, PoissonRatio = 0.3 };
            var fem = new LinearFemAnalysis2DUniformHardcoded(numElementsX, numElementsY, material,
                LinearFemAnalysis2DUniformHardcoded.BoundaryConditions.Cantilever2LoadCases);
            var filter = new MeshIndependentSensitivityFilter2DUniform(numElementsX, numElementsY, filterAreaRadius);
            var simp = new TopologySimpLinear2D(fem, filter, volumeFraction);
            simp.PenalizationExponent = penalty;
            var logger = new ObjectiveFunctionLogger();
            simp.Logger = logger;

            // Run the optimization
            simp.Initialize();
            (double compliance, Vector densities) = simp.Optimize();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(
                TopologySimp99LinesTests.ComplianceHistoryForCantileverWith2LoadCases);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-11));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_cantilever_loadcases_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Vector densitiesExpected = reader.ReadFile(densitiesPath).Reshape(false); 
            Assert.True(densitiesExpected.Equals(densities, 1E-11));
        }

        [Fact]
        private static void TestMbbBeamGeneral()
        {
            // Define the fem analysis and the filter
            double thickness = 1.0;
            var material = new ElasticMaterial2D_v2(StressState2D.PlaneStress) { YoungModulus = 1.0, PoissonRatio = 0.3 };
            var dynamicProperties = new DynamicMaterial(1.0, 0.0, 0.0);

            // Model with 1 subdomain
            var model = new Model_v2();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Generate mesh
            int numElementsX = 60, numElementsY = 20;
            double lengthX = numElementsX;
            double lengthY = numElementsY;
            var mesher = new UniformMeshGenerator2D_v2(0, 0, lengthX, lengthY, numElementsX, numElementsY);
            (IReadOnlyList<Node_v2> nodes, IReadOnlyList<CellConnectivity_v2> connectivity) = mesher.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < nodes.Count; ++n) model.NodesDictionary.Add(n, nodes[n]);

            // Add Quad4 elements to the model
            var factory = new ContinuumElement2DFactory(thickness, material, dynamicProperties);
            for (int e = 0; e < connectivity.Count; ++e)
            {
                ContinuumElement2D element = factory.CreateElement(connectivity[e].CellType, connectivity[e].Vertices);
                var elementWrapper = new Element_v2() { ID = e, ElementType = element };
                foreach (Node_v2 node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(e, elementWrapper);
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }

            // Roller boundary condition at left edge
            double tol = 1E-10; //TODO: this should be chosen w.r.t. the element size along X
            foreach (var node in model.Nodes.Where(node => Math.Abs(node.X) <= tol))
            {
                node.Constraints.Add(new Constraint() { DOF = DOFType.X, Amount = 0.0 });
            }

            // Roller boundary condition at bottom right corner
            var bottomRightNode = model.Nodes.Where(
                node => (Math.Abs(node.X - lengthX) <= tol) && (Math.Abs(node.Y - lengthY) <= tol)); //TODO: this means Y is downwards to match the filer. Change this (and the load)
            Assert.True(bottomRightNode.Count() == 1);
            bottomRightNode.First().Constraints.Add(new Constraint() { DOF = DOFType.Y, Amount = 0.0 });

            // Apply concentrated load at top left corner
            double load = 1.0;
            var topLeftNode = model.Nodes.Where(
                node => (Math.Abs(node.X) <= tol) && (Math.Abs(node.Y) <= tol));
            Assert.True(topLeftNode.Count() == 1);
            model.Loads.Add(new Load_v2() { Amount = load, Node = topLeftNode.First(), DOF = DOFType.Y });

            // Define the solver
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);

            // Define the fem analysis and the filter
            double filterAreaRadius = 1.5;
            var fem = new LinearFemAnalysis2DGeneral(model, solver);
            var filter = new MeshIndependentSensitivityFilter2DUniform(numElementsX, numElementsY, filterAreaRadius);

            // Run the test
            TestMbbBeam(fem, filter);
        }

        [Fact]
        private static void TestMbbBeamHardCoded()
        {
            // Parameters
            int numElementsX = 60, numElementsY = 20;
            double filterAreaRadius = 1.5;

            // Define the optimization
            var material = new ElasticMaterial2D(StressState2D.PlaneStress) { YoungModulus = 1.0, PoissonRatio = 0.3 };
            var fem = new LinearFemAnalysis2DUniformHardcoded(numElementsX, numElementsY, material,
                LinearFemAnalysis2DUniformHardcoded.BoundaryConditions.MbbBeam);
            var filter = new MeshIndependentSensitivityFilter2DUniform(numElementsX, numElementsY, filterAreaRadius);

            TestMbbBeam(fem, filter);
        }

        private static void TestMbbBeam(ILinearFemAnalysis fem, IDensityFilter filter)
        {
            // Parameters
            double volumeFraction = 0.5, penalty = 3.0;

            // Define the optimization
            var simp = new TopologySimpLinear2D(fem, filter, volumeFraction);
            simp.PenalizationExponent = penalty;
            var logger = new ObjectiveFunctionLogger();
            simp.Logger = logger;

            // Run the optimization
            simp.Initialize();
            (double compliance, Vector densities) = simp.Optimize();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(TopologySimp99LinesTests.ComplianceHistoryForMbbBeam);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-10));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_MBBbeam_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Vector densitiesExpected = reader.ReadFile(densitiesPath).Reshape(false);
            Assert.True(densitiesExpected.Equals(densities, 1E-9));
        }
    }
}
