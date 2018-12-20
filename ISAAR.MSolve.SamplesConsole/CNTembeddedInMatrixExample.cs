﻿using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using Xunit;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.FEM.Embedding;
using System.Linq;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using System.IO;

namespace ISAAR.MSolve.SamplesConsole
{
    public class CNTembeddedInMatrixExample
    {
        public static void EmbeddedCNTinMatrix_NewtonRaphson()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Choose model
            EmbeddedModelBuilder.EmbeddedExample(model);
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);

            // Skyline Solver
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int totalDOFs = model.TotalDOFs;
            int increments = 100;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[3601][DOFType.X],
            model.NodalDOFsDictionary[3601][DOFType.Y],
            model.NodalDOFsDictionary[3601][DOFType.Z]});

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            var monitorDOFdisplacement = linearSystems[1].Solution[10709];
            var expectedDOFdisplacement = 97.5308525603;
            Assert.Equal(monitorDOFdisplacement, expectedDOFdisplacement, 3);
        }

        public static void EmbeddedCNTinMatrix_DisplacementControl()
        {

        }

        public static class EmbeddedModelBuilder
        {
            public static void EmbeddedExample(Model model)
            {
                HostElementsBuilder(model);
                EmbeddedElementsBuilder(model);
                //var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key <= 2500000).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 2500000).Select(kv => kv.Value), true);
                var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key <= 2500).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 2500).Select(kv => kv.Value), true);
            }

            public static void HostElementsBuilder(Model model)
            {
                string workingDirectory = @"E:\GEORGE_DATA\DESKTOP\input files"; //"..\..\..\Resources\Beam3DInputFiles";
                string MatrixGeometryFileName = "MATRIX_3D-L_x=5-L_y=5-L_z=100-5x5x100-Geometry_MSolve.inp"; //"MATRIX_3D-L_x=5-L_y=5-L_z=100-50x50x1000-Geometry.inp";
                string MatrixGonnectivityFileName = "MATRIX_3D-L_x=5-L_y=5-L_z=100-5x5x100-ConnMatr_MSolve.inp"; //"MATRIX_3D-L_x=5-L_y=5-L_z=100-50x50x1000-ConnMatr.inp";
                int MatrixNodes = File.ReadLines(workingDirectory + '\\' + MatrixGeometryFileName).Count();
                int MatrixElements = File.ReadLines(workingDirectory + '\\' + MatrixGonnectivityFileName).Count();

                // Nodes Geometry                
                using (TextReader reader = File.OpenText(workingDirectory + '\\' + MatrixGeometryFileName))
                {
                    for (int i = 0; i < MatrixNodes; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int nodeID = int.Parse(bits[0]);
                        double nodeX = double.Parse(bits[1]);
                        double nodeY = double.Parse(bits[2]);
                        double nodeZ = double.Parse(bits[3]);
                        model.NodesDictionary.Add(nodeID, new Node { ID = nodeID, X = nodeX, Y = nodeY, Z = nodeZ });
                    }
                }

                // Boundary Conditions
                for (int iNode = 1; iNode <= 36; iNode++)
                {
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.X });
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Y });
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Z });
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotX });
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotY });
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
                }

                // Create Material
                ElasticMaterial3D solidMaterial = new ElasticMaterial3D()
                {
                    YoungModulus = 3.76,
                    PoissonRatio = 0.3779,
                };

                // Generate elements
                using (TextReader reader = File.OpenText(workingDirectory + '\\' + MatrixGonnectivityFileName))
                {
                    for (int i = 0; i < MatrixElements; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int elementID = int.Parse(bits[0]);
                        int node1 = int.Parse(bits[1]);
                        int node2 = int.Parse(bits[2]);
                        int node3 = int.Parse(bits[3]);
                        int node4 = int.Parse(bits[4]);
                        int node5 = int.Parse(bits[5]);
                        int node6 = int.Parse(bits[6]);
                        int node7 = int.Parse(bits[7]);
                        int node8 = int.Parse(bits[8]);
                        // Hexa8NL element definition
                        Element hexa8NLelement = new Element()
                        {
                            ID = elementID,
                            ElementType = new Hexa8NonLinear(solidMaterial, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
                        };
                        // Add nodes to the created element
                        hexa8NLelement.AddNode(model.NodesDictionary[node1]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node2]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node3]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node4]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node5]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node6]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node7]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node8]);
                        // Add Hexa element to the element and subdomains dictionary of the model
                        model.ElementsDictionary.Add(hexa8NLelement.ID, hexa8NLelement);
                        model.SubdomainsDictionary[1].ElementsDictionary.Add(hexa8NLelement.ID, hexa8NLelement);
                    }                    
                }

                // Add nodal load values at the top nodes of the model
                for (int iNode = 3601; iNode <= 3606; iNode++) //(int iNode = 2603551; iNode < 2603601; iNode++)
                {
                    model.Loads.Add(new Load() { Amount = 1.6666667, Node = model.NodesDictionary[iNode], DOF = DOFType.Y });
                }
                //model.Loads.Add(new Load() { Amount = 5, Node = model.NodesDictionary[3601], DOF = DOFType.Y });
                //model.Loads.Add(new Load() { Amount = 5, Node = model.NodesDictionary[3606], DOF = DOFType.Y });
            }

            public static void EmbeddedElementsBuilder(Model model)
            {
                // define mechanical properties
                double youngModulus = 16710.0;
                double shearModulus = 8080.0;
                double poissonRatio = 0.034; //(youngModulus / (2 * shearModulus)) - 1;
                double area = 5.594673861218848d - 003;  // CNT(20,20)-LinearEBE-TBT-L = 10nm
                double inertiaY = 2.490804749753243D - 006; //1058.55;
                double inertiaZ = 2.490804749753243D - 006; // 1058.55;
                double torsionalInertia = 4.981609499506486D - 006; //496.38;
                double effectiveAreaY = area;
                double effectiveAreaZ = area;
                string workingDirectory = @"E:\GEORGE_DATA\DESKTOP\input files"; //"..\..\..\Resources\Beam3DInputFiles";
                string CNTgeometryFileName = "CNT-8-8-L=100-Geometry-2.inp"; //"EmbeddedCNT-8-8-L=100-h=0-k=1-EBE-L=1-NumberOfCNTs=1-Geometry_beam.inp"; //
                string CNTconnectivityFileName = "CNT-8-8-L=100-ConnMatr-2.inp"; //"EmbeddedCNT-8-8-L=100-h=0-k=1-EBE-L=1-NumberOfCNTs=1-ConnMatr_beam.inp"; //
                int CNTNodes = File.ReadLines(workingDirectory + '\\' + CNTgeometryFileName).Count();
                int CNTElems = File.ReadLines(workingDirectory + '\\' + CNTconnectivityFileName).Count();

                // Geometry
                using (TextReader reader = File.OpenText(workingDirectory + '\\' + CNTgeometryFileName))
                {
                    for (int i = 0; i < CNTNodes; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int nodeID = int.Parse(bits[0]) + 3636; //2603601;
                        double nodeX = double.Parse(bits[1]);
                        double nodeY = double.Parse(bits[2]);
                        double nodeZ = double.Parse(bits[3]);
                        model.NodesDictionary.Add(nodeID, new Node { ID = nodeID, X = nodeX, Y = nodeY, Z = nodeZ });
                    }
                }

                // Create new 3D material
                ElasticMaterial3D beamMaterial = new ElasticMaterial3D
                {
                    YoungModulus = youngModulus,
                    PoissonRatio = poissonRatio,
                };

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                // element nodes

                using (TextReader reader = File.OpenText(workingDirectory + '\\' + CNTconnectivityFileName))
                {
                    for (int i = 0; i < CNTElems; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int elementID = int.Parse(bits[0]) + 2500; //2500000;
                        int node1 = int.Parse(bits[1]) + 3636; //2603601;
                        int node2 = int.Parse(bits[2]) + 3636; //2603601;
                        // element nodes
                        IList<Node> elementNodes = new List<Node>();
                        elementNodes.Add(model.NodesDictionary[node1]);
                        elementNodes.Add(model.NodesDictionary[node2]);
                        // create element
                        var beam_1 = new Beam3DCorotationalQuaternion(elementNodes, beamMaterial, 7.85, beamSection);
                        var beamElement = new Element { ID = elementID, ElementType = beam_1 };
                        // Add nodes to the created element
                        beamElement.AddNode(model.NodesDictionary[node1]);
                        beamElement.AddNode(model.NodesDictionary[node2]);
                        // beam stiffness matrix
                        // var a = beam_1.StiffnessMatrix(beamElement);
                        // Add beam element to the element and subdomains dictionary of the model
                        model.ElementsDictionary.Add(beamElement.ID, beamElement);
                        model.SubdomainsDictionary[1].ElementsDictionary.Add(beamElement.ID, beamElement);
                    }                    
                }
            }            
        }
    }
}
