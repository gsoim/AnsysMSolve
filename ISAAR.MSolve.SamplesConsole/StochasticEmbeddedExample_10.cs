using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Problems;
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
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.FEM.Postprocessing;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using static ISAAR.MSolve.SamplesConsole.StochasticEmbeddedExample_10.Run2a_Elastic;

namespace ISAAR.MSolve.SamplesConsole
{
    public class StochasticEmbeddedExample_10
    {
        public static class Run2a_Elastic
        {
            private const int subdomainID = 0;
            private const int hostElements = 1;
            private const int hostNodes = 8;
            private const int embeddedElements = 1;
            private const int embeddedNodes = 2;
            private const double nodalLoad = +10.0; // +1000.0;//
            private const double nodalDisplacement = +10.0;
            private const int increments = 10;
            private const int monitorNode = 3;

            public static void SingleMatrix_NewtonRaphson_Stochastic(int noStochasticSimulation)
            {
                // Model creation
                var model = new Model();

                // Subdomains
                //model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = 1 });
                model.SubdomainsDictionary.Add(0, new Subdomain(0));

                // Choose model
                EbeEmbeddedModelBuilder.SingleMatrixBuilder_Stochastic(model, noStochasticSimulation);

                // Boundary Conditions - [Left-End]
                for (int iNode = 5; iNode <= 8; iNode++)
                {
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX});
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY});
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
                }                               

                // Loading Conditions - [Right-End] - {2 nodes}
                for (int iNode = 3; iNode <= 4; iNode++)
                {
                    model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[iNode], DOF = StructuralDof.TranslationY });
                }

                // Choose linear equation system solver
                //var solverBuilder = new SkylineSolver.Builder();
                //SkylineSolver solver = solverBuilder.BuildSolver(model);
                var solverBuilder = new SuiteSparseSolver.Builder();
                SuiteSparseSolver solver = solverBuilder.BuildSolver(model);

                // Choose the provider of the problem -> here a structural problem
                var provider = new ProblemStructural(model, solver);

                // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer     
                var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments)
                {
                    MaxIterationsPerIncrement = 100,
                    NumIterationsForMatrixRebuild = 1,
                    ResidualTolerance = 5E-03
                };

                LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();

                // Choose parent analyzer -> Parent: Static
                var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

                // Request output
                

                // Run the analysis
                parentAnalyzer.Initialize();
                parentAnalyzer.Solve();
            }

            public static void EBEembeddedInMatrix_NewtonRaphson_Stochastic(int noStochasticSimulation)
            {
                //VectorExtensions.AssignTotalAffinityCount();

                // Model creation
                var model = new Model();

                // Subdomains
                //model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = 1 });
                model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

                // Choose model
                EbeEmbeddedModelBuilder.FullyBondedEmbeddedBuilder_Stochastic(model, noStochasticSimulation);

                // Boundary Conditions - [Left - End]
                for (int iNode = 5; iNode <= 8; iNode++)
                {
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX});
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
                }

                // Loading Conditions - [Right-End] - {2 nodes}
                for (int iNode = 3; iNode <= 4; iNode++)
                {
                    model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[iNode], DOF = StructuralDof.TranslationY });
                }

                // Choose linear equation system solver
                var solverBuilder = new SkylineSolver.Builder();
                SkylineSolver solver = solverBuilder.BuildSolver(model);
                //var solverBuilder = new SuiteSparseSolver.Builder();
                //SuiteSparseSolver solver = solverBuilder.BuildSolver(model);

                // Choose the provider of the problem -> here a structural problem
                var provider = new ProblemStructural(model, solver);

                // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer     
                var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments)
                {
                    MaxIterationsPerIncrement = 100,
                    NumIterationsForMatrixRebuild = 1,
                    ResidualTolerance = 5E-03
                };

                LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();

                // Choose parent analyzer -> Parent: Static
                var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

                // Request output
                //string currentOutputFileName = "Run2a-Stochastic-CNT-Results.txt";
                //string extension = Path.GetExtension(currentOutputFileName);
                //string pathName = outputDirectory;
                //string fileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(currentOutputFileName));
                //string outputFile = string.Format("{0}_{1}{2}", fileNameOnly, noStochasticSimulation, extension);
                //var logger = new TotalLoadsDisplacementsPerIncrementLog(model.SubdomainsDictionary[subdomainID], increments,
                //    model.NodesDictionary[monitorNode], monitorDof, outputFile);
                //childAnalyzer.IncrementalLogs.Add(subdomainID, logger);

                // Run the analysis
                parentAnalyzer.Initialize();
                parentAnalyzer.Solve();
            }

            public static class EbeEmbeddedModelBuilder
            {


                private static int hostElements;
                private static int hostNodes;

                public static void SingleMatrixBuilder_Stochastic(Model model, int i)
                {
                    HostElements(model);
                }

                public static void FullyBondedEmbeddedBuilder_Stochastic(Model model, int i)
                {
                    HostElements(model);
                    EmbeddedElements_Stochastic(model, i);
                    var embeddedGrouping = EmbeddedBeam3DGrouping.CreateFullyBonded(model, model.ElementsDictionary
                        .Where(x => x.Key <= hostElements).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > hostElements)
                        .Select(kv => kv.Value), true);
                }

                private static void HostElements(Model model)
                {
                    string workingDirectory = Path.Combine(Directory.GetCurrentDirectory(), "InputFilesEmbBeam");
                    string MatrixGeometryFileName = "MATRIX_3D-L_x=10-L_y=10-L_z=20-1x1x1-Geometry_MSolve.inp";
                    string MatrixConnectivityFileName = "MATRIX_3D-L_x=10-L_y=10-L_z=20-1x1x1-ConnMatr_MSolve.inp";
                    int matrixNodes = File.ReadLines(workingDirectory + '\\' + MatrixGeometryFileName).Count();
                    int matrixElements = File.ReadLines(workingDirectory + '\\' + MatrixConnectivityFileName).Count();
                    hostElements = matrixElements;
                    hostNodes = matrixNodes;
                    // Nodes Geometry
                    using (TextReader reader = File.OpenText(workingDirectory + '\\' + MatrixGeometryFileName))
                    {
                        for (int i = 0; i < matrixNodes; i++)
                        {
                            string text = reader.ReadLine();
                            string[] bits = text.Split(',');
                            int nodeID = int.Parse(bits[0]);
                            double nodeX = double.Parse(bits[1]);
                            double nodeY = double.Parse(bits[2]);
                            double nodeZ = double.Parse(bits[3]);
                            model.NodesDictionary.Add(nodeID, new Node(nodeID, nodeX, nodeY, nodeZ));
                        }
                    }

                    // Create Elastic Material
                    var solidMaterial = new ElasticMaterial3D()
                    {
                        YoungModulus = 4.00,
                        PoissonRatio = 0.40,
                    };

                    // Generate elements
                    using (TextReader reader = File.OpenText(workingDirectory + '\\' + MatrixConnectivityFileName))
                    {
                        for (int i = 0; i < matrixElements; i++)
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

                            var hexa8NLelement = new Element()
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
                            //model.SubdomainsDictionary[0].ElementsDictionary.Add(hexa8NLelement.ID, hexa8NLelement);
                            model.SubdomainsDictionary[0].Elements.Add(hexa8NLelement);
                        }
                    }
                }

                private static void EmbeddedElements_Stochastic(Model model, int noStochasticSimulation)
                {
                    // define mechanical properties
                    double youngModulus = 1.0; // 5490; // 
                    double shearModulus = 1.0; // 871; // 
                    double poissonRatio = (youngModulus / (2 * shearModulus)) - 1; //2.15; // 0.034;
                    double area = 1776.65;  // CNT(20,20)-LinearEBE-TBT-L = 10nm
                    double inertiaY = 1058.55;
                    double inertiaZ = 1058.55;
                    double torsionalInertia = 496.38;
                    double effectiveAreaY = area;
                    double effectiveAreaZ = area;
                    string workingDirectory = Path.Combine(Directory.GetCurrentDirectory(), "InputFilesEmbBeam");

                    string CNTgeometryFileName = "nodes.txt";
                    string CNTconnectivityFileName = "connectivity.txt";

                    string fileNameOnlyCNTgeometryFileName = Path.Combine(workingDirectory, Path.GetFileNameWithoutExtension(CNTgeometryFileName));
                    string fileNameOnlyCNTconnectivityFileName = Path.Combine(workingDirectory, Path.GetFileNameWithoutExtension(CNTconnectivityFileName));
                    string extension = Path.GetExtension(CNTgeometryFileName);
                    string extension_2 = Path.GetExtension(CNTconnectivityFileName);

                    string currentCNTgeometryFileName = string.Format("{0}_{1}{2}", fileNameOnlyCNTgeometryFileName, noStochasticSimulation, extension);

                    //string currentCNTconnectivityFileName = string.Format("{0}_{1}{2}", fileNameOnlyCNTconnectivityFileName, noStochasticSimulation, extension_2);
                    string currentCNTconnectivityFileName = string.Format("{0}{1}", fileNameOnlyCNTconnectivityFileName, extension_2);

                    int CNTNodes = File.ReadLines(currentCNTgeometryFileName).Count();
                    int CNTElems = File.ReadLines(currentCNTconnectivityFileName).Count();

                    // Geometry
                    using (TextReader reader = File.OpenText(currentCNTgeometryFileName))
                    {
                        for (int i = 0; i < CNTNodes; i++)
                        {
                            string text = reader.ReadLine();
                            string[] bits = text.Split(',');
                            int nodeID = int.Parse(bits[0]) + hostNodes; // matrixNodes
                            double nodeX = double.Parse(bits[1]);
                            double nodeY = double.Parse(bits[2]);
                            double nodeZ = double.Parse(bits[3]);
                            model.NodesDictionary.Add(nodeID, new Node(nodeID, nodeX, nodeY, nodeZ));
                        }
                    }

                    // Create new 3D material
                    var beamMaterial = new ElasticMaterial3D
                    {
                        YoungModulus = youngModulus,
                        PoissonRatio = poissonRatio,
                    };

                    // Create new Beam3D section and element
                    var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                    // element nodes

                    using (TextReader reader = File.OpenText(currentCNTconnectivityFileName))
                    {
                        for (int i = 0; i < CNTElems; i++)
                        {
                            string text = reader.ReadLine();
                            string[] bits = text.Split(',');
                            int elementID = int.Parse(bits[0]) + hostElements; // matrixElements
                            int node1 = int.Parse(bits[1]) + hostNodes; // matrixNodes
                            int node2 = int.Parse(bits[2]) + hostNodes; // matrixNodes
                                                                        // element nodes
                            var elementNodes = new List<Node>();
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
                            //model.SubdomainsDictionary[0].ElementsDictionary.Add(beamElement.ID, beamElement);
                            model.SubdomainsDictionary[0].Elements.Add(beamElement);
                        }
                    }
                }


            }
        }

    }
}
