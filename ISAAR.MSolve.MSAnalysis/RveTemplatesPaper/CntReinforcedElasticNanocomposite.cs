using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using Troschuetz.Random;

namespace ISAAR.MSolve.MSAnalysis.RveTemplatesPaper
{
    public class CntReinforcedElasticNanocomposite //: IdegenerateRVEbuilder
    {
        int hexa1 = 5;
        int hexa2 = 5;
        int hexa3 = 5;

        double L01 = 100;
        double L02 = 100;
        double L03 = 100;
        IIsotropicContinuumMaterial3D matrixMaterial;


        // cnt paramaters
        IIsotropicContinuumMaterial3D CntMaterial;
        int numberOfCnts;
        // define mechanical properties
        double youngModulus = 1.0; // 5490; // 
        double shearModulus = 1.0; // 871; // 
        double poissonRatio;  //2.15; // 0.034;
        double area = 1776.65;  // CNT(20,20)-LinearEBE-TBT-L = 10nm
        double inertiaY = 1058.55;
        double inertiaZ = 1058.55;
        double torsionalInertia = 496.38;
        double effectiveAreaY; 
        double effectiveAreaZ;
        int subdomainID = 0;

        public CntReinforcedElasticNanocomposite(IIsotropicContinuumMaterial3D matrixMaterial, int numberOfCnts)
        {
            this.matrixMaterial = matrixMaterial;

            this.CntMaterial = new ElasticMaterial3D()
            {
                YoungModulus = youngModulus,
                PoissonRatio = (youngModulus / (2 * shearModulus)) - 1
            };

            double effectiveAreaY = area;
            double effectiveAreaZ = area;

            this.numberOfCnts = numberOfCnts;
        }

        public Tuple<Model, Dictionary<int, Node>, double> GetModelAndBoundaryNodes()
        {
            (int[] NodeIds, double[,] node_coords) = GetHexaRveNodesData();
            int [,] elementConnectivity= GetHexaRveConnectivity();

            Model model= new Model();
            AddHexaElements(model,NodeIds,node_coords, elementConnectivity);
            var hostElements = model.Elements.Count();
            var hostNodes = model.Nodes.Count();


            (int[] cntNodeIds, double[,] cntNodecoords, int[,] cntElementConnectivity) = GetCntBeamsNodesData(hostNodes, hostElements);
            AddCntBeamElements(model, cntNodeIds, cntNodecoords, cntElementConnectivity);
            
            var embeddedGrouping = EmbeddedBeam3DGrouping.CreateFullyBonded(model, model.ElementsDictionary
                        .Where(x => x.Key <= hostElements).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > hostElements)
                        .Select(kv => kv.Value), true);

            throw new NotImplementedException();

        }

        private void AddCntBeamElements(Model model, int[] cntNodeIds, double[,] cntNodeCoordinates, int[,] cntElementConnectivity)
        {
            var hostElements = model.Elements.Count;
            for (int i = 0; i < cntNodeIds.Length; i++)
            {
                var nodeID = cntNodeIds[i];
                var coordX = cntNodeCoordinates[i, 0];
                var coordY = cntNodeCoordinates[i, 1];
                var coordZ = cntNodeCoordinates[i, 2];

                model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: coordX, y: coordY, z: coordZ));
            }
            
            var beamSection =
                new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

            for (int i = 0; i < cntElementConnectivity.GetLength(0); i++)
            {
                var elementNodes = new List<Node>
                {
                    model.NodesDictionary[cntElementConnectivity[i, 0]],
                    model.NodesDictionary[cntElementConnectivity[i, 1]]
                };

                var beam_1 = new Beam3DCorotationalQuaternion(elementNodes, CntMaterial, 7.85, beamSection);
                var beamElement = new Element { ID = i+hostElements, ElementType = beam_1 };

                beamElement.AddNode(elementNodes[0]);
                beamElement.AddNode(elementNodes[1]);

                model.ElementsDictionary.Add(beamElement.ID, beamElement);
                model.SubdomainsDictionary[0].Elements.Add(beamElement);
            }
            
        }

        private (int[] cntNodeIds, double[,] cntNodeCoordinates, int[,] cntElementConnectivity) GetCntBeamsNodesData(int hostNodes, int hostElements)
        {
            var cntNodeIds = new int[numberOfCnts*2];
            var cntNodeCoordinates = new double[numberOfCnts * 2, 3];
            var cntElementConnectivity = new int[numberOfCnts, 2];

            var trandom = new TRandom();
            for (int i = 0; i < numberOfCnts; i++)
            {
                var randomY = trandom.ContinuousUniform(0, L02);
                var randomZ = trandom.ContinuousUniform(0, L03);

                cntNodeIds[2 * i] = hostNodes + 2 * i;
                cntNodeIds[2 * i + 1] = hostNodes + 2 * i + 1;

                cntNodeCoordinates[2 * i, 0] = 1;
                cntNodeCoordinates[2 * i, 1] = randomY;
                cntNodeCoordinates[2 * i, 2] = randomZ;

                cntNodeCoordinates[2 * i + 1, 0] = L01-1;
                cntNodeCoordinates[2 * i + 1, 1] = randomY;
                cntNodeCoordinates[2 * i + 1, 2] = randomZ;

                cntElementConnectivity[i, 0] = cntNodeIds[2 * i];
                cntElementConnectivity[i, 1] = cntNodeIds[2 * i + 1];
            }


            return (cntNodeIds, cntNodeCoordinates, cntElementConnectivity);
        }

        private void AddHexaElements(Model model, int[] nodeIds, double[,] node_coords, int[,] elementConnectivity)
        {
            for (int i1 = 0; i1 < nodeIds.Length; i1++)
            {
                int nodeID = nodeIds[i1];
                double nodeCoordX = node_coords[i1, 0];
                double nodeCoordY = node_coords[i1, 1];
                double nodeCoordZ = node_coords[i1, 2];

                model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
            }
            
            var factoryOutter = new ContinuumElement3DFactory(matrixMaterial, null);
            for (int i1 = 0; i1 < elementConnectivity.GetLength(0); i1++)
            {
                List<Node> nodeSet = new List<Node>();
                for (int j = 0; j < 8; j++)
                {
                    int nodeID = elementConnectivity[i1, j];
                    nodeSet.Add((Node)model.NodesDictionary[nodeID]);
                }

                Element e1 = new Element()
                {
                    ID = elementConnectivity[i1, 0],
                    ElementType = factoryOutter.CreateElement(CellType.Hexa8, nodeSet) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };

                for (int j = 0; j < 8; j++)
                {
                    int nodeID = elementConnectivity[i1, j];
                    e1.NodesDictionary.Add(nodeID, model.NodesDictionary[nodeID]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
            }
        }

        private int[,] GetHexaRveConnectivity()
        {
            int numberOfElements = hexa1 * hexa2 * hexa3;
            int[,] elementConnectivity = new int[numberOfElements,8];

            int counterElement = 0;
            for (int i = 0; i < hexa3; i++)
            {
                for (int j = 0; j < hexa2; j++)
                {
                    for (int k = 0; k < hexa1; k++)
                    {
                        elementConnectivity[counterElement, 0] = k + (hexa1 + 1) * j + (hexa1 + 1) * (hexa2 + 1) * i;
                        elementConnectivity[counterElement, 1] = elementConnectivity[counterElement, 0] + 1;
                        elementConnectivity[counterElement, 2] = elementConnectivity[counterElement, 0] + (hexa1 + 1);
                        elementConnectivity[counterElement, 3] = elementConnectivity[counterElement, 2] + 1;

                        elementConnectivity[counterElement, 4] = k + (hexa1 + 1) * j + (hexa1 + 1) * (hexa2 + 1) * (i+1);
                        elementConnectivity[counterElement, 5] = elementConnectivity[counterElement, 4] + 1;
                        elementConnectivity[counterElement, 6] = elementConnectivity[counterElement, 4] + (hexa1 + 1);
                        elementConnectivity[counterElement, 7] = elementConnectivity[counterElement, 6] + 1;
                        
                        counterElement++;
                    }
                }
            }

            return elementConnectivity;
        }

        private (int[] NodeIds, double[,] node_coords) GetHexaRveNodesData()
        {
            int numberOfNodes = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            int[] NodeIds = new int[numberOfNodes];
            double[,] node_coords = new double[numberOfNodes, 3];

            int counterNode=0;

            for (int  i = 0;  i <= hexa3;  i++)
            {
                for (int j = 0; j <= hexa2; j++)
                {
                    for (int k = 0; k <= hexa1; k++)
                    {
                        NodeIds[counterNode] = counterNode;
                        node_coords[counterNode, 0] = k * (L01 / hexa1) - L01 * 0.5;
                        node_coords[counterNode, 1] = j * (L02 / hexa2) - L02 * 0.5;
                        node_coords[counterNode, 2] = i * (L03 / hexa3) - L03 * 0.5;
                        counterNode++;
                    }
                }
            }

            return (NodeIds, node_coords);
        }




    }
}
