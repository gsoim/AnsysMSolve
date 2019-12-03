using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads
{

    public class CompositeMaterialModeluilder : IdegenerateRVEbuilder
    {// test 3d
        double E_outter, ni_outter, E_inner, ni_inner, L01, L02, L03 ;
        double boundarySearchTol;

        //TODO: input material to be cloned.  
        public CompositeMaterialModeluilder(double E_outter, double ni_outter, double E_inner, double ni_inner, double L01, double L02, double L03, double boundarySearchTol = 1e-09)
        {
            this.E_outter = E_outter;
            this.ni_outter = ni_outter;
            this.E_inner = E_inner;
            this.ni_inner = ni_inner;
            this.L01 = L01;
            this.L02 = L02;
            this.L03 = L03;
            this.boundarySearchTol = boundarySearchTol;
        }

        public IRVEbuilder Clone(int a)
        {
            throw new NotImplementedException();
        }

        public Tuple<Model, Dictionary<int, Node>, double> GetModelAndBoundaryNodes()
        {
            Model model = new Model();
            
            model.SubdomainsDictionary[0] = new Subdomain(0);
       
            var (Outter_elements_Node_data, Inner_elements_Node_data, node_coords, NodeIds, boundaryNodesIds) =
                GetModelCreationData();

            for (int i1 = 0; i1 < NodeIds.GetLength(0); i1++)
            {
                int nodeID = NodeIds[i1];
                double nodeCoordX = node_coords[i1, 0];
                double nodeCoordY = node_coords[i1, 1];
                double nodeCoordZ = node_coords[i1, 2];

                model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
            }


            ElasticMaterial3D outerMaterial = new ElasticMaterial3D()
            {
                YoungModulus = E_outter,
                PoissonRatio = ni_outter,
            };

            ElasticMaterial3D innerMaterial = new ElasticMaterial3D()
            {
                YoungModulus = E_inner,
                PoissonRatio = ni_inner,
            };

            //define outer elements group
            for (int i1 = 0; i1 < Outter_elements_Node_data.GetLength(0); i1++)
            {
                Element e1 = new Element()
                {
                    ID = Outter_elements_Node_data[i1,0],
                    ElementType = new Hexa8NonLinear(outerMaterial, GaussLegendre3D.GetQuadratureWithOrder(2,2,2)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };

                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(Outter_elements_Node_data[i1, j+1], model.NodesDictionary[Outter_elements_Node_data[i1, j+1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
            }

            //define outer elements group
            for (int i1 = 0; i1 < Inner_elements_Node_data.GetLength(0); i1++)
            {
                Element e1 = new Element()
                {
                    ID = Inner_elements_Node_data[i1, 0],
                    ElementType = new Hexa8NonLinear(outerMaterial, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };

                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(Inner_elements_Node_data[i1, j + 1], model.NodesDictionary[Inner_elements_Node_data[i1, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
            }


            for (int i2 = 0; i2 < Outter_elements_Node_data.GetLength(0); i2++)
            {
                int subdomainID = 0;
                Element element = model.ElementsDictionary[Outter_elements_Node_data[i2, 0]];
                model.SubdomainsDictionary[subdomainID].Elements.Add(element);
            }

            for (int i2 = 0; i2 < Inner_elements_Node_data.GetLength(0); i2++)
            {
                int subdomainID = 0;
                Element element = model.ElementsDictionary[Inner_elements_Node_data[i2, 0]];
                model.SubdomainsDictionary[subdomainID].Elements.Add(element);
            }

                       
            //for (int i1 = 0; i1 < constraintIds.GetLength(0); i1++)
            //{
            //    model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            //    model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            //    model.NodesDictionary[constraintIds[i1]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });

            //}
            //// Load
            //model.Loads.Add(new Load() { Node = model.NodesDictionary[63], DOF = StructuralDof.TranslationZ, Amount = 1.0 });




            return  new Tuple<Model, Dictionary<int, Node>, double>(model, new Dictionary<int, Node>(), L01*L02*L03); 
        }

        public Dictionary<Node, IList<IDofType>> GetModelRigidBodyNodeConstraints(Model model)
        {
            throw new NotImplementedException();
            // ulopoihsh omoiws me return FEMMeshBuilder.GetConstraintsOfDegenerateRVEForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
        }





        private static (int[,], int[,], double[,], int[], int[]) GetModelCreationData()
        {
            var Outter_elements_Node_data = new int[,] {
            { 1,10,9,11,12,14,13,15,16},
            { 1,10,9,11,12,14,13,15,16} };

            var Inner_elements_Node_data = new int[,] {
            { 1,10,9,11,12,14,13,15,16},
            { 1,10,9,11,12,14,13,15,16} };


            double[,] node_coords = new double[,]
                {{-50.0000000000000000,-50.0000000000000000,50.0000000000000000},
                    {-50.0000000000000000,-50.0000000000000000,-50.0000000000000000},
                    {-50.0000000000000000,50.0000000000000000,50.0000000000000000},
                    {-50.0000000000000000,50.0000000000000000,-50.0000000000000000},
                    {50.0000000000000000,-50.0000000000000000,50.0000000000000000},
                    {50.0000000000000000,-50.0000000000000000,-50.0000000000000000},
                    {50.0000000000000000,50.0000000000000000,50.0000000000000000}};
                    

            int[] NodeIds = new int[7]{1,
                26,
                51,
                76,
                77,
                10,
                17};

            var boundaryNodesIds = new int[3] {1,
                26,
                51};

            return (Outter_elements_Node_data, Inner_elements_Node_data, node_coords, NodeIds, boundaryNodesIds);
        }

             

        
    }
}
