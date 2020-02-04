using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.FEM.Postprocessing
{
    public class ParaviewEmbedded3D
    {
        private Model _model;
        private IVectorView _solution;
        private string _filename;

        public ParaviewEmbedded3D(Model model, IVectorView solution, string filename)
        {
            _model = model;
            _solution = solution;
            _filename = filename;
        }

        public void CreateParaviewFile()
        {
            WriteParaviewFile3D();
        }

        
        private void WriteParaviewFile3D()
        {
            WriteParaviewHostElements();
            WriteParaviewEmbeddedElements();
            //WriteParaviewEmbeddedCohesiveElements();
        }

        private void WriteParaviewHostElements()
        {
            int[] conn = new int[] { 6, 7, 4, 5, 2, 3, 0, 1 };
            var elements = _model.Elements.Where(e => e.ElementType is Hexa8NonLinear).ToList();
            var nodes = new List<INode>();
            elements.ForEach(e => nodes.AddRange(e.Nodes));
            nodes = nodes.Distinct().ToList();

            var numberOfPoints = nodes.Count;
            var numberOfCells = elements.Count;

            int numberOfVerticesPerCell = 0;
            int paraviewCellCode = 0;

            numberOfVerticesPerCell = 8;

            paraviewCellCode = 12;

            using (StreamWriter outputFile = new StreamWriter(Path.Combine(Directory.GetCurrentDirectory(), $"{_filename}_Paraview_Host.vtu")))
            {
                outputFile.WriteLine("<?xml version=\"1.0\"?>");
                outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">");
                outputFile.WriteLine("<UnstructuredGrid>");

                outputFile.WriteLine($"<Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
                outputFile.WriteLine($"<PointData Vectors=\"U\">");
                outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"U\" format=\"ascii\" NumberOfComponents=\"3\">");

                for (int i = 0; i < numberOfPoints; i++)
                {
                    var node = nodes[i];
                    var dx = 0;
                    var dy = 0;
                    var dz = 0;

                    outputFile.WriteLine($"{dx} {dy} {dz}");
                }

                outputFile.WriteLine("</DataArray>");

                outputFile.WriteLine("</PointData>");
                outputFile.WriteLine("<Points>");
                outputFile.WriteLine("<DataArray type=\"Float32\" NumberOfComponents=\"3\">");

                for (int i = 0; i < numberOfPoints; i++)
                    outputFile.WriteLine($"{nodes[i].X} {nodes[i].Y} {nodes[i].Z}");

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("</Points>");
                outputFile.WriteLine("<Cells>");
                outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"connectivity\">");

                for (int i = 0; i < numberOfCells; i++)
                {
                    for (int j = 0; j < elements[i].Nodes.Count; j++)
                    {
                        var elementNode = elements[i].Nodes[conn[j]];
                        var indexOfElementNode = nodes.IndexOf(elementNode);
                        outputFile.Write($"{indexOfElementNode} ");
                    }
                    outputFile.WriteLine("");
                }

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"offsets\">");

                var offset = 0;
                for (int i = 0; i < numberOfCells; i++)
                {
                    offset += numberOfVerticesPerCell;
                    outputFile.WriteLine(offset);
                }

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"types\">");

                for (int i = 0; i < numberOfCells; i++)
                    outputFile.WriteLine(paraviewCellCode);

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("</Cells>");
                outputFile.WriteLine("</Piece>");
                outputFile.WriteLine("</UnstructuredGrid>");
                outputFile.WriteLine("</VTKFile>");
            }
        }

        private void WriteParaviewEmbeddedElements()
        {
            int[] conn = new int[] { 0, 1 };
            var elements = _model.Elements.Where(e => e.ElementType is Beam3DCorotationalQuaternion).ToList();
            var nodes = new List<INode>();
            elements.ForEach(e => nodes.AddRange(e.Nodes));
            nodes = nodes.Distinct().ToList();

            var numberOfPoints = nodes.Count;
            var numberOfCells = elements.Count;

            int numberOfVerticesPerCell = 0;
            int paraviewCellCode = 0;

            numberOfVerticesPerCell = 2;

            paraviewCellCode = 3;

            Dictionary<INode, double[]> nodalDisplacements = new Dictionary<INode, double[]>();
            foreach (var element in elements)
            {
                //var embeddedBeamDisplacements = _model.SubdomainsDictionary[0].CalculateElementDisplacements(element, _solution);
                //var beamDisplacements = element.ElementType.DofEnumerator.GetTransformedDisplacementsVector(embeddedBeamDisplacements);
                if (!nodalDisplacements.ContainsKey(element.Nodes[0]))
                {
                    nodalDisplacements.Add(element.Nodes[0], new double[] { 0, 0, 0 });
                }
                if (!nodalDisplacements.ContainsKey(element.Nodes[1]))
                {
                    nodalDisplacements.Add(element.Nodes[1], new double[] { 0, 0, 0 });
                }
            }

            using (StreamWriter outputFile = new StreamWriter(Path.Combine(Directory.GetCurrentDirectory(), $"{_filename}_Paraview_Embedded.vtu")))
            {
                outputFile.WriteLine("<?xml version=\"1.0\"?>");
                outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">");
                outputFile.WriteLine("<UnstructuredGrid>");

                outputFile.WriteLine($"<Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
                outputFile.WriteLine($"<PointData Vectors=\"U\">");
                outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"U\" format=\"ascii\" NumberOfComponents=\"3\">");

                for (int i = 0; i < numberOfPoints; i++)
                {
                    var node = nodes[i];
                    var dx = nodalDisplacements[node][0];
                    var dy = nodalDisplacements[node][1];
                    var dz = nodalDisplacements[node][2];

                    outputFile.WriteLine($"{dx} {dy} {dz}");
                }

                outputFile.WriteLine("</DataArray>");

                outputFile.WriteLine("</PointData>");
                outputFile.WriteLine("<Points>");
                outputFile.WriteLine("<DataArray type=\"Float32\" NumberOfComponents=\"3\">");

                for (int i = 0; i < numberOfPoints; i++)
                    outputFile.WriteLine($"{nodes[i].X} {nodes[i].Y} {nodes[i].Z}");

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("</Points>");
                outputFile.WriteLine("<Cells>");
                outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"connectivity\">");

                for (int i = 0; i < numberOfCells; i++)
                {
                    for (int j = 0; j < elements[i].Nodes.Count; j++)
                    {
                        var elementNode = elements[i].Nodes[conn[j]];
                        var indexOfElementNode = nodes.IndexOf(elementNode);
                        outputFile.Write($"{indexOfElementNode} ");
                    }
                    outputFile.WriteLine("");
                }

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"offsets\">");

                var offset = 0;
                for (int i = 0; i < numberOfCells; i++)
                {
                    offset += numberOfVerticesPerCell;
                    outputFile.WriteLine(offset);
                }

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"types\">");

                for (int i = 0; i < numberOfCells; i++)
                    outputFile.WriteLine(paraviewCellCode);

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("</Cells>");
                outputFile.WriteLine("</Piece>");
                outputFile.WriteLine("</UnstructuredGrid>");
                outputFile.WriteLine("</VTKFile>");
            }
        }

        //private void WriteParaviewEmbeddedCohesiveElements()
        //{
        //    int[] conn = new int[] { 2, 3 };  // new int[] { 0, 1 };
        //    var elements = _model.Elements.Where(e => e.ElementType is CohesiveBeam3DToBeam3D).ToList();
            
        //    var hostElements = _model.Elements.Where(e => e.ElementType is Hexa8NonLinear_v2).ToList();
        //    var hostNodes = new List<INode>();
        //    hostElements.ForEach(e => hostNodes.AddRange(e.Nodes));
        //    var numberOfHostNodes = hostNodes.Distinct().ToList().Count;

        //    var beamElements = _model.Elements.Where(e => e.ElementType is Beam3DCorotationalQuaternion_v2).ToList();
        //    var beamNodes = new List<INode>();
        //    beamElements.ForEach(e => beamNodes.AddRange(e.Nodes));
        //    var numberOfBeamNodes = beamNodes.Distinct().ToList().Count;

        //    var nodes = new List<INode>();
        //    elements.ForEach(e => nodes.AddRange(e.Nodes));
        //    nodes = nodes.Where(n => n.ID > (numberOfHostNodes + numberOfBeamNodes)).Distinct().ToList();

        //    var numberOfPoints = nodes.Count;
        //    var numberOfCells = elements.Count;

        //    int numberOfVerticesPerCell = 0;
        //    int paraviewCellCode = 0;

        //    numberOfVerticesPerCell = 2;

        //    paraviewCellCode = 3;

        //    Dictionary<INode, double[]> nodalDisplacements = new Dictionary<INode, double[]>();
        //    foreach (var element in elements)
        //    {
        //        var embeddedBeamDisplacements = _model.SubdomainsDictionary[0].CalculateElementDisplacements(element, _solution);
        //        var beamDisplacements = element.ElementType.DofEnumerator.GetTransformedDisplacementsVector(embeddedBeamDisplacements);
        //        if (!nodalDisplacements.ContainsKey(element.Nodes[2]))
        //        {
        //            nodalDisplacements.Add(element.Nodes[2], new double[] { beamDisplacements[0], beamDisplacements[1], beamDisplacements[2] });
        //        }
        //        if (!nodalDisplacements.ContainsKey(element.Nodes[3]))
        //        {
        //            nodalDisplacements.Add(element.Nodes[3], new double[] { beamDisplacements[6], beamDisplacements[7], beamDisplacements[8] });
        //        }
        //    }

        //    using (StreamWriter outputFile = new StreamWriter(Path.Combine(Directory.GetCurrentDirectory(), $"{_filename}_Paraview_Cohesive.vtu")))
        //    {
        //        outputFile.WriteLine("<?xml version=\"1.0\"?>");
        //        outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">");
        //        outputFile.WriteLine("<UnstructuredGrid>");

        //        outputFile.WriteLine($"<Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
        //        outputFile.WriteLine($"<PointData Vectors=\"U\">");
        //        outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"U\" format=\"ascii\" NumberOfComponents=\"3\">");

        //        for (int i = 0; i < numberOfPoints; i++)
        //        {
        //            var node = nodes[i];
        //            var dx = nodalDisplacements[node][0];
        //            var dy = nodalDisplacements[node][1];
        //            var dz = nodalDisplacements[node][2];

        //            outputFile.WriteLine($"{dx} {dy} {dz}");
        //        }


        //        //for (int i = 2; i < (2 * numberOfPoints); i = i + 4)
        //        //{
        //        //    for (int j = 0; j <= 1; j++)
        //        //    {
        //        //        var node = nodes[i + j];
        //        //        var dx = nodalDisplacements[node][0];
        //        //        var dy = nodalDisplacements[node][1];
        //        //        var dz = nodalDisplacements[node][2];

        //        //        outputFile.WriteLine($"{dx} {dy} {dz}");
        //        //    }                    
        //        //}

        //        outputFile.WriteLine("</DataArray>");

        //        outputFile.WriteLine("</PointData>");
        //        outputFile.WriteLine("<Points>");
        //        outputFile.WriteLine("<DataArray type=\"Float32\" NumberOfComponents=\"3\">");

        //        for (int i = 0; i < numberOfPoints; i++)
        //            outputFile.WriteLine($"{nodes[i].X} {nodes[i].Y} {nodes[i].Z}");

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("</Points>");
        //        outputFile.WriteLine("<Cells>");
        //        outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"connectivity\">");

        //        //for (int i = 0; i < numberOfCells; i++)
        //        //{
        //        //    for (int j = 0; j < 2; j++)
        //        //    {
        //        //        var elementNode = elements[i].Nodes[conn[j]];
        //        //        var indexOfElementNode = nodes.IndexOf(elementNode);
        //        //        outputFile.Write($"{indexOfElementNode} ");
        //        //    }
        //        //    outputFile.WriteLine("");
        //        //}

        //        for (int i = 0; i < numberOfCells; i++)
        //        {
        //            foreach (var node in elements[i].Nodes.Where(x => x.ID > (numberOfHostNodes + numberOfBeamNodes)))
        //            {
        //                var indexOfElementNode = nodes.IndexOf(node);
        //                outputFile.Write($"{indexOfElementNode} ");
        //                outputFile.WriteLine("");
        //            }
                    
        //            outputFile.WriteLine("");
        //        }

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"offsets\">");

        //        var offset = 0;
        //        for (int i = 0; i < numberOfCells; i++)
        //        {
        //            offset += numberOfVerticesPerCell;
        //            outputFile.WriteLine(offset);
        //        }

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"types\">");

        //        for (int i = 0; i < numberOfCells; i++)
        //            outputFile.WriteLine(paraviewCellCode);

        //        outputFile.WriteLine("</DataArray>");
        //        outputFile.WriteLine("</Cells>");
        //        outputFile.WriteLine("</Piece>");
        //        outputFile.WriteLine("</UnstructuredGrid>");
        //        outputFile.WriteLine("</VTKFile>");
        //    }
        //}
    }
}
