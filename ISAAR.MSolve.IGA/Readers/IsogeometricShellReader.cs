using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Structural;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Readers
{
    public enum GeometricalFormulation
	{
		Linear,
		NonLinear,
		SectionNonLinear
	}

	/// <summary>
	/// Reader for custom isogeometric shell model files.
	/// </summary>
	public class IsogeometricShellReader
	{
		private readonly string _filename;
        private readonly GeometricalFormulation _formulation;
        private readonly IShellMaterial _material;
        private readonly IShellSectionMaterial _sectionaMaterial;

        private enum Attributes
		{
			numberofdimensions,
			numberofpatches, numberofinterfaces,
			patchid, interfaceid,
			degreeksi, degreeheta, degreezeta,
			numberofcpksi, numberofcpheta, numberofcpzeta,
			knotvaluevectorksi, knotvaluevectorheta, knotvaluevectorzeta,
			patchcpid,
			cpcoord, thickness, material, end
		}

		/// <summary>
		/// Defines a custom isogeometric shell file reader.
		/// </summary>
		/// <param name="modelCreator">An <see cref="ModelCreator"/> object responsible for generating the model.</param>
		/// <param name="filename">The name of the file to be read.</param>
		public IsogeometricShellReader(GeometricalFormulation formulation,
            string filename, IShellMaterial material=null,IShellSectionMaterial sectionMaterial=null )
		{
			_filename = filename;
            _formulation = formulation;
            _material = material;
            _sectionaMaterial = sectionMaterial;
        }

		private Dictionary<int, int[]> ControlPointIDsDictionary = new Dictionary<int, int[]>();
        private int NumberOfDimensions;
        private double Thickness;
        private int DegreeKsi;
        private int DegreeHeta;
        private int NumberOfControlPointsKsi;
        private int NumberOfControlPointsHeta;
        private Vector KnotValueVectorKsi;
        private Vector KnotValueVectorHeta;
        private Dictionary<int, ControlPoint> ControlPointsDictionary=new Dictionary<int, ControlPoint>();

        public Model GenerateModelFromFile()
        {
            ReadShellModelFromFile();
            return GenerateModel();
        }

        private void ReadShellModelFromFile()
		{
			char[] delimeters = { ' ', '=', '\t' };
			Attributes? name = null;

			int patchID = -1;
			int numberOfValues = 0;
			int[] localControlPointIDs;
			int counterElementID = 0;
			int counterCPID;

			string[] text = File.ReadAllLines(_filename);

			for (int i = 0; i < text.Length; i++)
			{
				string[] line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
				if (line.Length == 0)
				{
					continue;
				}
				try
				{
					name = (Attributes)Enum.Parse(typeof(Attributes), line[0].ToLower());
				}
				catch (Exception exception)
				{
					throw new KeyNotFoundException($"Variable name {line[0]} is not found. {exception.Message}");
				}

				switch (name)
				{
					case Attributes.numberofdimensions:
						this.NumberOfDimensions = Int32.Parse(line[1]);
						break;

					case Attributes.thickness:
						this.Thickness = Double.Parse(line[1], CultureInfo.InvariantCulture);
						break;

					case Attributes.numberofpatches:
						break;

					case Attributes.patchid:
						patchID = Int32.Parse(line[1]);
						break;

					case Attributes.degreeksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Ksi of a patch must be defined after the patchID");
						this.DegreeKsi = Int32.Parse(line[1]);
						break;

					case Attributes.degreeheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Heta of a patch must be defined after the patchID");
						this.DegreeHeta = Int32.Parse(line[1]);
						break;

					case Attributes.numberofcpksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Ksi of a patch must be defined after the patchID");
						this.NumberOfControlPointsKsi = Int32.Parse(line[1]);
						break;

					case Attributes.numberofcpheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Heta of a patch must be defined after the patchID");
						this.NumberOfControlPointsHeta = Int32.Parse(line[1]);
						break;

					case Attributes.knotvaluevectorksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Ksi of a patch must be defined after the patchID");
						if (this.DegreeKsi == 0 || this.NumberOfControlPointsKsi == 0)
							throw new ArgumentOutOfRangeException("Degree Ksi and number of Control Points per axis Ksi must be defined before Knot Value Vector Ksi.");
						numberOfValues = this.DegreeKsi + this.NumberOfControlPointsKsi + 1;
						double[] KnotValueVectorKsi = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
							KnotValueVectorKsi[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						this.KnotValueVectorKsi = Vector.CreateFromArray(KnotValueVectorKsi);
						break;

					case Attributes.knotvaluevectorheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Heta of a patch must be defined after the patchID");
						if (this.DegreeHeta == 0 || this.NumberOfControlPointsHeta == 0)
							throw new ArgumentOutOfRangeException("Degree Heta and number of Control Points per axis Heta must be defined before Knot Value Vector Heta.");
						numberOfValues = this.DegreeHeta + this.NumberOfControlPointsHeta + 1;
						double[] KnotValueVectorHeta = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
							KnotValueVectorHeta[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						this.KnotValueVectorHeta = Vector.CreateFromArray(KnotValueVectorHeta);
						break;

					case Attributes.patchcpid:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Control Points ID of a patch must be defined after the patchID");
						int numberOfPatchCP = this.NumberOfControlPointsKsi *
											  this.NumberOfControlPointsHeta;
						localControlPointIDs = new int[numberOfPatchCP];
						for (int j = 0; j < numberOfPatchCP; j++)
						{
							localControlPointIDs[j] = Int32.Parse(line[j + 1]);
						}
						ControlPointIDsDictionary.Add(patchID, localControlPointIDs);
						break;

					case Attributes.cpcoord:
						var numberOfControlPoints = Int32.Parse(line[1]);
						for (int j = 0; j < numberOfControlPoints; j++)
						{
							i++;
							line = text[i].Split(delimeters);
							int controlPointGlobalID = Int32.Parse(line[0]);
							double x = double.Parse(line[1], CultureInfo.InvariantCulture);
							double y = double.Parse(line[2], CultureInfo.InvariantCulture);
							double z = double.Parse(line[3], CultureInfo.InvariantCulture);
							double w = double.Parse(line[4], CultureInfo.InvariantCulture);
							ControlPoint controlPoint = new ControlPoint()
							{ ID = controlPointGlobalID, X = x, Y = y, Z = z, WeightFactor = w };
							this.ControlPointsDictionary.Add(controlPointGlobalID, controlPoint);
						}
						break;

					case Attributes.end:
						return;
				}
			}
		}

        private Model GenerateModel()
        {
            var model= new Model();
            model.PatchesDictionary.Add(0, new Patch());
            foreach (var controlPoint in ControlPointsDictionary)
            {
                model.ControlPointsDictionary.Add(controlPoint.Key, controlPoint.Value);
                model.PatchesDictionary[0].ControlPoints.Add(controlPoint.Value);
            }

            CreateEdges2D(model);
            CreateNURBSShells(model, _formulation);
            var counterElementID = 0;
            foreach (var element in model.PatchesDictionary[0].Elements)
            {
                model.ElementsDictionary.Add(counterElementID, element);
                counterElementID++;
            }
            return model;
        }



        private void CreateNURBSShells(Model model, GeometricalFormulation formulation)
        {
            #region Knots

            Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];

            var knots = CreateShellKnots(singleKnotValuesKsi, singleKnotValuesHeta);

            #endregion Knots

            #region Elements

            CreateShellElements(model,singleKnotValuesKsi, singleKnotValuesHeta, knots, formulation);

            #endregion Elements
        }


        private void CreateShellElements(Model model,Vector singleKnotValuesKsi, Vector singleKnotValuesHeta, List<Knot> knots, GeometricalFormulation formulation)
        {
            Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
            Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
            if (numberOfElementsKsi * numberOfElementsHeta == 0)
            {
                throw new ArgumentException("Number of Elements should be defined before Element Connectivity");
            }

            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                for (int j = 0; j < numberOfElementsHeta; j++)
                {
                    IList<Knot> knotsOfElement = new List<Knot>
                    {
                        knots[i * singleKnotValuesHeta.Length + j],
                        knots[i * singleKnotValuesHeta.Length + j + 1],
                        knots[(i + 1) * singleKnotValuesHeta.Length + j],
                        knots[(i + 1) * singleKnotValuesHeta.Length + j + 1]
                    };

                    int multiplicityElementKsi = 0;
                    if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
                    {
                        multiplicityElementKsi = (int)multiplicityKsi[i + 1] - this.DegreeKsi;
                    }

                    int multiplicityElementHeta = 0;
                    if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
                    {
                        multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
                    }

                    int nurbsSupportKsi = this.DegreeKsi + 1;
                    int nurbsSupportHeta = this.DegreeHeta + 1;

                    IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                    for (int k = 0; k < nurbsSupportKsi; k++)
                    {
                        for (int l = 0; l < nurbsSupportHeta; l++)
                        {
                            int controlPointID = (i + multiplicityElementKsi) * this.NumberOfControlPointsHeta +
                                                 (j + multiplicityElementHeta) + k * this.NumberOfControlPointsHeta + l;
                            elementControlPoints.Add(model.ControlPointsDictionary[controlPointID]);
                        }
                    }

                    int elementID = i * numberOfElementsHeta + j;
                    Element element = null;

                    var gauss= new GaussQuadrature();
                    var gaussPoints = gauss.CalculateElementGaussPoints(DegreeKsi, DegreeHeta, knotsOfElement).ToArray();
                    var parametricGaussPointKsi = new double[DegreeKsi + 1];
                    for (int m = 0; m < DegreeKsi + 1; m++)
                    {
                        parametricGaussPointKsi[m] = gaussPoints[m * (DegreeHeta + 1)].Ksi;
                    }

                    var parametricGaussPointHeta = new double[DegreeHeta + 1];
                    for (int m = 0; m < DegreeHeta + 1; m++)
                    {
                        parametricGaussPointHeta[m] = gaussPoints[m].Heta;
                    }
                    var nurbs = new Nurbs2D(DegreeKsi, KnotValueVectorKsi.CopyToArray(),
                        DegreeHeta, KnotValueVectorHeta.CopyToArray(),
                        elementControlPoints.ToArray(), parametricGaussPointKsi, parametricGaussPointHeta);
                    switch (formulation)
                    {
                        case GeometricalFormulation.Linear:
                            element = new Element
                            {
                                ID = elementID,
                                Patch = model.PatchesDictionary[0],
                                ElementType = new KirchhoffLoveShellElement(_sectionaMaterial, nurbs, gaussPoints, Thickness)
                            };
                            element.AddKnots(knotsOfElement);
                            element.AddControlPoints(elementControlPoints);
                            break;
                        case GeometricalFormulation.NonLinear:
                            element = new NurbsKirchhoffLoveShellElementNL(_material,
                                knotsOfElement,nurbs, elementControlPoints, model.PatchesDictionary[0],Thickness, DegreeKsi, DegreeHeta)
                            {
                                ID = elementID,
                                ElementType = new NurbsKirchhoffLoveShellElementNL(_material,
                                    knotsOfElement, nurbs, elementControlPoints, model.PatchesDictionary[0], Thickness, DegreeKsi, DegreeHeta)
                                {
                                    ID = elementID,
                                }
                            };
                            break;
                        case GeometricalFormulation.SectionNonLinear:
                            throw new NotImplementedException();
                            //element = new NurbsKirchhoffLoveShellElementSectionNL(_material,knotsOfElement.ToArray(), elementControlPoints, nurbs, model.PatchesDictionary[0],Thickness,DegreeKsi, DegreeHeta)
                            // {
                            //    ID = elementID,
                            //    ElementType = new NurbsKirchhoffLoveShellElementSectionNL(_material, knotsOfElement.ToArray(), elementControlPoints, nurbs, model.PatchesDictionary[0], Thickness, DegreeKsi, DegreeHeta)
                            //};
                            //break;
                    }
                    model.PatchesDictionary[0].Elements.Add(element);
                }
            }
        }


        private static List<Knot> CreateShellKnots(Vector singleKnotValuesKsi, Vector singleKnotValuesHeta)
        {
            List<Knot> knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = 0.0 });
                    id++;
                }
            }

            return knots;
        }


        private void CreateEdges2D(Model model)
        {
            #region EdgeLeft

            Edge edgeLeft = new Edge(NumberOfControlPointsHeta, 1)
            {
                ID = 0,
                Degree = this.DegreeHeta,
                KnotValueVector = this.KnotValueVectorHeta,
                NumberOfControlPoints = this.NumberOfControlPointsHeta,
                Patch = model.PatchesDictionary[0]
            };
            int counter = 0;
            for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
            {
                edgeLeft.ControlPointsDictionary.Add(counter++, model.ControlPointsDictionary[i]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(0, edgeLeft);

            #endregion 

            #region EdgeRight

            var edgeRight = new Edge(NumberOfControlPointsHeta, 1)
            {
                ID = 1,
                Degree = this.DegreeHeta,
                KnotValueVector = this.KnotValueVectorHeta,
                NumberOfControlPoints = this.NumberOfControlPointsHeta,
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
            {
                edgeRight.ControlPointsDictionary.Add(counter++,
                    model.ControlPointsDictionary[i + this.NumberOfControlPointsHeta * (this.NumberOfControlPointsKsi - 1)]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(1, edgeRight);

            #endregion EdgeLeft

            #region EdgeBottom

            Edge edgeBottom = new Edge(NumberOfControlPointsHeta, 1)
            {
                ID = 2,
                Degree = this.DegreeKsi,
                KnotValueVector = this.KnotValueVectorKsi,
                NumberOfControlPoints = this.NumberOfControlPointsKsi,
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
            {
                edgeBottom.ControlPointsDictionary.Add(counter++, model.ControlPointsDictionary[i * this.NumberOfControlPointsHeta]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(2, edgeBottom);

            #endregion EdgeBottom

            #region EdgeUp

            var edgeUp = new Edge(NumberOfControlPointsHeta, 1)
            {
                ID = 3,
                Degree = this.DegreeKsi,
                KnotValueVector = this.KnotValueVectorKsi,
                NumberOfControlPoints = this.NumberOfControlPointsKsi,
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
            {
                edgeUp.ControlPointsDictionary.Add(counter++,
                    model.ControlPointsDictionary[i * this.NumberOfControlPointsHeta + this.NumberOfControlPointsHeta - 1]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(3, edgeUp);

            #endregion EdgeUp
        }

    }
}
