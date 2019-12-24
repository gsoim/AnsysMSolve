using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Continuum;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Readers
{
    /// <summary>
	/// Reader for custom isogeometric model files.
	/// </summary>
	public class IsogeometricReader
	{
		private readonly string _filename;
		private int patchID = -1;
        private int NumberOfDimensions;
        private double Thickness;
        private int NumberOfPatches;
        private IContinuumMaterial2D _material2D;
        private IContinuumMaterial3D _material3D;
        private Dictionary<int, int> DegreeKsiDictionary= new Dictionary<int, int>();
        private Dictionary<int, int> DegreeHetaDictionary= new Dictionary<int, int>();
        private Dictionary<int, int> DegreeZetaDictionary= new Dictionary<int, int>();
        private Dictionary<int, int> NumberOfControlPointsKsiDictionary=new Dictionary<int, int>();
        private Dictionary<int, int> NumberOfControlPointsHetaDictionary=new Dictionary<int, int>();
        private Dictionary<int, int> NumberOfControlPointsZetaDictionary=new Dictionary<int, int>();
        private Dictionary<int, double[]> KnotValueVectorsKsiDictionary=new Dictionary<int, double[]>();
        private Dictionary<int, double[]> KnotValueVectorsHetaDictionary= new Dictionary<int, double[]>();
        private Dictionary<int, double[]> KnotValueVectorsZetaDictionary=new Dictionary<int, double[]>();
        private Dictionary<int, int[]> ControlPointIDsDictionary=new Dictionary<int, int[]>();
        private int NumberOfControlPoints;
        private Dictionary<int, ControlPoint> ControlPointsDictionary=new Dictionary<int, ControlPoint>();

        public int NumberOfInterfaces { get; private set; }

        /// <summary>
        /// Defines a custom isogeometric file reader.
        /// </summary>
        /// <param name="modelCreator">An isogeometric <see cref="Model"/>.</param>
        /// <param name="filename">The name of the file to be read.</param>
        public IsogeometricReader(string filename, IContinuumMaterial2D material2D=null, IContinuumMaterial3D material3D=null)
		{
			_filename = filename;
            _material2D = material2D;
            _material3D = material3D;
        }

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

        public Model CreateModelFromFile()
        {
            ReadModelFromFile();
            return CreatePatchData();
        }
		
		private void ReadModelFromFile()
		{
			char[] delimeters = { ' ', '=', '\t' };
			Attributes? name = null;

			int numberOfValues = 0;
			int[] localControlPointIDs;

			var text = System.IO.File.ReadAllLines(_filename);

			for (var i = 0; i < text.Length; i++)
			{
				var line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
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
						NumberOfDimensions = int.Parse(line[1]);
						break;

					case Attributes.thickness:
						Thickness = double.Parse(line[1], CultureInfo.InvariantCulture);
						break;

					case Attributes.numberofpatches:
						NumberOfPatches = int.Parse(line[1]);
						break;

					case Attributes.numberofinterfaces:
						NumberOfInterfaces = int.Parse(line[1]);
						break;

					case Attributes.material:
						break;

					case Attributes.patchid:
						patchID = int.Parse(line[1]);
						break;

					case Attributes.degreeksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Ksi of a patch must be defined after the patchID");
						DegreeKsiDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.degreeheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Heta of a patch must be defined after the patchID");
						DegreeHetaDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.degreezeta:
						if (NumberOfDimensions == 2)
							throw new ArgumentOutOfRangeException("You must not define degree Zeta in case of 2D");
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Zeta of a patch must be defined after the patchID");
						DegreeZetaDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.numberofcpksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Ksi of a patch must be defined after the patchID");
						NumberOfControlPointsKsiDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.numberofcpheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Heta of a patch must be defined after the patchID");
						NumberOfControlPointsHetaDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.numberofcpzeta:
						if (NumberOfDimensions == 2)
							throw new ArgumentOutOfRangeException("You must not define number of Control Points Zeta in case of 2D");
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Zeta of a patch must be defined after the patchID");
						NumberOfControlPointsZetaDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.knotvaluevectorksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Ksi of a patch must be defined after the patchID");
						if (DegreeKsiDictionary[patchID] == 0 || NumberOfControlPointsKsiDictionary[patchID] == 0)
						{
							throw new ArgumentOutOfRangeException("Degree Ksi and number of Control Points per axis Ksi must be defined before Knot Value Vector Ksi.");
						}
						numberOfValues = NumberOfControlPointsKsiDictionary[patchID] + DegreeKsiDictionary[patchID] + 1;
						var KnotValueVectorKsi = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
						{
							KnotValueVectorKsi[j] = double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						}
						KnotValueVectorsKsiDictionary.Add(patchID, KnotValueVectorKsi);
						break;

					case Attributes.knotvaluevectorheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Heta of a patch must be defined after the patchID");
						if (DegreeHetaDictionary[patchID] == 0 || NumberOfControlPointsHetaDictionary[patchID] == 0)
						{
							throw new ArgumentOutOfRangeException("Degree Heta and number of Control Points per axis Heta must be defined before Knot Value Vector Heta.");
						}
						numberOfValues = NumberOfControlPointsHetaDictionary[patchID] + DegreeHetaDictionary[patchID] + 1;
						double[] KnotValueVectorHeta = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
						{
							KnotValueVectorHeta[j] = double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						}
						KnotValueVectorsHetaDictionary.Add(patchID, KnotValueVectorHeta);
						break;

					case Attributes.knotvaluevectorzeta:
						if (NumberOfDimensions == 2)
							throw new ArgumentOutOfRangeException("You must not define Knot Value Vector Zeta in case of 2D");
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Zeta of a patch must be defined after the patchID");
						if (DegreeZetaDictionary[patchID] == 0 || NumberOfControlPointsZetaDictionary[patchID] == 0)
						{
							throw new ArgumentOutOfRangeException("Degree Zeta and number of Control Points per axis Zeta must be defined before Knot Value Vector Zeta.");
						}
						numberOfValues = NumberOfControlPointsZetaDictionary[patchID] + DegreeZetaDictionary[patchID] + 1;
						double[] KnotValueVectorZeta = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
						{
							KnotValueVectorZeta[j] = double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						}
						KnotValueVectorsZetaDictionary.Add(patchID, KnotValueVectorZeta);
						break;

					case Attributes.patchcpid:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Control Points ID of a patch must be defined after the patchID");
						int numberOfPatchCP = (NumberOfDimensions == 3) ? NumberOfControlPointsKsiDictionary[patchID] *
							NumberOfControlPointsHetaDictionary[patchID] *
							NumberOfControlPointsZetaDictionary[patchID] :
							NumberOfControlPointsKsiDictionary[patchID] *
							NumberOfControlPointsHetaDictionary[patchID];
						localControlPointIDs = new int[numberOfPatchCP];
						for (int j = 0; j < numberOfPatchCP; j++)
						{
							localControlPointIDs[j] = int.Parse(line[j + 1]);
						}
						ControlPointIDsDictionary.Add(patchID, localControlPointIDs);
						break;

					case Attributes.cpcoord:
						NumberOfControlPoints = Int32.Parse(line[1]);
						for (int j = 0; j < NumberOfControlPoints; j++)
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
							ControlPointsDictionary.Add(controlPointGlobalID, controlPoint);
						}
						break;

					case Attributes.end:
						return;
				}
			}
		}



        private static List<Knot> CreateKnots2D(Vector singleKnotValuesKsi, Vector singleKnotValuesHeta)
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

        private static List<Knot> CreateKnots3D(Vector singleKnotValuesKsi, Vector singleKnotValuesHeta, Vector singleKnotValuesZeta)
        {
            List<Knot> knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    for (int k = 0; k < singleKnotValuesZeta.Length; k++)
                    {
                        knots.Add(new Knot()
                        {
                            ID = id,
                            Ksi = singleKnotValuesKsi[i],
                            Heta = singleKnotValuesHeta[j],
                            Zeta = singleKnotValuesZeta[k]
                        });
                        id++;
                    }
                }
            }

            return knots;
        }


        private void BuildEdgesDictionary(Model model)
        {
            if (this.NumberOfDimensions == 2)
            {
                CreateEdges2D(model);
            }
            else
            {
                CreateEdges3D(model);
            }
        }

        private void BuildFacesDictionary(Model model)
        {
            if (this.NumberOfDimensions == 3)
            {
                CreateFaces3D(model);
            }
        }

        private void CreateEdges2D(Model model)
        {
            #region EdgeLeft

            var edgeLeft = new Edge(NumberOfControlPointsHetaDictionary[0],1)
            {
                ID = 0,
                Degree = this.DegreeHetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsHetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            int counter = 0;
            for (int i = 0; i < NumberOfControlPointsHetaDictionary[0]; i++)
            {
                edgeLeft.ControlPointsDictionary.Add(counter++, ControlPointsDictionary[i]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(0, edgeLeft);

            #endregion 

            #region EdgeRight

            var edgeRight = new Edge(NumberOfControlPointsHetaDictionary[0],1)
            {
                ID = 1,
                Degree = DegreeHetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsHetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsHetaDictionary[0]; i++)
            {
                edgeRight.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[
                        i + NumberOfControlPointsHetaDictionary[0] * (NumberOfControlPointsKsiDictionary[0] - 1)]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(1, edgeRight);

            #endregion EdgeLeft

            #region EdgeBottom

            Edge edgeBottom = new Edge(NumberOfControlPointsHetaDictionary[0],1)
            {
                ID = 2,
                Degree = DegreeKsiDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsKsiDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                edgeBottom.ControlPointsDictionary.Add(counter++, ControlPointsDictionary[i * NumberOfControlPointsHetaDictionary[0]]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(2, edgeBottom);

            #endregion EdgeBottom

            #region EdgeUp

            Edge edgeUp = new Edge(NumberOfControlPointsHetaDictionary[0],1)
            {
                ID = 3,
                Degree = DegreeKsiDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsKsiDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                edgeUp.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsHetaDictionary[0] + NumberOfControlPointsHetaDictionary[0] - 1]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(3, edgeUp);

            #endregion EdgeUp
        }


        private void CreateEdges3D(Model model)
        {
            #region Edge1

            Edge edge1 = new Edge(NumberOfControlPointsHetaDictionary[0], NumberOfControlPointsZetaDictionary[0])
            {
                ID = 0,
                Degree = DegreeZetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsZetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            int counter = 0;
            for (int i = 0; i < NumberOfControlPointsZetaDictionary[0]; i++)
            {
                edge1.ControlPointsDictionary.Add(counter++, ControlPointsDictionary[i]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(0, edge1);

            #endregion Edge1

            #region Edge2

            Edge edge2 = new Edge(NumberOfControlPointsHetaDictionary[0], NumberOfControlPointsZetaDictionary[0])
            {
                ID = 1,
                Degree = DegreeZetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsZetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsZetaDictionary[0]; i++)
            {
                edge2.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i + NumberOfControlPointsZetaDictionary[0]* (NumberOfControlPointsHetaDictionary[0]- 1)]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(1, edge2);

            #endregion Edge2

            #region Edge3

            Edge edge3 = new Edge(NumberOfControlPointsHetaDictionary[0], NumberOfControlPointsZetaDictionary[0])
            {
                ID = 2,
                Degree = DegreeHetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsHetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsHetaDictionary[0]; i++)
            {
                edge3.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsZetaDictionary[0]]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(2, edge3);

            #endregion Edge3

            #region Edge4

            Edge edge4 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0])
            {
                ID = 3,
                Degree = DegreeHetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsHetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsHetaDictionary[0]; i++)
            {
                edge4.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsZetaDictionary[0]+ NumberOfControlPointsZetaDictionary[0]- 1]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(3, edge4);

            #endregion Edge4

            #region Edge5

            Edge edge5 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0])
            {
                ID = 4,
                Degree = DegreeZetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsZetaDictionary[0]
            };
            int offset = NumberOfControlPointsZetaDictionary[0]* NumberOfControlPointsHetaDictionary[0]*
                         (NumberOfControlPointsKsiDictionary[0] - 1);
            edge5.Patch = model.PatchesDictionary[0];
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsZetaDictionary[0]; i++)
            {
                edge5.ControlPointsDictionary.Add(counter++, ControlPointsDictionary[i + offset]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(4, edge5);

            #endregion Edge5

            #region Edge6

            Edge edge6 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0])
            {
                ID = 5,
                Degree = DegreeZetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsZetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsZetaDictionary[0]; i++)
            {
                edge6.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i + NumberOfControlPointsZetaDictionary[0] * (NumberOfControlPointsHetaDictionary[0]- 1) + offset]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(5, edge6);

            #endregion Edge6

            #region Edge7

            Edge edge7 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0])
            {
                ID = 6,
                Degree = DegreeHetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsHetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsHetaDictionary[0]; i++)
            {
                edge7.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsZetaDictionary[0]+ offset]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(6, edge7);

            #endregion Edge7

            #region Edge8

            Edge edge8 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0])
            {
                ID = 7,
                Degree = DegreeHetaDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsHetaDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsHetaDictionary[0]; i++)
            {
                edge8.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsZetaDictionary[0]+ NumberOfControlPointsZetaDictionary[0]- 1 + offset]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(7, edge8);

            #endregion Edge8

            #region Edge9

            Edge edge9 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0])
            {
                ID = 8,
                Degree = DegreeKsiDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsKsiDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                edge9.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsZetaDictionary[0]* NumberOfControlPointsHetaDictionary[0]]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(8, edge9);

            #endregion Edge9

            #region Edge10

            Edge edge10 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0])
            {
                ID = 9,
                Degree = DegreeKsiDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsKsiDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                edge10.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsZetaDictionary[0]* NumberOfControlPointsHetaDictionary[0]
                                  + NumberOfControlPointsZetaDictionary[0]- 1]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(9, edge10);

            #endregion Edge10

            #region Edge11

            Edge edge11 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0])
            {
                ID = 10,
                Degree = DegreeKsiDictionary[0],
                KnotValueVector = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]),
                NumberOfControlPoints = NumberOfControlPointsKsiDictionary[0],
                Patch = model.PatchesDictionary[0]
            };
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                edge11.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsZetaDictionary[0]* NumberOfControlPointsHetaDictionary[0]+
                                  NumberOfControlPointsZetaDictionary[0]* (NumberOfControlPointsHetaDictionary[0]- 1)]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(10, edge11);

            #endregion Edge11

            #region Edge12

            Edge edge12 = new Edge(NumberOfControlPointsHetaDictionary[0],
                NumberOfControlPointsZetaDictionary[0]);
            edge12.ID = 11;
            edge12.Degree = DegreeKsiDictionary[0];
            edge12.KnotValueVector = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]);
            edge12.NumberOfControlPoints = NumberOfControlPointsKsiDictionary[0];
            edge12.Patch = model.PatchesDictionary[0];
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                edge12.ControlPointsDictionary.Add(counter++,
                    ControlPointsDictionary[i * NumberOfControlPointsZetaDictionary[0]* NumberOfControlPointsHetaDictionary[0]+
                                  NumberOfControlPointsZetaDictionary[0]* (NumberOfControlPointsHetaDictionary[0]- 1) +
                                  NumberOfControlPointsZetaDictionary[0]- 1]);
            }

            model.PatchesDictionary[0].EdgesDictionary.Add(11, edge12);

            #endregion Edge12
        }

        private void CreateElements2D(Model model, Vector singleKnotValuesKsi, Vector singleKnotValuesHeta, List<Knot> knots)
        {
            Vector multiplicityKsi = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]).RemoveDuplicatesFindMultiplicity()[1];
            Vector multiplicityHeta = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]).RemoveDuplicatesFindMultiplicity()[1];

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
                    if (multiplicityKsi[i + 1] - DegreeKsiDictionary[0]> 0)
                    {
                        multiplicityElementKsi = (int)multiplicityKsi[i + 1] - DegreeKsiDictionary[0];
                    }

                    int multiplicityElementHeta = 0;
                    if (multiplicityHeta[j + 1] - DegreeHetaDictionary[0]> 0)
                    {
                        multiplicityElementHeta = (int)multiplicityHeta[j + 1] - DegreeHetaDictionary[0];
                    }

                    int nurbsSupportKsi = DegreeKsiDictionary[0] + 1;
                    int nurbsSupportHeta = DegreeHetaDictionary[0]+ 1;

                    IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                    for (int k = 0; k < nurbsSupportKsi; k++)
                    {
                        for (int l = 0; l < nurbsSupportHeta; l++)
                        {
                            int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHetaDictionary[0]+
                                                 (j + multiplicityElementHeta) + k * NumberOfControlPointsHetaDictionary[0]+ l;
                            elementControlPoints.Add(this.ControlPointsDictionary[controlPointID]);
                        }
                    }

                    int elementID = i * numberOfElementsHeta + j;


                    var gauss = new GaussQuadrature();
                    var gaussPoints = gauss.CalculateElementGaussPoints(DegreeKsiDictionary[0], 
                        DegreeHetaDictionary[0], knotsOfElement).ToArray();
                    var parametricGaussPointKsi = new double[DegreeKsiDictionary[0] + 1];
                    for (int m = 0; m < DegreeKsiDictionary[0] + 1;m++)
                    {
                        parametricGaussPointKsi[m] = gaussPoints[m * (DegreeHetaDictionary[0] + 1)].Ksi;
                    }

                    var parametricGaussPointHeta = new double[DegreeHetaDictionary[0] + 1];
                    for (int m = 0; m < DegreeHetaDictionary[0] + 1; m++)
                    {
                        parametricGaussPointHeta[m] = gaussPoints[m].Heta;
                    }
                    var nurbs = new Nurbs2D(DegreeKsiDictionary[0], KnotValueVectorsKsiDictionary[0],
                        DegreeHetaDictionary[0], KnotValueVectorsHetaDictionary[0],
                        elementControlPoints.ToArray(), parametricGaussPointKsi, parametricGaussPointHeta);
                    Element element = new Element
                    {
                        ID = elementID,
                        Patch = model.PatchesDictionary[0],
                        ElementType = new ContinuumElement2D(_material2D, nurbs, gaussPoints, Thickness)
                    };
                    element.AddKnots(knotsOfElement);
                    element.AddControlPoints(elementControlPoints.ToList());
                    model.PatchesDictionary[0].Elements.Add(element);
                    model.ElementsDictionary.Add(elementID, element);
                }
            }
        }

        private void CreateElements3D(Model model,List<Knot> knots)
        {
            Vector singlesKnotValuesKsi = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]).RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityKsi = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]).RemoveDuplicatesFindMultiplicity()[1];
            Vector singlesKnotValuesHeta = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]).RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityHeta = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]).RemoveDuplicatesFindMultiplicity()[1];
            Vector singlesKnotValuesZeta = Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]).RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityZeta = Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]).RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singlesKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singlesKnotValuesHeta.Length - 1;
            int numberOfElementsZeta = singlesKnotValuesZeta.Length - 1;

            if (numberOfElementsKsi * numberOfElementsHeta * numberOfElementsZeta == 0)
            {
                throw new ArgumentException("Number of Elements should be defined before Element Connectivity");
            }

            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                for (int j = 0; j < numberOfElementsHeta; j++)
                {
                    for (int k = 0; k < numberOfElementsZeta; k++)
                    {
                        IList<Knot> knotsOfElement = new List<Knot>
                        {
                            knots[
                            i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length +
                            k],
                            knots[
                            i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length +
                            k + 1],
                            knots[
                            i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            (j + 1) * singlesKnotValuesZeta.Length + k],
                            knots[
                            i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            (j + 1) * singlesKnotValuesZeta.Length + k + 1],
                            knots[
                            (i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            j * singlesKnotValuesZeta.Length + k],
                            knots[
                            (i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            j * singlesKnotValuesZeta.Length + k + 1],
                            knots[
                            (i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            (j + 1) * singlesKnotValuesZeta.Length + k],
                            knots[
                            (i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            (j + 1) * singlesKnotValuesZeta.Length + k + 1]
                        };

                        int multiplicityElementKsi = 0;
                        if (multiplicityKsi[i + 1] - this.DegreeKsiDictionary[0] > 0)
                        {
                            multiplicityElementKsi = (int)multiplicityKsi[i + 1] - DegreeKsiDictionary[0];
                        }

                        int multiplicityElementHeta = 0;
                        if (multiplicityHeta[j + 1] - this.DegreeHetaDictionary[0] > 0)
                        {
                            multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHetaDictionary[0];
                        }

                        int multiplicityElementZeta = 0;
                        if (multiplicityZeta[k + 1] - this.DegreeZetaDictionary[0] > 0)
                        {
                            multiplicityElementZeta = (int)multiplicityZeta[k + 1] - this.DegreeZetaDictionary[0];
                        }

                        int nurbsSupportKsi = this.DegreeKsiDictionary[0] + 1;
                        int nurbsSupportHeta = this.DegreeHetaDictionary[0] + 1;
                        int nurbsSupportZeta = this.DegreeZetaDictionary[0] + 1;

                        IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                        for (int l = 0; l < nurbsSupportKsi; l++)
                        {
                            for (int m = 0; m < nurbsSupportHeta; m++)
                            {
                                for (int n = 0; n < nurbsSupportZeta; n++)
                                {
                                    int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHetaDictionary[0]*
                                                         NumberOfControlPointsZetaDictionary[0]+
                                                         (j + multiplicityElementHeta) * NumberOfControlPointsZetaDictionary[0]+
                                                         (k + multiplicityElementZeta) +
                                                         l * NumberOfControlPointsHetaDictionary[0]* NumberOfControlPointsZetaDictionary[0]+
                                                         m * NumberOfControlPointsZetaDictionary[0]+ n;

                                    elementControlPoints.Add(ControlPointsDictionary[controlPointID]);
                                }
                            }
                        }

                        int elementID = i * numberOfElementsHeta * numberOfElementsZeta + j * numberOfElementsZeta + k;
                        var gauss = new GaussQuadrature();
                        var gaussPoints = gauss.CalculateElementGaussPoints(DegreeKsiDictionary[0],
                            DegreeHetaDictionary[0], DegreeZetaDictionary[0], knotsOfElement);
                        var nurbs= new Nurbs3D(NumberOfControlPointsKsiDictionary[0], 
                            NumberOfControlPointsHetaDictionary[0],NumberOfControlPointsZetaDictionary[0],
                            DegreeKsiDictionary[0], DegreeHetaDictionary[0],DegreeZetaDictionary[0],
                            KnotValueVectorsKsiDictionary[0], KnotValueVectorsHetaDictionary[0],
                            KnotValueVectorsZetaDictionary[0], elementControlPoints.ToArray(), gaussPoints);

                        Element element = new Element
                        {
                            ID = elementID,
                            Patch = model.PatchesDictionary[0],
                            ElementType = new ContinuumElement3D(_material3D, nurbs, gaussPoints)
                        };
                        element.AddKnots(knotsOfElement);
                        element.AddControlPoints(elementControlPoints.ToList());
                        model.PatchesDictionary[0].Elements.Add(element);
                        model.PatchesDictionary[0].Elements.Add(element);
                        model.ElementsDictionary.Add(elementID, element);
                    }
                }
            }
        }

        private void CreateFaces3D(Model model)
        {
            #region FaceLeft

            Face faceLeft = new Face();
            faceLeft.Degrees[0] = DegreeHetaDictionary[0];
            faceLeft.Degrees[1] = DegreeZetaDictionary[0];
            faceLeft.KnotValueVectors.Add(0, Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]));
            faceLeft.KnotValueVectors.Add(1, Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]));
            faceLeft.Patch = model.PatchesDictionary[0];
            int counter = 0;
            for (int i = 0; i < NumberOfControlPointsHetaDictionary[0]; i++)
            {
                for (int j = 0; j < NumberOfControlPointsZetaDictionary[0]; j++)
                {
                    faceLeft.ControlPointsDictionary.Add(counter++,
                        ControlPointsDictionary[j + NumberOfControlPointsZetaDictionary[0]* i]);
                }
            }

            model.PatchesDictionary[0].FacesDictionary.Add(0, faceLeft);

            #endregion 

            #region FaceLeft

            Face faceRight = new Face();
            faceRight.Degrees[0] = DegreeHetaDictionary[0];
            faceRight.Degrees[1] = DegreeZetaDictionary[0];
            faceRight.KnotValueVectors.Add(0, Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]));
            faceRight.KnotValueVectors.Add(1, Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]));
            faceRight.Patch = model.PatchesDictionary[0];
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsHetaDictionary[0]; i++)
            {
                for (int j = 0; j < NumberOfControlPointsZetaDictionary[0]; j++)
                {
                    faceRight.ControlPointsDictionary.Add(counter++,
                        ControlPointsDictionary[j + NumberOfControlPointsZetaDictionary[0]* i +
                                      NumberOfControlPointsHetaDictionary[0]* NumberOfControlPointsZetaDictionary[0]*
                                      (NumberOfControlPointsKsiDictionary[0]- 1)]);
                }
            }

            model.PatchesDictionary[0].FacesDictionary.Add(1, faceRight);

            #endregion FaceLeft

            #region FaceBottom

            Face faceBottom = new Face();
            faceBottom.Degrees[0] = DegreeKsiDictionary[0];
            faceBottom.Degrees[1] = DegreeHetaDictionary[0];
            faceBottom.KnotValueVectors.Add(0, Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]));
            faceBottom.KnotValueVectors.Add(1, Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]));
            faceBottom.Patch = model.PatchesDictionary[0];
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                for (int j = 0; j < NumberOfControlPointsHetaDictionary[0]; j++)
                {
                    faceBottom.ControlPointsDictionary.Add(counter++,
                        ControlPointsDictionary[j * NumberOfControlPointsZetaDictionary[0]+
                                      i * NumberOfControlPointsZetaDictionary[0]* NumberOfControlPointsHetaDictionary[0]]);
                }
            }

            model.PatchesDictionary[0].FacesDictionary.Add(2, faceBottom);

            #endregion FaceBottom

            #region FaceUp

            Face faceUp = new Face();
            faceUp.Degrees[0] = DegreeKsiDictionary[0];
            faceUp.Degrees[1] = DegreeHetaDictionary[0];
            faceUp.KnotValueVectors.Add(0, Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]));
            faceUp.KnotValueVectors.Add(1, Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]));
            faceUp.Patch = model.PatchesDictionary[0];
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                for (int j = 0; j < NumberOfControlPointsHetaDictionary[0]; j++)
                {
                    faceUp.ControlPointsDictionary.Add(counter++,
                        ControlPointsDictionary[NumberOfControlPointsZetaDictionary[0]- 1 + j * NumberOfControlPointsZetaDictionary[0]+
                                      i * NumberOfControlPointsZetaDictionary[0]* NumberOfControlPointsHetaDictionary[0]]);
                }
            }

            model.PatchesDictionary[0].FacesDictionary.Add(3, faceUp);

            #endregion FaceUp

            #region FaceFront

            Face faceFront = new Face();
            faceFront.Degrees[0] = DegreeKsiDictionary[0];
            faceFront.Degrees[1] = DegreeZetaDictionary[0];
            faceFront.KnotValueVectors.Add(0, Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]));
            faceFront.KnotValueVectors.Add(1, Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]));
            faceFront.Patch = model.PatchesDictionary[0];
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                for (int j = 0; j < NumberOfControlPointsZetaDictionary[0]; j++)
                {
                    faceFront.ControlPointsDictionary.Add(counter++,
                        ControlPointsDictionary[j + i * NumberOfControlPointsHetaDictionary[0]* NumberOfControlPointsZetaDictionary[0]]);
                }
            }

            model.PatchesDictionary[0].FacesDictionary.Add(4, faceFront);

            #endregion FaceFront

            #region FaceBack

            Face faceBack = new Face();
            faceBack.Degrees[0] = DegreeKsiDictionary[0];
            faceBack.Degrees[1] = DegreeZetaDictionary[0];
            faceBack.KnotValueVectors.Add(0, Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]));
            faceBack.KnotValueVectors.Add(1, Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]));
            faceBack.Patch = model.PatchesDictionary[0];
            counter = 0;
            for (int i = 0; i < NumberOfControlPointsKsiDictionary[0]; i++)
            {
                for (int j = 0; j < NumberOfControlPointsZetaDictionary[0]; j++)
                {
                    faceBack.ControlPointsDictionary.Add(counter++,
                        ControlPointsDictionary[j + i * NumberOfControlPointsHetaDictionary[0]* NumberOfControlPointsZetaDictionary[0]+
                                      NumberOfControlPointsZetaDictionary[0]* (NumberOfControlPointsHetaDictionary[0]- 1)]);
                }
            }

            model.PatchesDictionary[0].FacesDictionary.Add(5, faceBack);

            #endregion FaceBack
        }

        private void CreateNURBSElements2D(Model model)
        {
            #region Knots

            Vector singleKnotValuesKsi = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]).RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesHeta = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]).RemoveDuplicatesFindMultiplicity()[0];

            var knots = CreateKnots2D(singleKnotValuesKsi, singleKnotValuesHeta);

            #endregion Knots

            #region Elements

            CreateElements2D(model, singleKnotValuesKsi, singleKnotValuesHeta, knots);

            #endregion Elements
        }

        private void CreateNURBSElements3D(Model model)
        {
            #region Knots

            Vector singleKnotValuesKsi = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[0]).RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesHeta = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[0]).RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesZeta = Vector.CreateFromArray(KnotValueVectorsZetaDictionary[0]).RemoveDuplicatesFindMultiplicity()[0];

            var knots = CreateKnots3D(singleKnotValuesKsi, singleKnotValuesHeta, singleKnotValuesZeta);

            #endregion Knots

            #region Elements

            CreateElements3D(model, knots);

            #endregion Elements
        }

        private void CreatePatchData2D(Model model)
        {
            model.PatchesDictionary.Add(0, new Patch());
            model.PatchesDictionary[0].NumberOfDimensions = 2;
            foreach (var controlPoint in ControlPointsDictionary)
                model.PatchesDictionary[0].ControlPoints.Add(controlPoint.Value);
            CreateNURBSElements2D(model);
            BuildEdgesDictionary(model);
        }

        private void CreatePatchData3D(Model model)
        {
            model.PatchesDictionary.Add(0, new Patch());
            model.PatchesDictionary[0].NumberOfDimensions = 3;
            foreach (var controlPoint in ControlPointsDictionary)
                model.PatchesDictionary[0].ControlPoints.Add(controlPoint.Value);
            CreateNURBSElements3D(model);
            BuildEdgesDictionary(model);
            BuildFacesDictionary(model);
        }

        public Model CreatePatchData()
        {
            var model= new Model();
            foreach (var controlPoint in ControlPointsDictionary)
            {
                model.ControlPointsDictionary.Add(controlPoint.Key, controlPoint.Value);
            }

            if (this.NumberOfDimensions == 2)
                CreatePatchData2D(model);
            else
                CreatePatchData3D(model);

            return model;
        }
    }


}
