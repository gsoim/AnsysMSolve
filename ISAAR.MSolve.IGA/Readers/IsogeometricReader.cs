namespace MGroup.IGA.Readers
{
	using System;
	using System.Collections.Generic;
	using System.Globalization;

	using MGroup.IGA.Entities;
	using MGroup.Materials;

	/// <summary>
	/// Reader for custom isogeometric model files.
	/// </summary>
	public class IsogeometricReader
	{
		private readonly string _filename;
		private readonly ModelCreator _modelCreator;
		private int patchID = -1;

		/// <summary>
		/// Defines a custom isogeometric file reader.
		/// </summary>
		/// <param name="modelCreator">An isogeometric <see cref="Model"/>.</param>
		/// <param name="filename">The name of the file to be read.</param>
		public IsogeometricReader(ModelCreator modelCreator, string filename)
		{
			_modelCreator = modelCreator;
			_filename = filename;
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

		/// <summary>
		/// Creates model from custom isogeometric file.
		/// </summary>
		public void CreateModelFromFile()
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
						_modelCreator.NumberOfDimensions = int.Parse(line[1]);
						break;

					case Attributes.thickness:
						_modelCreator.Thickness = double.Parse(line[1], CultureInfo.InvariantCulture);
						break;

					case Attributes.numberofpatches:
						_modelCreator.NumberOfPatches = int.Parse(line[1]);
						break;

					case Attributes.numberofinterfaces:
						_modelCreator.Model.NumberOfInterfaces = int.Parse(line[1]);
						break;

					case Attributes.material:
						if (_modelCreator.NumberOfDimensions == 2)
						{
							_modelCreator.Material = new ElasticMaterial2D((line[4] == "plstress") ? StressState2D.PlaneStress : StressState2D.PlaneStrain) { YoungModulus = Double.Parse(line[2], CultureInfo.InvariantCulture), PoissonRatio = Double.Parse(line[3], CultureInfo.InvariantCulture) };
						}
						else
						{
							_modelCreator.Material = new ElasticMaterial3D { YoungModulus = double.Parse(line[2], CultureInfo.InvariantCulture), PoissonRatio = Double.Parse(line[3], CultureInfo.InvariantCulture) };
						}
						break;

					case Attributes.patchid:
						patchID = int.Parse(line[1]);
						break;

					case Attributes.degreeksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Ksi of a patch must be defined after the patchID");
						_modelCreator.DegreeKsiDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.degreeheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Heta of a patch must be defined after the patchID");
						_modelCreator.DegreeHetaDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.degreezeta:
						if (_modelCreator.NumberOfDimensions == 2)
							throw new ArgumentOutOfRangeException("You must not define degree Zeta in case of 2D");
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Zeta of a patch must be defined after the patchID");
						_modelCreator.DegreeZetaDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.numberofcpksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Ksi of a patch must be defined after the patchID");
						_modelCreator.NumberOfControlPointsKsiDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.numberofcpheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Heta of a patch must be defined after the patchID");
						_modelCreator.NumberOfControlPointsHetaDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.numberofcpzeta:
						if (_modelCreator.NumberOfDimensions == 2)
							throw new ArgumentOutOfRangeException("You must not define number of Control Points Zeta in case of 2D");
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Zeta of a patch must be defined after the patchID");
						_modelCreator.NumberOfControlPointsZetaDictionary.Add(patchID, int.Parse(line[1]));
						break;

					case Attributes.knotvaluevectorksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Ksi of a patch must be defined after the patchID");
						if (_modelCreator.DegreeKsiDictionary[patchID] == 0 || _modelCreator.NumberOfControlPointsKsiDictionary[patchID] == 0)
						{
							throw new ArgumentOutOfRangeException("Degree Ksi and number of Control Points per axis Ksi must be defined before Knot Value Vector Ksi.");
						}
						numberOfValues = _modelCreator.NumberOfControlPointsKsiDictionary[patchID] + _modelCreator.DegreeKsiDictionary[patchID] + 1;
						var KnotValueVectorKsi = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
						{
							KnotValueVectorKsi[j] = double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						}
						_modelCreator.KnotValueVectorsKsiDictionary.Add(patchID, KnotValueVectorKsi);
						break;

					case Attributes.knotvaluevectorheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Heta of a patch must be defined after the patchID");
						if (_modelCreator.DegreeHetaDictionary[patchID] == 0 || _modelCreator.NumberOfControlPointsHetaDictionary[patchID] == 0)
						{
							throw new ArgumentOutOfRangeException("Degree Heta and number of Control Points per axis Heta must be defined before Knot Value Vector Heta.");
						}
						numberOfValues = _modelCreator.NumberOfControlPointsHetaDictionary[patchID] + _modelCreator.DegreeHetaDictionary[patchID] + 1;
						double[] KnotValueVectorHeta = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
						{
							KnotValueVectorHeta[j] = double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						}
						_modelCreator.KnotValueVectorsHetaDictionary.Add(patchID, KnotValueVectorHeta);
						break;

					case Attributes.knotvaluevectorzeta:
						if (_modelCreator.NumberOfDimensions == 2)
							throw new ArgumentOutOfRangeException("You must not define Knot Value Vector Zeta in case of 2D");
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Zeta of a patch must be defined after the patchID");
						if (_modelCreator.DegreeZetaDictionary[patchID] == 0 || _modelCreator.NumberOfControlPointsZetaDictionary[patchID] == 0)
						{
							throw new ArgumentOutOfRangeException("Degree Zeta and number of Control Points per axis Zeta must be defined before Knot Value Vector Zeta.");
						}
						numberOfValues = _modelCreator.NumberOfControlPointsZetaDictionary[patchID] + _modelCreator.DegreeZetaDictionary[patchID] + 1;
						double[] KnotValueVectorZeta = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
						{
							KnotValueVectorZeta[j] = double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						}
						_modelCreator.KnotValueVectorsZetaDictionary.Add(patchID, KnotValueVectorZeta);
						break;

					case Attributes.patchcpid:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Control Points ID of a patch must be defined after the patchID");
						int numberOfPatchCP = (_modelCreator.NumberOfDimensions == 3) ? _modelCreator.NumberOfControlPointsKsiDictionary[patchID] *
							_modelCreator.NumberOfControlPointsHetaDictionary[patchID] *
							_modelCreator.NumberOfControlPointsZetaDictionary[patchID] :
							_modelCreator.NumberOfControlPointsKsiDictionary[patchID] *
							_modelCreator.NumberOfControlPointsHetaDictionary[patchID];
						localControlPointIDs = new int[numberOfPatchCP];
						for (int j = 0; j < numberOfPatchCP; j++)
						{
							localControlPointIDs[j] = int.Parse(line[j + 1]);
						}
						_modelCreator.ControlPointIDsDictionary.Add(patchID, localControlPointIDs);
						break;

					case Attributes.cpcoord:
						_modelCreator.NumberOfControlPoints = Int32.Parse(line[1]);
						for (int j = 0; j < _modelCreator.NumberOfControlPoints; j++)
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
							_modelCreator.ControlPointsDictionary.Add(controlPointGlobalID, controlPoint);
						}
						break;

					case Attributes.end:
						_modelCreator.CreateModelData();
						return;
				}
			}
		}
	}
}
