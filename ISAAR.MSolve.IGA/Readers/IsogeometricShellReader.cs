namespace MGroup.IGA.Readers
{
	using System;
	using System.Collections.Generic;
	using System.Globalization;

	using MGroup.IGA.Entities;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Materials;

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
		private readonly Model _model;
		private readonly string _filename;

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
		public IsogeometricShellReader(Model model, string filename)
		{
			_model = model;
			_filename = filename;
		}

		private Dictionary<int, int[]> ControlPointIDsDictionary = new Dictionary<int, int[]>();

		/// <summary>
		/// Creates model from custom isogeometric shell file.
		/// </summary>
		public void CreateShellModelFromFile(GeometricalFormulation formulation)
		{
			char[] delimeters = { ' ', '=', '\t' };
			Attributes? name = null;

			int patchID = -1;
			int numberOfValues = 0;
			int[] localControlPointIDs;
			int counterElementID = 0;
			int counterCPID;

			string[] text = System.IO.File.ReadAllLines(_filename);

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
						_model.PatchesDictionary[patchID].NumberOfDimensions = int.Parse(line[1]);
						break;

					case Attributes.thickness:
						_model.PatchesDictionary[patchID].Thickness = double.Parse(line[1], CultureInfo.InvariantCulture);
						break;

					case Attributes.numberofpatches:
						break;

					case Attributes.material:
						_model.PatchesDictionary[patchID].Material = new ElasticMaterial2D(StressState2D.PlaneStrain) { YoungModulus = double.Parse(line[2], CultureInfo.InvariantCulture), PoissonRatio = double.Parse(line[3], CultureInfo.InvariantCulture) };
						break;

					case Attributes.patchid:
						patchID = int.Parse(line[1]);
						_model.PatchesDictionary.Add(patchID, new Patch());
						break;

					case Attributes.degreeksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Ksi of a patch must be defined after the patchID");
						_model.PatchesDictionary[patchID].DegreeKsi = int.Parse(line[1]);
						break;

					case Attributes.degreeheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Degree Heta of a patch must be defined after the patchID");
						_model.PatchesDictionary[patchID].DegreeHeta = int.Parse(line[1]);
						break;

					case Attributes.numberofcpksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Ksi of a patch must be defined after the patchID");
						_model.PatchesDictionary[patchID].NumberOfControlPointsKsi = int.Parse(line[1]);
						break;

					case Attributes.numberofcpheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Number of Control Points Heta of a patch must be defined after the patchID");
						_model.PatchesDictionary[patchID].NumberOfControlPointsHeta = int.Parse(line[1]);
						break;

					case Attributes.knotvaluevectorksi:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Ksi of a patch must be defined after the patchID");
						if (_model.PatchesDictionary[patchID].DegreeKsi == 0 || _model.PatchesDictionary[patchID].NumberOfControlPointsKsi == 0)
							throw new ArgumentOutOfRangeException("Degree Ksi and number of Control Points per axis Ksi must be defined before Knot Value Vector Ksi.");
						numberOfValues = _model.PatchesDictionary[patchID].DegreeKsi + _model.PatchesDictionary[patchID].NumberOfControlPointsKsi + 1;
						double[] KnotValueVectorKsi = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
							KnotValueVectorKsi[j] = double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						_model.PatchesDictionary[patchID].KnotValueVectorKsi = Vector.CreateFromArray(KnotValueVectorKsi);
						break;

					case Attributes.knotvaluevectorheta:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("KnotValue Vector Heta of a patch must be defined after the patchID");
						if (_model.PatchesDictionary[patchID].DegreeHeta == 0 || _model.PatchesDictionary[patchID].NumberOfControlPointsHeta == 0)
							throw new ArgumentOutOfRangeException("Degree Heta and number of Control Points per axis Heta must be defined before Knot Value Vector Heta.");
						numberOfValues = _model.PatchesDictionary[patchID].DegreeHeta + _model.PatchesDictionary[patchID].NumberOfControlPointsHeta + 1;
						double[] KnotValueVectorHeta = new double[numberOfValues];
						for (int j = 0; j < numberOfValues; j++)
							KnotValueVectorHeta[j] = double.Parse(line[j + 1], CultureInfo.InvariantCulture);
						_model.PatchesDictionary[patchID].KnotValueVectorHeta = Vector.CreateFromArray(KnotValueVectorHeta);
						break;

					case Attributes.patchcpid:
						if (patchID == -1)
							throw new ArgumentOutOfRangeException("Control Points ID of a patch must be defined after the patchID");
						int numberOfPatchCP = _model.PatchesDictionary[patchID].NumberOfControlPointsKsi *
											  _model.PatchesDictionary[patchID].NumberOfControlPointsHeta;
						localControlPointIDs = new int[numberOfPatchCP];
						for (int j = 0; j < numberOfPatchCP; j++)
						{
							localControlPointIDs[j] = int.Parse(line[j + 1]);
						}
						ControlPointIDsDictionary.Add(patchID, localControlPointIDs);
						break;

					case Attributes.cpcoord:
						var numberOfControlPoints = int.Parse(line[1]);
						for (int j = 0; j < numberOfControlPoints; j++)
						{
							i++;
							line = text[i].Split(delimeters);
							int controlPointGlobalID = int.Parse(line[0]);
							double x = double.Parse(line[1], CultureInfo.InvariantCulture);
							double y = double.Parse(line[2], CultureInfo.InvariantCulture);
							double z = double.Parse(line[3], CultureInfo.InvariantCulture);
							double w = double.Parse(line[4], CultureInfo.InvariantCulture);
							ControlPoint controlPoint = new ControlPoint()
							{ ID = controlPointGlobalID, X = x, Y = y, Z = z, WeightFactor = w };
							_model.ControlPointsDictionary.Add(controlPointGlobalID, controlPoint);
						}
						break;

					case Attributes.end:
						for (int j = 0; j < ControlPointIDsDictionary[patchID].Length; j++)
							((List<ControlPoint>)_model.PatchesDictionary[patchID].ControlPoints).Add(_model.ControlPointsDictionary[ControlPointIDsDictionary[patchID][j]]);

						_model.PatchesDictionary[patchID].CreateNurbsShell(formulation);
						foreach (var element in _model.PatchesDictionary[patchID].Elements)
							_model.ElementsDictionary.Add(counterElementID++, element);
						return;
				}
			}
		}
	}
}
