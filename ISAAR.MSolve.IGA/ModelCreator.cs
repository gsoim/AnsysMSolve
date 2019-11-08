using System;
using System.Collections.Generic;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA
{
    /// <summary>
	/// Supportive class for creating model with a file reader.
	/// </summary>
	public class ModelCreator
	{
		/// <summary>
		/// Defines a <see cref="ModelCreator"/>.
		/// </summary>
		/// <param name="model">An isogeometric <see cref="Model"/>.</param>
		public ModelCreator(Model model)
		{
			this.Model = model;
		}

		/// <summary>
		/// Dictionary ControlPoint IDs per Patch.
		/// </summary>
		public Dictionary<int, int[]> ControlPointIDsDictionary { get; } = new Dictionary<int, int[]>();

		/// <summary>
		/// Dictionary ControlPoint per Patch.
		/// </summary>
		public Dictionary<int, ControlPoint> ControlPointsDictionary { get; } = new Dictionary<int, ControlPoint>();

		/// <summary>
		/// Dictionary Degree Heta per Patch.
		/// </summary>
		public Dictionary<int, int> DegreeHetaDictionary { get; } = new Dictionary<int, int>();

		/// <summary>
		/// Dictionary Degree Ksi per Patch.
		/// </summary>
		public Dictionary<int, int> DegreeKsiDictionary { get; } = new Dictionary<int, int>();

		/// <summary>
		/// Dictionary Degree Zeta per Patch.
		/// </summary>
		public Dictionary<int, int> DegreeZetaDictionary { get; } = new Dictionary<int, int>();

		/// <summary>
		/// Dictionary Knot Value Vector Heta per Patch.
		/// </summary>
		public Dictionary<int, double[]> KnotValueVectorsHetaDictionary { get; } = new Dictionary<int, double[]>();

		/// <summary>
		/// Dictionary Knot Value Vector Ksi per Patch.
		/// </summary>
		public Dictionary<int, double[]> KnotValueVectorsKsiDictionary { get; } = new Dictionary<int, double[]>();

		/// <summary>
		/// Dictionary Knot Value Vector Zeta per Patch.
		/// </summary>
		public Dictionary<int, double[]> KnotValueVectorsZetaDictionary { get; } = new Dictionary<int, double[]>();

		/// <summary>
		/// Model material.
		/// </summary>
		public IFiniteElementMaterial Material { get; set; }

		/// <summary>
		/// Model to be constructed.
		/// </summary>
		public Model Model { get; private set; }

		/// <summary>
		/// Total number of Control Points for the Model.
		/// </summary>
		public int NumberOfControlPoints { get; set; }

		/// <summary>
		/// Dictionary number of control points per axis Heta per Patch.
		/// </summary>
		public Dictionary<int, int> NumberOfControlPointsHetaDictionary { get; } = new Dictionary<int, int>();

		/// <summary>
		/// Dictionary number of control points per axis Ksi per Patch.
		/// </summary>
		public Dictionary<int, int> NumberOfControlPointsKsiDictionary { get; } = new Dictionary<int, int>();

		/// <summary>
		/// Dictionary number of control points per axis Zeta per Patch.
		/// </summary>
		public Dictionary<int, int> NumberOfControlPointsZetaDictionary { get; } = new Dictionary<int, int>();

		/// <summary>
		/// Dimensionality of the problem.
		/// </summary>
		public int NumberOfDimensions { get; set; }

		/// <summary>
		/// Total number of Patches.
		/// </summary>
		public int NumberOfPatches { get; set; }

		/// <summary>
		/// Boolean that is used to split C0 continuity into new patches.
		/// </summary>
		public bool SplitToContinuousPatches { get; set; } = false;

		/// <summary>
		/// Thickness of 2D/ Shell elements.
		/// </summary>
		public double Thickness { get; set; }

		/// <summary>
		/// Creates model data.
		/// </summary>
		public void CreateModelData()
		{
			if (SplitToContinuousPatches)
			{
				AssignContinuousPatchesToModel();
			}
			else
			{
				AssignPatchesToModel();
			}
		}

		private void AssignContinuousPatchesToModel()
		{
			int counterPatch = 0;
			for (int patchID = 0; patchID < NumberOfPatches; patchID++)
			{
				#region FindSubPatches

				var tupleKsi = DetectSubPatches(Vector.CreateFromArray(KnotValueVectorsKsiDictionary[patchID]), DegreeKsiDictionary[patchID]);
				int subpatchesKsi = tupleKsi.Item1;
				Dictionary<int, Vector> subKnotVectorsKsi = tupleKsi.Item2;

				var tupleHeta = DetectSubPatches(Vector.CreateFromArray(KnotValueVectorsHetaDictionary[patchID]), DegreeHetaDictionary[patchID]);
				int subpatchesHeta = tupleHeta.Item1;
				Dictionary<int, Vector> subKnotVectorsHeta = tupleHeta.Item2;

				int subpatchesZeta = 1;
				Tuple<int, Dictionary<int, Vector>> tupleZeta = new Tuple<int, Dictionary<int, Vector>>(1, null);
				Dictionary<int, Vector> subKnotVectorsZeta = new Dictionary<int, Vector>();
				if (this.NumberOfDimensions == 3)
				{
					tupleZeta = DetectSubPatches(Vector.CreateFromArray(KnotValueVectorsZetaDictionary[patchID]), DegreeZetaDictionary[patchID]);
					subpatchesZeta = tupleZeta.Item1;
					subKnotVectorsZeta = tupleZeta.Item2;
				}

				#endregion FindSubPatches

				var controlPointIDs = FindControlPointsOfEachSubPatch(patchID, tupleKsi, tupleHeta, tupleZeta);

				for (int i = 0; i < subpatchesKsi; i++)
				{
					for (int j = 0; j < subpatchesHeta; j++)
					{
						for (int k = 0; k < subpatchesZeta; k++)
						{
							Patch patch = new Patch()
							{
								NumberOfDimensions = this.NumberOfDimensions,
								DegreeKsi = DegreeKsiDictionary[patchID],
								DegreeHeta = DegreeHetaDictionary[patchID],
								DegreeZeta = (NumberOfDimensions == 2) ? 0 : DegreeZetaDictionary[patchID],
								NumberOfControlPointsKsi = subKnotVectorsKsi[i].Length - DegreeKsiDictionary[patchID] - 1,
								NumberOfControlPointsHeta = subKnotVectorsHeta[j].Length - DegreeHetaDictionary[patchID] - 1,
								NumberOfControlPointsZeta = (NumberOfDimensions == 2) ? 0 : subKnotVectorsZeta[k].Length - DegreeZetaDictionary[patchID] - 1,
								Material = this.Material,
								Thickness = (NumberOfDimensions == 2) ? this.Thickness : 0,
								KnotValueVectorKsi = subKnotVectorsKsi[i],
								KnotValueVectorHeta = subKnotVectorsHeta[j],
								KnotValueVectorZeta = (NumberOfDimensions == 2) ? null : subKnotVectorsZeta[k]
							};

							for (int m = 0; m < controlPointIDs[i, j, k].Length; m++)
							{
								((List<ControlPoint>)patch.ControlPoints).Add(ControlPointsDictionary[ControlPointIDsDictionary[patchID][controlPointIDs[i, j, k][m]]]);
							}
							patch.CreatePatchData();
							this.Model.PatchesDictionary.Add(counterPatch++, patch);
						}
					}
				}
			}
		}

		private void AssignPatchesToModel()
		{
			int counterElementID = 0;
			int counterCPID = 0;
			for (int patchID = 0; patchID < NumberOfPatches; patchID++)
			{
				Patch patch = new Patch()
				{
					NumberOfDimensions = this.NumberOfDimensions,
					DegreeKsi = DegreeKsiDictionary[patchID],
					DegreeHeta = DegreeHetaDictionary[patchID],
					DegreeZeta = (NumberOfDimensions == 2) ? 0 : DegreeZetaDictionary[patchID],
					NumberOfControlPointsKsi = NumberOfControlPointsKsiDictionary[patchID],
					NumberOfControlPointsHeta = NumberOfControlPointsHetaDictionary[patchID],
					NumberOfControlPointsZeta = (NumberOfDimensions == 2) ? 0 : NumberOfControlPointsZetaDictionary[patchID],
					Material = this.Material,
					Thickness = (NumberOfDimensions == 2) ? this.Thickness : 0,
					KnotValueVectorKsi = Vector.CreateFromArray(KnotValueVectorsKsiDictionary[patchID]),
					KnotValueVectorHeta = Vector.CreateFromArray(KnotValueVectorsHetaDictionary[patchID]),
					KnotValueVectorZeta = (NumberOfDimensions == 2) ? null : Vector.CreateFromArray(KnotValueVectorsZetaDictionary[patchID]),
				};

				for (int j = 0; j < ControlPointIDsDictionary[patchID].Length; j++)
					((List<ControlPoint>)patch.ControlPoints).Add(ControlPointsDictionary[ControlPointIDsDictionary[patchID][j]]);
				patch.CreatePatchData();
				foreach (var element in patch.Elements)
					Model.ElementsDictionary.Add(counterElementID++, element);
				foreach (var controlPoint in patch.ControlPoints)
					Model.ControlPointsDictionary.Add(counterCPID++, controlPoint);

				this.Model.PatchesDictionary.Add(patchID, patch);
			}
		}

		private int[][] CalculateAxisControlPointsIDs(Tuple<int, Dictionary<int, Vector>> tupleAxis, Dictionary<int, int> DegreeDictionary, int patchID)
		{
			int numberOfSubpatches = tupleAxis.Item1;
			Dictionary<int, Vector> subKnotVectorsKsi = tupleAxis.Item2;
			int[][] axisControlPointIDs = new int[numberOfSubpatches][];
			int counterCPKsi = 0;
			for (int i = 0; i < numberOfSubpatches; i++)
			{
				int length = subKnotVectorsKsi[i].Length - DegreeDictionary[patchID] - 1;
				axisControlPointIDs[i] = new int[length];
				for (int j = 0; j < length; j++)
				{
					axisControlPointIDs[i][j] = counterCPKsi - i;
					counterCPKsi++;
				}
			}
			return axisControlPointIDs;
		}

		private Tuple<int, Dictionary<int, Vector>> DetectSubPatches(Vector knotValueVector, int degree)
		{
			Dictionary<int, Vector> SubKnotVectors = new Dictionary<int, Vector>();
			Vector[] result = knotValueVector.RemoveDuplicatesFindMultiplicity();
			Vector singleValues = result[0];
			Vector multiplicity = result[1];

			int initialKnotVectorPosition = 0;
			int endingKnotVectorPosition = 0;
			int counterPatch = 0;
			for (int i = 0; i < singleValues.Length - 1; i++)
			{
				if (multiplicity[i + 1] - multiplicity[i] + 1 == degree)
				{
					endingKnotVectorPosition = i + (int)multiplicity[i];
					if (initialKnotVectorPosition == 0)
					{
						int length = endingKnotVectorPosition - initialKnotVectorPosition + 1 + degree;
						Vector subKnotVector = Vector.CreateZero(length);
						for (int j = 0; j < endingKnotVectorPosition - initialKnotVectorPosition + 1; j++)
							subKnotVector[j] = knotValueVector[initialKnotVectorPosition + j];
						for (int j = endingKnotVectorPosition - initialKnotVectorPosition + 1; j < length; j++)
							subKnotVector[j] = knotValueVector[endingKnotVectorPosition];
						SubKnotVectors.Add(counterPatch++, subKnotVector);
						initialKnotVectorPosition = endingKnotVectorPosition;
					}
					else
					{
						int length = endingKnotVectorPosition - initialKnotVectorPosition + 2 + degree;
						Vector subKnotVector = Vector.CreateZero(length);
						subKnotVector[0] = knotValueVector[initialKnotVectorPosition];
						for (int j = 1; j < endingKnotVectorPosition - initialKnotVectorPosition + 2; j++)
							subKnotVector[j] = knotValueVector[initialKnotVectorPosition + j - 1];
						for (int j = endingKnotVectorPosition - initialKnotVectorPosition + 2; j < length; j++)
							subKnotVector[j] = knotValueVector[endingKnotVectorPosition];
						SubKnotVectors.Add(counterPatch++, subKnotVector);
						initialKnotVectorPosition = endingKnotVectorPosition;
					}
				}
				else if (i == singleValues.Length - 2 && multiplicity[i + 1] - multiplicity[i] + 1 != degree)
				{
					endingKnotVectorPosition = knotValueVector.Length - 1;
					if (initialKnotVectorPosition == 0)
					{
						SubKnotVectors.Add(counterPatch++, knotValueVector);
					}
					else
					{
						int length = endingKnotVectorPosition - initialKnotVectorPosition + 2;
						Vector subKnotVector = Vector.CreateZero(length);
						subKnotVector[0] = knotValueVector[initialKnotVectorPosition];
						for (int j = 1; j < endingKnotVectorPosition - initialKnotVectorPosition + 2; j++)
							subKnotVector[j] = knotValueVector[initialKnotVectorPosition + j - 1];
						SubKnotVectors.Add(counterPatch++, subKnotVector);
					}
				}
				else if (initialKnotVectorPosition == singleValues.Length - degree)
				{
					endingKnotVectorPosition = knotValueVector.Length - 1;
					int length = endingKnotVectorPosition - initialKnotVectorPosition + 2;
					Vector subKnotVector = Vector.CreateZero(length);
					subKnotVector[0] = knotValueVector[initialKnotVectorPosition];
					for (int j = 1; j < endingKnotVectorPosition - initialKnotVectorPosition + 2; j++)
						subKnotVector[j] = knotValueVector[initialKnotVectorPosition + j - 1];
					SubKnotVectors.Add(counterPatch++, subKnotVector);
					break;
				}
			}
			return new Tuple<int, Dictionary<int, Vector>>(counterPatch, SubKnotVectors);
		}

		private int[,,][] FindControlPointsOfEachSubPatch(int patchID, Tuple<int, Dictionary<int, Vector>> tupleKsi, Tuple<int, Dictionary<int, Vector>> tupleHeta, Tuple<int, Dictionary<int, Vector>> tupleZeta)
		{
			int[,,][] controlPointIDs;
			if (NumberOfDimensions == 2)
			{
				int[][] axisKsiControlPointIDs = CalculateAxisControlPointsIDs(tupleKsi, DegreeKsiDictionary, patchID);
				int[][] axisHetaControlPointIDs = CalculateAxisControlPointsIDs(tupleHeta, DegreeHetaDictionary, patchID);

				int numberOfSubpatchesKsi = axisKsiControlPointIDs.Length;
				int numberOfSubpatchesHeta = axisHetaControlPointIDs.Length;

				controlPointIDs = new int[numberOfSubpatchesKsi, numberOfSubpatchesHeta, 1][];

				for (int i = 0; i < numberOfSubpatchesKsi; i++)
				{
					for (int j = 0; j < numberOfSubpatchesHeta; j++)
					{
						int numberOfSubPatchCP = (tupleKsi.Item2[i].Length - DegreeKsiDictionary[patchID] - 1) *
							(tupleHeta.Item2[j].Length - DegreeHetaDictionary[patchID] - 1);
						controlPointIDs[i, j, 0] = new int[numberOfSubPatchCP];
						for (int m = 0; m < axisKsiControlPointIDs[i].Length; m++)
						{
							for (int n = 0; n < axisHetaControlPointIDs[j].Length; n++)
							{
								int indexLocal = m * axisHetaControlPointIDs[j].Length + n;
								int indexGlobal = axisKsiControlPointIDs[i][m] * NumberOfControlPointsHetaDictionary[patchID] +
									axisHetaControlPointIDs[j][n];
								controlPointIDs[i, j, 0][indexLocal] = indexGlobal;
							}
						}
					}
				}
			}
			else
			{
				int[][] axisKsiControlPointIDs = CalculateAxisControlPointsIDs(tupleKsi, DegreeKsiDictionary, patchID);
				int[][] axisHetaControlPointIDs = CalculateAxisControlPointsIDs(tupleHeta, DegreeHetaDictionary, patchID);
				int[][] axisZetaControlPointIDs = CalculateAxisControlPointsIDs(tupleZeta, DegreeZetaDictionary, patchID);

				int numberOfSubpatchesKsi = axisKsiControlPointIDs.Length;
				int numberOfSubpatchesHeta = axisHetaControlPointIDs.Length;
				int numberOfSubpatchesZeta = (NumberOfDimensions == 3) ? axisZetaControlPointIDs.Length : 1;

				controlPointIDs = new int[numberOfSubpatchesKsi, numberOfSubpatchesHeta, numberOfSubpatchesZeta][];

				for (int i = 0; i < numberOfSubpatchesKsi; i++)
				{
					for (int j = 0; j < numberOfSubpatchesHeta; j++)
					{
						for (int k = 0; k < numberOfSubpatchesZeta; k++)
						{
							int numberOfSubPatchCP = (tupleKsi.Item2[i].Length - DegreeKsiDictionary[patchID] - 1) *
								(tupleHeta.Item2[j].Length - DegreeHetaDictionary[patchID] - 1) *
								(tupleZeta.Item2[k].Length - DegreeZetaDictionary[patchID] - 1);
							controlPointIDs[i, j, k] = new int[numberOfSubPatchCP];
							for (int m = 0; m < axisKsiControlPointIDs[i].Length; m++)
							{
								for (int n = 0; n < axisHetaControlPointIDs[j].Length; n++)
								{
									for (int p = 0; p < axisZetaControlPointIDs[k].Length; p++)
									{
										int indexLocal = m * axisHetaControlPointIDs[j].Length * axisZetaControlPointIDs[k].Length +
											n * axisZetaControlPointIDs[k].Length + p;
										int indexGlobal = axisKsiControlPointIDs[i][m] * NumberOfControlPointsHetaDictionary[patchID] * NumberOfControlPointsZetaDictionary[patchID] +
											axisHetaControlPointIDs[j][n] * NumberOfControlPointsZetaDictionary[patchID] + axisZetaControlPointIDs[k][p];
										controlPointIDs[i, j, k][indexLocal] = indexGlobal;
									}
								}
							}
						}
					}
				}
			}

			return controlPointIDs;
		}
	}
}
