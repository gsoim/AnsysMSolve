namespace MGroup.IGA.Entities
{
	using System;
	using System.Collections.Generic;

	using MGroup.IGA.Elements;
	using MGroup.IGA.Entities.Loads;
	using MGroup.IGA.Interfaces;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.MSolve.Discretization.FreedomDegrees;

	/// <summary>
	/// Surface boundary entity.
	/// </summary>
	public class Face : Boundary
	{
		/// <summary>
		/// <see cref="int"/> array with size two that contains the polynomial degrees of the two parametric directions of the <see cref="Face"/>.
		/// </summary>
		public int[] Degrees { get; } = new int[2];

		/// <summary>
		/// Dictionary that contains the Knot Value vectors for the two parametric directions of the <see cref="Face"/>.
		/// </summary>
		public Dictionary<int, Vector> KnotValueVectors { get; } = new Dictionary<int, Vector>();

		/// <summary>
		/// The patch that contains the <see cref="Face"/>.
		/// </summary>
		public Patch Patch { get; set; }

		/// <summary>
		/// List with the boundary conditions applied to the face.
		/// </summary>
		public List<IBoundaryCondition> BoundaryConditions { get; } = new List<IBoundaryCondition>();

		/// <summary>
		/// Dictionary that contains the degrees of freedom that belong to the face.
		/// </summary>
		public Dictionary<int, Dictionary<IDofType, int>> ControlPointDOFsDictionary { get; } = new Dictionary<int, Dictionary<IDofType, int>>();

		/// <summary>
		/// Dictionary that contains the <see cref="ControlPoint"/>s of the <see cref="Face"/>.
		/// </summary>
		public Dictionary<int, ControlPoint> ControlPointsDictionary { get; } = new Dictionary<int, ControlPoint>();

		/// <summary>
		/// Dictionary of the elements of the <see cref="Face"/>.
		/// </summary>
		public Dictionary<int, Element> ElementsDictionary { get; } = new Dictionary<int, Element>();

		/// <summary>
		/// List with the loading conditions applied to the face.
		/// </summary>
		public List<LoadingCondition> LoadingConditions { get; } = new List<LoadingCondition>();

		/// <summary>
		/// Calculate the loads applied to the face.
		/// </summary>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load.</returns>
		public Dictionary<int, double> CalculateLoads()
		{
			Dictionary<int, double> faceLoad = new Dictionary<int, double>();
			foreach (LoadingCondition loading in LoadingConditions)
			{
				Dictionary<int, double> load = CalculateLoadingCondition(loading);
				foreach (int dof in load.Keys)
				{
					if (faceLoad.ContainsKey(dof))
					{
						faceLoad[dof] += load[dof];
					}
					else
					{
						faceLoad.Add(dof, load[dof]);
					}
				}
			}

			return faceLoad;
		}

		private Dictionary<int, double> CalculateLoadingCondition(LoadingCondition loading)
		{
			LoadProvider provider = new LoadProvider();
			if (ElementsDictionary.Count == 0) CreateFaceElements();
			Dictionary<int, double> load = new Dictionary<int, double>();
			if (loading is NeumannBoundaryCondition)
			{
				foreach (Element element in ElementsDictionary.Values)
					foreach (int dof in provider.LoadNeumann(element, this, loading as NeumannBoundaryCondition).Keys)
					{
						if (load.ContainsKey(dof))
						{
							load[dof] += provider.LoadNeumann(element, this, loading as NeumannBoundaryCondition)[dof];
						}
						else
						{
							load.Add(dof, provider.LoadNeumann(element, this, loading as NeumannBoundaryCondition)[dof]);
						}
					}
			}
			else if (loading is PressureBoundaryCondition)
			{
				foreach (Element element in ElementsDictionary.Values)
					foreach (int dof in provider.LoadPressure(element, this, loading as PressureBoundaryCondition).Keys)
					{
						if (load.ContainsKey(dof))
						{
							load[dof] += provider.LoadPressure(element, this, loading as PressureBoundaryCondition)[dof];
						}
						else
						{
							load.Add(dof, provider.LoadPressure(element, this, loading as PressureBoundaryCondition)[dof]);
						}
					}
			}
			return load;
		}

		private void CreateFaceElements()
		{
			#region Knots

			Vector singleKnotValuesKsi = KnotValueVectors[0].RemoveDuplicatesFindMultiplicity()[0];
			Vector singleKnotValuesHeta = KnotValueVectors[1].RemoveDuplicatesFindMultiplicity()[0];

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

			#endregion Knots

			#region Elements

			Vector multiplicityKsi = KnotValueVectors[0].RemoveDuplicatesFindMultiplicity()[1];
			Vector multiplicityHeta = KnotValueVectors[1].RemoveDuplicatesFindMultiplicity()[1];

			int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
			int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
			if (numberOfElementsKsi * numberOfElementsHeta == 0)
			{
				throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
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
					if (multiplicityKsi[i + 1] - this.Degrees[0] > 0)
					{
						multiplicityElementKsi = (int)multiplicityKsi[i + 1] - this.Degrees[0];
					}

					int multiplicityElementHeta = 0;
					if (multiplicityHeta[j + 1] - this.Degrees[1] > 0)
					{
						multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.Degrees[1];
					}

					int nurbsSupportKsi = this.Degrees[0] + 1;
					int nurbsSupportHeta = this.Degrees[1] + 1;

					int NumberOfControlPointsHeta = KnotValueVectors[1].Length - Degrees[1] - 1;

					IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

					for (int k = 0; k < nurbsSupportKsi; k++)
					{
						for (int l = 0; l < nurbsSupportHeta; l++)
						{
							int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHeta +
								(j + multiplicityElementHeta) + k * NumberOfControlPointsHeta + l;

							elementControlPoints.Add(this.ControlPointsDictionary[controlPointID]);
						}
					}
					int elementID = i * numberOfElementsHeta + j;
					Element element = new NurbsElement2D()
					{
						ID = elementID,
						ElementType = new NurbsElement2D(),
						Patch = Patch,
						Model = Patch.Elements[0].Model
					};
					element.AddKnots(knotsOfElement);
					element.AddControlPoints(elementControlPoints);
					this.ElementsDictionary.Add(elementID, element);
					//this.PatchesDictionary[1].ElementsDictionary.Add(element.ID, element);
				}
			}

			#endregion Elements
		}
	}
}
