using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Continuum;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Entities
{
    public class Edge : Boundary
	{
        private readonly int _numberOfCpHeta;
        private readonly int _numberOfCpZeta;

        public Edge(int numberOfCpHeta, int NumberOfCpZeta)
        {
            _numberOfCpHeta = numberOfCpHeta;
            _numberOfCpZeta = NumberOfCpZeta;
        }
		/// <summary>
		/// List of BoundaryCondition applied to the edge.
		/// </summary>
		public List<IBoundaryCondition> BoundaryConditions { get; } = new List<IBoundaryCondition>();

		/// <summary>
		/// Dictionary that contains degrees of freedom that belong to the edge.
		/// </summary>
		public Dictionary<int, Dictionary<IDofType, int>> ControlPointDOFsDictionary { get; } = new Dictionary<int, Dictionary<IDofType, int>>();

		/// <summary>
		/// Dictionary that contains <see cref="ControlPoint"/> that belong to the edge.
		/// </summary>
		public Dictionary<int, ControlPoint> ControlPointsDictionary { get; } = new Dictionary<int, ControlPoint>();

		/// <summary>
		/// Polynomial degree used for the representation of the edge.
		/// </summary>
		public int Degree { get; set; }

		/// <summary>
		/// Dictionary containing elements of the edge.
		/// </summary>
		public Dictionary<int, Element> ElementsDictionary { get; } = new Dictionary<int, Element>();

		/// <summary>
		/// Edge ID.
		/// </summary>
		public int ID { get; set; }

		/// <summary>
		/// Knot Value Vector that describes the surface.
		/// </summary>
		public Vector KnotValueVector { get; set; }

		/// <summary>
		/// List of loading conditions applied to the surface.
		/// </summary>
		public List<LoadingCondition> LoadingConditions { get; } = new List<LoadingCondition>();

		/// <summary>
		/// Multitude of <see cref="ControlPoint"/> of the edge.
		/// </summary>
		public int NumberOfControlPoints { get; set; }

		/// <summary>
		/// The patch that the edge belong to.
		/// </summary>
		public Patch Patch { get; set; }

		/// <summary>
		/// Calculates the loads applied to the edge.
		/// </summary>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the loads.</returns>
		public Dictionary<int, double> CalculateLoads()
		{
			Dictionary<int, double> edgeLoad = new Dictionary<int, double>();
			foreach (LoadingCondition loading in LoadingConditions)
			{
				Dictionary<int, double> load = CalculateLoadingCondition(loading);
				foreach (int dof in load.Keys)
				{
					if (edgeLoad.ContainsKey(dof))
					{
						edgeLoad[dof] += load[dof];
					}
					else
					{
						edgeLoad.Add(dof, load[dof]);
					}
				}
			}

			return edgeLoad;
		}

		private Dictionary<int, double> CalculateLoadingCondition(LoadingCondition loading)
		{
			var provider = new LoadProvider();
			if (ElementsDictionary.Count == 0) CreateEdgeElements();
			var load = new Dictionary<int, double>();
			switch (loading)
			{
				case NeumannBoundaryCondition condition:
					CalculateNeumann(provider, condition, load);
					break;

				case PressureBoundaryCondition condition:
					CalculatePressure(provider, condition, load);
					break;
			}
			return load;
		}

		private void CalculateNeumann(LoadProvider provider, NeumannBoundaryCondition condition, Dictionary<int, double> load)
		{
			foreach (var element in ElementsDictionary.Values)
			{
				var loadNeumann = provider.LoadNeumann(element, this, condition);
				foreach (int dof in loadNeumann.Keys)
				{
					if (load.ContainsKey(dof))
					{
						load[dof] += loadNeumann[dof];
					}
					else
					{
						load.Add(dof, loadNeumann[dof]);
					}
				}
			}
		}

		private void CalculatePressure(LoadProvider provider, PressureBoundaryCondition condition, Dictionary<int, double> load)
		{
			foreach (Element element in ElementsDictionary.Values)
			{
				foreach (int dof in provider.LoadPressure(element, this, condition).Keys)
				{
					if (load.ContainsKey(dof))
					{
						load[dof] += provider.LoadPressure(element, this, condition)[dof];
					}
					else
					{
						load.Add(dof, provider.LoadPressure(element, this, condition)[dof]);
					}
				}
			}
		}

		private void CreateEdgeElements()
		{
			#region Knots

			Vector singleKnotValuesKsi = KnotValueVector.RemoveDuplicatesFindMultiplicity()[0];

			List<Knot> knots = new List<Knot>();

			int id = 0;
			for (int i = 0; i < singleKnotValuesKsi.Length; i++)
			{
				knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = 0.0, Zeta = 0.0 });
				id++;
			}

			#endregion Knots

			#region Elements

			Vector multiplicityKsi = KnotValueVector.RemoveDuplicatesFindMultiplicity()[1];

			int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
			if (numberOfElementsKsi == 0)
			{
				throw new ArgumentNullException("Number of Elements should be defined before Element Connectivity");
			}

			for (int i = 0; i < numberOfElementsKsi; i++)
			{
				IList<Knot> knotsOfElement = new List<Knot>
				{
					knots[i],
					knots[i + 1]
				};

				int multiplicityElementKsi = 0;
				if (multiplicityKsi[i + 1] - this.Degree > 0)
				{
					multiplicityElementKsi = (int)multiplicityKsi[i + 1] - this.Degree;
				}

				int nurbsSupportKsi = this.Degree + 1;

				IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

				for (int k = 0; k < nurbsSupportKsi; k++)
				{
					int controlPointID = i + multiplicityElementKsi + k;
					elementControlPoints.Add(this.ControlPointsDictionary[controlPointID]);
				}

				int elementID = i;
                var gaussQuadrature= new GaussQuadrature();
                var gaussPoints = gaussQuadrature.CalculateElementGaussPoints(
                    Degree, knotsOfElement).ToArray();
				var edgeCP=CalculateEdgeControlPoints(this, elementControlPoints, _numberOfCpHeta, _numberOfCpZeta);
                var nurbs = new Nurbs1D(Degree, KnotValueVector.CopyToArray(), edgeCP.ToArray(), gaussPoints);
                var element = new Element
				{
					ID = elementID,
					ElementType = new ContinuumElement1D(nurbs, gaussPoints,_numberOfCpHeta, _numberOfCpZeta),
					Patch = this.Patch,
					Model = Patch.Elements[0].Model
				};

				element.AddKnots(knotsOfElement);
				element.AddControlPoints(elementControlPoints);
				this.ElementsDictionary.Add(elementID, element);
			}

			#endregion Elements
		}


        public List<ControlPoint> CalculateEdgeControlPoints(Edge edge, IList<ControlPoint> controlPoints,
            int numberOfCPHeta, int numberOfCPZeta)
        {
            var edgeCP = new List<ControlPoint>();

            foreach (var controlPoint in controlPoints)
            {
                if (edge.Patch.NumberOfDimensions == 2)
                {
                    edgeCP.Add(new ControlPoint()
                    {
                        ID = (edge.ID < 2)
                            ? controlPoint.ID % numberOfCPHeta
							: controlPoint.ID / numberOfCPHeta,
                        Ksi = controlPoint.Ksi,
                        Heta = controlPoint.Heta,
                        Zeta = controlPoint.Zeta,
                        X = controlPoint.X,
                        Y = controlPoint.Y,
                        Z = controlPoint.Z,
                        WeightFactor = controlPoint.WeightFactor,
                    });
                }
                else
                {
                    var id = FindAxisControlPointId3D(edge, controlPoint, numberOfCPHeta, numberOfCPZeta);
                    edgeCP.Add(new ControlPoint()
                    {
                        ID = id,
                        Ksi = controlPoint.Ksi,
                        Heta = controlPoint.Heta,
                        Zeta = controlPoint.Zeta,
                        X = controlPoint.X,
                        Y = controlPoint.Y,
                        Z = controlPoint.Z,
                        WeightFactor = controlPoint.WeightFactor,
                    });
                }
            }

            return edgeCP;
        }


        private int FindAxisControlPointId3D(Edge edge, ControlPoint controlPoint,
            int numberOfCPHeta, int numberOfCPZeta)
        {
            int ID = -1;
            switch (edge.ID)
            {
                case 1:
                case 2:
                case 5:
                case 6:
                    ID = controlPoint.ID % (numberOfCPHeta * numberOfCPZeta) %
                         numberOfCPZeta;
                    break;

                case 3:
                case 4:
                case 7:
                case 8:
                    ID = controlPoint.ID % (numberOfCPHeta * numberOfCPZeta) /
                         numberOfCPZeta;
                    break;

                case 9:
                case 10:
                case 11:
                case 12:
                    ID = controlPoint.ID / (numberOfCPHeta * numberOfCPZeta);
                    break;
            }

            return ID;
        }


	}
}
