using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    /// <summary>
	/// Patch entity of Isogeometric analysis that is similar to FEM Subdomain.
	/// </summary>
	public class Patch : ISubdomain
	{
		private readonly List<ControlPoint> controlPoints = new List<ControlPoint>();

		/// <summary>
		/// Boolean that implements equivalent property of <see cref="ISubdomain"/>.
		/// </summary>
		public bool ConnectivityModified { get; set; } = true;

		/// <summary>
		/// Ordering of constrained Control Points.
		/// </summary>
		public ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }

		/// <summary>
		/// Table containing constrained degrees of freedom and their value in a patch.
		/// </summary>
		public Table<INode, IDofType, double> Constraints { get; } = new Table<INode, IDofType, double>();

		/// <summary>
		/// Table containing control point loads per dof.
		/// </summary>
		public Table<ControlPoint, IDofType, double> ControlPointLoads { get; set; }

		/// <summary>
		/// Returns a list of the Control Points of the Patch.
		/// </summary>
		public IReadOnlyList<ControlPoint> ControlPoints => controlPoints;

		/// <summary>
		/// Degree of patch parametric axis Heta.
		/// </summary>
		public int DegreeHeta { get; set; }

		/// <summary>
		/// Degree of patch parametric axis Ksi.
		/// </summary>
		public int DegreeKsi { get; set; }

		/// <summary>
		/// Degree of patch parametric axis Zeta.
		/// </summary>
		public int DegreeZeta { get; set; }

		/// <summary>
		/// Dictionary containing the edges of a patch.
		/// </summary>
		public Dictionary<int, Edge> EdgesDictionary { get; } = new Dictionary<int, Edge>();

		/// <summary>
		/// Return a list with the elements of the patch.
		/// </summary>
		IReadOnlyList<IElement> ISubdomain.Elements => Elements;

		/// <summary>
		/// A list containing the elements of a patch.
		/// </summary>
		public List<Element> Elements { get; } = new List<Element>();

		/// <summary>
		/// Dictionary containing the faces of the patch.
		/// </summary>
		public Dictionary<int, Face> FacesDictionary { get; } = new Dictionary<int, Face>();

		/// <summary>
		/// Vector of the patch forces.
		/// </summary>
		public Vector Forces { get; set; }

		/// <summary>
		/// Ordering of patch free degrees of freedom.
		/// </summary>
		public ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

		/// <summary>
		/// Patch ID.
		/// </summary>
		public int ID { get; }

		/// <summary>
		/// Knot Value Vector of patch parametric axis Heta.
		/// </summary>
		public Vector KnotValueVectorHeta { get; set; }

		/// <summary>
		/// Knot Value Vector of patch parametric axis Ksi.
		/// </summary>
		public Vector KnotValueVectorKsi { get; set; }

		/// <summary>
		/// Knot Value Vector of patch parametric axis Zeta.
		/// </summary>
		public Vector KnotValueVectorZeta { get; set; }

		/// <summary>
		/// Patch material.
		/// </summary>
		public IFiniteElementMaterial Material { get; set; }

		/// <summary>
		/// Return a list of the Controls Points of the Patchas <see cref="INode"/>.
		/// </summary>
		IReadOnlyList<INode> ISubdomain.Nodes => controlPoints;

		/// <summary>
		/// Number of Control Points per patch parametrix axis Heta.
		/// </summary>
		public int NumberOfControlPointsHeta { get; set; }

		/// <summary>
		/// Number of Control Points per patch parametrix axis Ksi.
		/// </summary>
		public int NumberOfControlPointsKsi { get; set; }

		/// <summary>
		/// Number of Control Points per patch parametrix axis Zeta.
		/// </summary>
		public int NumberOfControlPointsZeta { get; set; }

		/// <summary>
		/// Dimensionality of the problem.
		/// </summary>
		public int NumberOfDimensions { get; set; }

		/// <summary>
		/// Boolean that implements equivalent property of <see cref="ISubdomain"/>.
		/// </summary>
		public bool StiffnessModified { get; set; }

		/// <summary>
		/// Thickness of 2D patch elements.
		/// </summary>
		public double Thickness { get; set; }

		/// <summary>
		/// Implements equivalent method of <see cref="ISubdomain"/>.
		/// </summary>
		public double[] CalculateElementDisplacements(IElement element, IVectorView globalDisplacementVector)
		{
			var elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
			SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
			return elementNodalDisplacements;
		}

		/// <summary>
		/// Implements equivalent method of <see cref="ISubdomain"/>.
		/// </summary>
		public double[] CalculateElementIncrementalConstraintDisplacements(IElement element, double constraintScalingFactor)
		{
			var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
			SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
			return elementNodalDisplacements;
		}

		/// <summary>
		/// Implements equivalent method of <see cref="ISubdomain"/>.
		/// </summary>
		public void ClearMaterialStresses()
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// Creates the patch data.
		/// </summary>
		public void CreatePatchData()
		{
			if (this.NumberOfDimensions == 2)
			{
				CreatePatchData2D();
			}
			else
			{
				CreatePatchData3D();
			}
		}

		/// <summary>
		/// Implements equivalent method of <see cref="ISubdomain"/>.
		/// </summary>
		public void ExtractConstraintsFromGlobal(Table<INode, IDofType, double> globalConstraints)
		{
			foreach (ControlPoint controlPoint in ControlPoints)
			{
				bool isControlPointConstrained = globalConstraints.TryGetDataOfRow(controlPoint,
					out IReadOnlyDictionary<IDofType, double> constraintsOfNode);
				if (isControlPointConstrained)
				{
					foreach (var dofDisplacementPair in constraintsOfNode)
					{
						Constraints[controlPoint, dofDisplacementPair.Key] = dofDisplacementPair.Value;
					}
				}
			}
		}

		/// <summary>
		/// Implements equivalent method of <see cref="ISubdomain"/>.
		/// </summary>
		public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
		{
			var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs);
			foreach (Element element in Elements)
			{
				double[] localSolution = CalculateElementDisplacements(element, solution);
				double[] localdSolution = CalculateElementDisplacements(element, dSolution);
				element.ElementType.CalculateStresses(element, localSolution, localdSolution);
				if (element.ElementType.MaterialModified)
					element.Patch.StiffnessModified = true;
				var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
				FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
			}

			return forces;
		}

		/// <summary>
		/// Implements equivalent method of <see cref="ISubdomain"/>.
		/// </summary>
		public void ResetMaterialsModifiedProperty()
		{
			this.StiffnessModified = false;
			foreach (Element element in Elements) element.ElementType.ResetMaterialModified();
		}

		/// <summary>
		/// Implements equivalent method of <see cref="ISubdomain"/>.
		/// </summary>
		public void SaveMaterialState()
		{
			foreach (Element element in Elements) element.ElementType.SaveMaterialState();
		}

		/// <summary>
		/// Implements equivalent method of <see cref="ISubdomain"/>.
		/// </summary>
		public void ScaleConstraints(double scalingFactor) => Constraints.ModifyValues((u) => scalingFactor * u);

		internal void CreateNurbsShell(GeometricalFormulation formulation)
		{
			BuildEdgesDictionary();
			CreateNURBSShells(formulation);
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

		private void BuildEdgesDictionary()
		{
			if (this.NumberOfDimensions == 2)
			{
				CreateEdges2D();
			}
			else
			{
				CreateEdges3D();
			}
		}

		private void BuildFacesDictionary()
		{
			if (this.NumberOfDimensions <= 2)
			{
				FacesDictionary.Clear();
			}
			else
			{
				CreateFaces3D();
			}
		}

		private void CreateEdges2D()
		{
			#region EdgeRight

			Edge edgeRight = new Edge();
			edgeRight.ID = 0;
			edgeRight.Degree = this.DegreeHeta;
			edgeRight.KnotValueVector = this.KnotValueVectorHeta;
			edgeRight.NumberOfControlPoints = this.NumberOfControlPointsHeta;
			edgeRight.Patch = this;
			int counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
			{
				edgeRight.ControlPointsDictionary.Add(counter++, ControlPoints[i]);
			}

			EdgesDictionary.Add(0, edgeRight);

			#endregion EdgeRight

			#region EdgeLeft

			Edge edgeLeft = new Edge();
			edgeLeft.ID = 1;
			edgeLeft.Degree = this.DegreeHeta;
			edgeLeft.KnotValueVector = this.KnotValueVectorHeta;
			edgeLeft.NumberOfControlPoints = this.NumberOfControlPointsHeta;
			edgeLeft.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
			{
				edgeLeft.ControlPointsDictionary.Add(counter++,
					ControlPoints[i + this.NumberOfControlPointsHeta * (this.NumberOfControlPointsKsi - 1)]);
			}

			EdgesDictionary.Add(1, edgeLeft);

			#endregion EdgeLeft

			#region EdgeBottom

			Edge edgeBottom = new Edge();
			edgeBottom.ID = 2;
			edgeBottom.Degree = this.DegreeKsi;
			edgeBottom.KnotValueVector = this.KnotValueVectorKsi;
			edgeBottom.NumberOfControlPoints = this.NumberOfControlPointsKsi;
			edgeBottom.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				edgeBottom.ControlPointsDictionary.Add(counter++, ControlPoints[i * this.NumberOfControlPointsHeta]);
			}

			EdgesDictionary.Add(2, edgeBottom);

			#endregion EdgeBottom

			#region EdgeUp

			Edge edgeUp = new Edge();
			edgeUp.ID = 3;
			edgeUp.Degree = this.DegreeKsi;
			edgeUp.KnotValueVector = this.KnotValueVectorKsi;
			edgeUp.NumberOfControlPoints = this.NumberOfControlPointsKsi;
			edgeUp.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				edgeUp.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsHeta + this.NumberOfControlPointsHeta - 1]);
			}

			EdgesDictionary.Add(3, edgeUp);

			#endregion EdgeUp
		}

		private void CreateEdges3D()
		{
			#region Edge1

			Edge edge1 = new Edge();
			edge1.ID = 0;
			edge1.Degree = this.DegreeZeta;
			edge1.KnotValueVector = this.KnotValueVectorZeta;
			edge1.NumberOfControlPoints = this.NumberOfControlPointsZeta;
			edge1.Patch = this;
			int counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
			{
				edge1.ControlPointsDictionary.Add(counter++, ControlPoints[i]);
			}

			EdgesDictionary.Add(0, edge1);

			#endregion Edge1

			#region Edge2

			Edge edge2 = new Edge();
			edge2.ID = 1;
			edge2.Degree = this.DegreeZeta;
			edge2.KnotValueVector = this.KnotValueVectorZeta;
			edge2.NumberOfControlPoints = this.NumberOfControlPointsZeta;
			edge2.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
			{
				edge2.ControlPointsDictionary.Add(counter++,
					ControlPoints[i + this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)]);
			}

			EdgesDictionary.Add(1, edge2);

			#endregion Edge2

			#region Edge3

			Edge edge3 = new Edge();
			edge3.ID = 2;
			edge3.Degree = this.DegreeHeta;
			edge3.KnotValueVector = this.KnotValueVectorHeta;
			edge3.NumberOfControlPoints = this.NumberOfControlPointsHeta;
			edge3.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
			{
				edge3.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsZeta]);
			}

			EdgesDictionary.Add(2, edge3);

			#endregion Edge3

			#region Edge4

			Edge edge4 = new Edge();
			edge4.ID = 3;
			edge4.Degree = this.DegreeHeta;
			edge4.KnotValueVector = this.KnotValueVectorHeta;
			edge4.NumberOfControlPoints = this.NumberOfControlPointsHeta;
			edge4.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
			{
				edge4.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsZeta + this.NumberOfControlPointsZeta - 1]);
			}

			EdgesDictionary.Add(3, edge4);

			#endregion Edge4

			#region Edge5

			Edge edge5 = new Edge();
			edge5.ID = 4;
			edge5.Degree = this.DegreeZeta;
			edge5.KnotValueVector = this.KnotValueVectorZeta;
			edge5.NumberOfControlPoints = this.NumberOfControlPointsZeta;
			int offset = this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta *
						 (this.NumberOfControlPointsKsi - 1);
			edge5.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
			{
				edge5.ControlPointsDictionary.Add(counter++, ControlPoints[i + offset]);
			}

			EdgesDictionary.Add(4, edge5);

			#endregion Edge5

			#region Edge6

			Edge edge6 = new Edge();
			edge6.ID = 5;
			edge6.Degree = this.DegreeZeta;
			edge6.KnotValueVector = this.KnotValueVectorZeta;
			edge6.NumberOfControlPoints = this.NumberOfControlPointsZeta;
			edge6.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
			{
				edge6.ControlPointsDictionary.Add(counter++,
					ControlPoints[i + this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1) + offset]);
			}

			EdgesDictionary.Add(5, edge6);

			#endregion Edge6

			#region Edge7

			Edge edge7 = new Edge();
			edge7.ID = 6;
			edge7.Degree = this.DegreeHeta;
			edge7.KnotValueVector = this.KnotValueVectorHeta;
			edge7.NumberOfControlPoints = this.NumberOfControlPointsHeta;
			edge7.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
			{
				edge7.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsZeta + offset]);
			}

			EdgesDictionary.Add(6, edge7);

			#endregion Edge7

			#region Edge8

			Edge edge8 = new Edge();
			edge8.ID = 7;
			edge8.Degree = this.DegreeHeta;
			edge8.KnotValueVector = this.KnotValueVectorHeta;
			edge8.NumberOfControlPoints = this.NumberOfControlPointsHeta;
			edge8.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
			{
				edge8.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsZeta + this.NumberOfControlPointsZeta - 1 + offset]);
			}

			EdgesDictionary.Add(7, edge8);

			#endregion Edge8

			#region Edge9

			Edge edge9 = new Edge();
			edge9.ID = 8;
			edge9.Degree = this.DegreeKsi;
			edge9.KnotValueVector = this.KnotValueVectorKsi;
			edge9.NumberOfControlPoints = this.NumberOfControlPointsKsi;
			edge9.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				edge9.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta]);
			}

			EdgesDictionary.Add(8, edge9);

			#endregion Edge9

			#region Edge10

			Edge edge10 = new Edge();
			edge10.ID = 9;
			edge10.Degree = this.DegreeKsi;
			edge10.KnotValueVector = this.KnotValueVectorKsi;
			edge10.NumberOfControlPoints = this.NumberOfControlPointsKsi;
			edge10.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				edge10.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta
								  + this.NumberOfControlPointsZeta - 1]);
			}

			EdgesDictionary.Add(9, edge10);

			#endregion Edge10

			#region Edge11

			Edge edge11 = new Edge();
			edge11.ID = 10;
			edge11.Degree = this.DegreeKsi;
			edge11.KnotValueVector = this.KnotValueVectorKsi;
			edge11.NumberOfControlPoints = this.NumberOfControlPointsKsi;
			edge11.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				edge11.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta +
								  this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)]);
			}

			EdgesDictionary.Add(10, edge11);

			#endregion Edge11

			#region Edge12

			Edge edge12 = new Edge();
			edge12.ID = 11;
			edge12.Degree = this.DegreeKsi;
			edge12.KnotValueVector = this.KnotValueVectorKsi;
			edge12.NumberOfControlPoints = this.NumberOfControlPointsKsi;
			edge12.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				edge12.ControlPointsDictionary.Add(counter++,
					ControlPoints[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta +
								  this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1) +
								  this.NumberOfControlPointsZeta - 1]);
			}

			EdgesDictionary.Add(11, edge12);

			#endregion Edge12
		}

		private void CreateElements2D(Vector singleKnotValuesKsi, Vector singleKnotValuesHeta, List<Knot> knots)
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
							elementControlPoints.Add(this.ControlPoints[controlPointID]);
						}
					}

					int elementID = i * numberOfElementsHeta + j;
					Element element = new NURBSElement2D()
					{
						ID = elementID,
						Patch = this,
						ElementType = new NURBSElement2D()
					};
					element.AddKnots(knotsOfElement);
					element.AddControlPoints(elementControlPoints.ToList());
					Elements.Add(element);
				}
			}
		}

		private void CreateElements3D(List<Knot> knots)
		{
			Vector singlesKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
			Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
			Vector singlesKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
			Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];
			Vector singlesKnotValuesZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];
			Vector multiplicityZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[1];

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
						if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
						{
							multiplicityElementKsi = (int)multiplicityKsi[i + 1] - DegreeKsi;
						}

						int multiplicityElementHeta = 0;
						if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
						{
							multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
						}

						int multiplicityElementZeta = 0;
						if (multiplicityZeta[k + 1] - this.DegreeZeta > 0)
						{
							multiplicityElementZeta = (int)multiplicityZeta[k + 1] - this.DegreeZeta;
						}

						int nurbsSupportKsi = this.DegreeKsi + 1;
						int nurbsSupportHeta = this.DegreeHeta + 1;
						int nurbsSupportZeta = this.DegreeZeta + 1;

						IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

						for (int l = 0; l < nurbsSupportKsi; l++)
						{
							for (int m = 0; m < nurbsSupportHeta; m++)
							{
								for (int n = 0; n < nurbsSupportZeta; n++)
								{
									int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHeta *
														 NumberOfControlPointsZeta +
														 (j + multiplicityElementHeta) * NumberOfControlPointsZeta +
														 (k + multiplicityElementZeta) +
														 l * NumberOfControlPointsHeta * NumberOfControlPointsZeta +
														 m * NumberOfControlPointsZeta + n;

									elementControlPoints.Add(ControlPoints[controlPointID]);
								}
							}
						}

						int elementID = i * numberOfElementsHeta * numberOfElementsZeta + j * numberOfElementsZeta + k;
						Element element = new NurbsElement3D()
						{
							ID = elementID,
							Patch = this,
							ElementType = new NurbsElement3D()
						};
						element.AddKnots(knotsOfElement);
						element.AddControlPoints(elementControlPoints.ToList());
						Elements.Add(element);
					}
				}
			}
		}

		private void CreateFaces3D()
		{
			#region FaceRight

			Face faceRight = new Face();
			faceRight.Degrees[0] = this.DegreeHeta;
			faceRight.Degrees[1] = this.DegreeZeta;
			faceRight.KnotValueVectors.Add(0, KnotValueVectorHeta);
			faceRight.KnotValueVectors.Add(1, KnotValueVectorZeta);
			faceRight.Patch = this;
			int counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
			{
				for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
				{
					faceRight.ControlPointsDictionary.Add(counter++,
						ControlPoints[j + this.NumberOfControlPointsZeta * i]);
				}
			}

			FacesDictionary.Add(0, faceRight);

			#endregion FaceRight

			#region FaceLeft

			Face faceLeft = new Face();
			faceLeft.Degrees[0] = this.DegreeHeta;
			faceLeft.Degrees[1] = this.DegreeZeta;
			faceLeft.KnotValueVectors.Add(0, KnotValueVectorHeta);
			faceLeft.KnotValueVectors.Add(1, KnotValueVectorZeta);
			faceLeft.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
			{
				for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
				{
					faceLeft.ControlPointsDictionary.Add(counter++,
						ControlPoints[j + this.NumberOfControlPointsZeta * i +
									  this.NumberOfControlPointsHeta * this.NumberOfControlPointsZeta *
									  (this.NumberOfControlPointsKsi - 1)]);
				}
			}

			FacesDictionary.Add(1, faceLeft);

			#endregion FaceLeft

			#region FaceBottom

			Face faceBottom = new Face();
			faceBottom.Degrees[0] = this.DegreeKsi;
			faceBottom.Degrees[1] = this.DegreeHeta;
			faceBottom.KnotValueVectors.Add(0, KnotValueVectorKsi);
			faceBottom.KnotValueVectors.Add(1, KnotValueVectorHeta);
			faceBottom.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				for (int j = 0; j < this.NumberOfControlPointsHeta; j++)
				{
					faceBottom.ControlPointsDictionary.Add(counter++,
						ControlPoints[j * this.NumberOfControlPointsZeta +
									  i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta]);
				}
			}

			FacesDictionary.Add(2, faceBottom);

			#endregion FaceBottom

			#region FaceUp

			Face faceUp = new Face();
			faceUp.Degrees[0] = this.DegreeKsi;
			faceUp.Degrees[1] = this.DegreeHeta;
			faceUp.KnotValueVectors.Add(0, KnotValueVectorKsi);
			faceUp.KnotValueVectors.Add(1, KnotValueVectorHeta);
			faceUp.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				for (int j = 0; j < this.NumberOfControlPointsHeta; j++)
				{
					faceUp.ControlPointsDictionary.Add(counter++,
						ControlPoints[this.NumberOfControlPointsZeta - 1 + j * this.NumberOfControlPointsZeta +
									  i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta]);
				}
			}

			FacesDictionary.Add(3, faceUp);

			#endregion FaceUp

			#region FaceFront

			Face faceFront = new Face();
			faceFront.Degrees[0] = this.DegreeKsi;
			faceFront.Degrees[1] = this.DegreeZeta;
			faceFront.KnotValueVectors.Add(0, KnotValueVectorKsi);
			faceFront.KnotValueVectors.Add(1, KnotValueVectorZeta);
			faceFront.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
				{
					faceFront.ControlPointsDictionary.Add(counter++,
						ControlPoints[j + i * this.NumberOfControlPointsHeta * this.NumberOfControlPointsZeta]);
				}
			}

			FacesDictionary.Add(4, faceFront);

			#endregion FaceFront

			#region FaceBack

			Face faceBack = new Face();
			faceBack.Degrees[0] = this.DegreeKsi;
			faceBack.Degrees[1] = this.DegreeZeta;
			faceBack.KnotValueVectors.Add(0, KnotValueVectorKsi);
			faceBack.KnotValueVectors.Add(1, KnotValueVectorZeta);
			faceBack.Patch = this;
			counter = 0;
			for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
			{
				for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
				{
					faceBack.ControlPointsDictionary.Add(counter++,
						ControlPoints[j + i * this.NumberOfControlPointsHeta * this.NumberOfControlPointsZeta +
									  this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)]);
				}
			}

			FacesDictionary.Add(5, faceBack);

			#endregion FaceBack
		}

		private void CreateNURBSElements2D()
		{
			#region Knots

			Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
			Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];

			var knots = CreateKnots2D(singleKnotValuesKsi, singleKnotValuesHeta);

			#endregion Knots

			#region Elements

			CreateElements2D(singleKnotValuesKsi, singleKnotValuesHeta, knots);

			#endregion Elements
		}

		private void CreateNURBSElements3D()
		{
			#region Knots

			Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
			Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
			Vector singleKnotValuesZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];

			var knots = CreateKnots3D(singleKnotValuesKsi, singleKnotValuesHeta, singleKnotValuesZeta);

			#endregion Knots

			#region Elements

			CreateElements3D(knots);

			#endregion Elements
		}

		private void CreateNURBSShells(GeometricalFormulation formulation)
		{
			#region Knots

			Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
			Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];

			var knots = CreateShellKnots(singleKnotValuesKsi, singleKnotValuesHeta);

			#endregion Knots

			#region Elements

			CreateShellElements(singleKnotValuesKsi, singleKnotValuesHeta, knots, formulation);

			#endregion Elements
		}

		private void CreatePatchData2D()
		{
			CreateNURBSElements2D();
			BuildEdgesDictionary();
		}

		private void CreatePatchData3D()
		{
			CreateNURBSElements3D();
			BuildEdgesDictionary();
			BuildFacesDictionary();
		}

		private void CreateShellElements(Vector singleKnotValuesKsi, Vector singleKnotValuesHeta, List<Knot> knots, GeometricalFormulation formulation)
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
							elementControlPoints.Add(ControlPoints[controlPointID]);
						}
					}

					int elementID = i * numberOfElementsHeta + j;
					Element element = null;
					switch (formulation)
					{
						case GeometricalFormulation.Linear:
							element = new NurbsKirchhoffLoveShellElement()
							{
								ID = elementID,
								Patch = this,
								ElementType = new NurbsKirchhoffLoveShellElement()
							};
							element.AddKnots(knotsOfElement);
							element.AddControlPoints(elementControlPoints);

							break;
						case GeometricalFormulation.NonLinear:
							element = new NurbsKirchhoffLoveShellElementNL(new ShellElasticMaterial2Dtransformationb()
							{
								YoungModulus = this.Material.YoungModulus,
								PoissonRatio = this.Material.PoissonRatio
							}, knotsOfElement, elementControlPoints, this, this.Thickness)
							{
								ID = elementID,
								Patch = this,
								Thickness = this.Thickness,
								ElementType = new NurbsKirchhoffLoveShellElementNL(new ShellElasticMaterial2Dtransformationb()
								{
									YoungModulus = this.Material.YoungModulus,
									PoissonRatio = this.Material.PoissonRatio
								}, knotsOfElement, elementControlPoints.ToList(), this, this.Thickness)
								{
									ID = elementID,
									Patch = this,
									Thickness = this.Thickness
								}
							};
							break;
						case GeometricalFormulation.SectionNonLinear:
							element = new NurbsKirchhoffLoveShellElementSectionNL(
								new ShellElasticSectionMaterial2D()
								{
									YoungModulus = this.Material.YoungModulus,
									PoissonRatio = this.Material.PoissonRatio,
									Thickness = this.Thickness
								}, knotsOfElement, elementControlPoints, this, this.Thickness)
							{
								ID = elementID,
								Patch = this,
								Thickness = this.Thickness,
								ElementType = new NurbsKirchhoffLoveShellElementSectionNL(new ShellElasticSectionMaterial2D()
								{
									YoungModulus = this.Material.YoungModulus,
									PoissonRatio = this.Material.PoissonRatio,
									Thickness = this.Thickness
								}, knotsOfElement, elementControlPoints, this, this.Thickness)
								{
									ID = elementID,
									Patch = this,
									Thickness = this.Thickness
								}
							};
								
							break;
					}
					Elements.Add(element);
				}
			}
		}

		private void DefineControlPointsFromElements()
		{
			var cpComparer = Comparer<ControlPoint>.Create((node1, node2) => node1.ID - node2.ID);
			var cpSet = new SortedSet<ControlPoint>(cpComparer);
			foreach (Element element in Elements)
			{
				foreach (ControlPoint node in element.ControlPoints) cpSet.Add(node);
			}
			controlPoints.AddRange(cpSet);
		}

		public IVector GetRHSFromSolutionWithInitialDisplacementsEffect(IVectorView solution, IVectorView dSolution, Dictionary<int, INode> boundaryNodes, Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements, int nIncrement, int totalIncrements)
		{
			throw new NotImplementedException();
		}
	}
}
