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
		/// Returns a list of the Control Points of the Patch.
		/// </summary>
		public List<ControlPoint> ControlPoints => controlPoints;

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
		/// Return a list of the Controls Points of the Patchas <see cref="INode"/>.
		/// </summary>
		IReadOnlyList<INode> ISubdomain.Nodes => controlPoints;
        
		/// <summary>
		/// Dimensionality of the problem.
		/// </summary>
		public int NumberOfDimensions { get; set; }

		/// <summary>
		/// Boolean that implements equivalent property of <see cref="ISubdomain"/>.
		/// </summary>
		public bool StiffnessModified { get; set; }


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
