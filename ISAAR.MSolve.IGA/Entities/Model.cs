using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Boundary;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Entities
{
    /// <summary>
	/// Model that contains all data needed for Isogeometric analysis.
	/// </summary>
	public class Model : IStructuralModel
	{
		private IGlobalFreeDofOrdering _globalDofOrdering;


		private IList<PenaltyDofPair> PenaltyDofPairs { get; } = new List<PenaltyDofPair>();

		public IList<ISurfaceLoad> SurfaceLoads { get; } = new List<ISurfaceLoad>();

		/// <summary>
		/// <see cref="Table{TRow,TColumn,TValue}"/> that contains the constrained degree of freedom and the value of their constrains.
		/// </summary>
		public Table<INode, IDofType, double> Constraints { get; private set; } =
			new Table<INode, IDofType, double>();

		/// <summary>
		/// Return an <see cref="IEnumerable{ControlPoint}"/> with the Control Points of the <see cref="Model"/>.
		/// </summary>
		public IEnumerable<ControlPoint> ControlPoints => ControlPointsDictionary.Values;

		/// <summary>
		/// Dictionary containing the Control Points of the model.
		/// </summary>
		public Dictionary<int, ControlPoint> ControlPointsDictionary { get; } = new Dictionary<int, ControlPoint>();

		/// <summary>
		/// List with the Elements of the <see cref="Model"/>.
		/// </summary>
		IReadOnlyList<IElement> IStructuralModel.Elements => ElementsDictionary.Values.ToList();

		/// <summary>
		/// List with the Elements of the <see cref="Model"/>.
		/// </summary>
		public IList<Element> Elements => ElementsDictionary.Values.ToList();

		/// <summary>
		/// Dictionary with the Elements of the <see cref="Model"/>.
		/// </summary>
		public Dictionary<int, Element> ElementsDictionary { get; } = new Dictionary<int, Element>();

		/// <summary>
		/// <see cref="IGlobalFreeDofOrdering"/> of the degrees of freedom of the <see cref="Model"/>.
		/// </summary>
		public IGlobalFreeDofOrdering GlobalDofOrdering
		{
			get => _globalDofOrdering;
			set
			{
				_globalDofOrdering = value;
				foreach (var patch in Patches)
				{
					patch.FreeDofOrdering = GlobalDofOrdering.SubdomainDofOrderings[patch];
				}
			}
		}

		private readonly List<PenaltyDofPair> penaltyBC = new List<PenaltyDofPair>();

		public void AddPenaltyConstrainedDofPair(PenaltyDofPair penaltyDofPair)
		{
			penaltyBC.Add(penaltyDofPair);
			var id = ElementsDictionary.Keys.Last() + 1;
			penaltyDofPair.ID = id;
			Element element = new PenaltyDofPair(penaltyDofPair.FirstPenaltyDof, penaltyDofPair.SecondPenaltyDof)
			{
				ID = id,
				Patch = PatchesDictionary[0],
				ElementType = penaltyDofPair,
			};

			PatchesDictionary[0].Elements.Add(element);
			ElementsDictionary.Add(id, element);
		}

		/// <summary>
		/// List containing the loads applied to the the <see cref="Model"/>.
		/// </summary>
		public IList<Load> Loads { get; private set; } = new List<Load>();

		/// <summary>
		/// List of <see cref="IMassAccelerationHistoryLoad"/> applied to the <see cref="Model"/>.
		/// </summary>
		public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; } =
			new List<IMassAccelerationHistoryLoad>();

		/// <summary>
		/// Return an <see cref="IReadOnlyList{ControlPoint}"/> with the Control Points of the <see cref="Model"/> as <see cref="INode"/>.
		/// </summary>
		IReadOnlyList<INode> IStructuralModel.Nodes => ControlPointsDictionary.Values.ToList();

		/// <summary>
		/// Number of interfaces between patches.
		/// </summary>
		public int NumberOfInterfaces { get; set; }

		/// <summary>
		/// Number of patches of the <see cref="Model"/>.
		/// </summary>
		public int NumberOfPatches { get; set; }

		/// <summary>
		/// List with the patches of the model.
		/// </summary>
		public IList<Patch> Patches => PatchesDictionary.Values.ToList();

		/// <summary>
		/// Dictionary with the patches of the model.
		/// </summary>
		public Dictionary<int, Patch> PatchesDictionary { get; } = new Dictionary<int, Patch>();

		/// <summary>
		/// Dictionary with the patches of the model returned as <see cref="ISubdomain"/>.
		/// </summary>
		IReadOnlyList<ISubdomain> IStructuralModel.Subdomains => PatchesDictionary.Values.ToList();

		/// <summary>
		/// List of time dependent loads added to the <see cref="Model"/>.
		/// </summary>
		public IList<ITimeDependentNodalLoad> TimeDependentNodalLoads { get; private set; } =
			new List<ITimeDependentNodalLoad>();

		/// <summary>
		/// Assigns nodal loads of the <see cref="Model"/>.
		/// </summary>
		/// <param name="distributeNodalLoads"><inheritdoc cref="NodalLoadsToSubdomainsDistributor"/></param>
		public void AssignLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		{
			foreach (var patch in PatchesDictionary.Values) patch.Forces.Clear();
			AssignNodalLoads(distributeNodalLoads);
			AssignBoundaryLoads();
			//Add possible penalty forces
		}

		/// <summary>
		/// Assigns mass acceleration loads of the time step to the <see cref="Model"/>.
		/// </summary>
		/// <param name="timeStep">An <see cref="int"/> denoting the number of the time step.</param>
		public void AssignMassAccelerationHistoryLoads(int timeStep)
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// Assigns nodal loads of the <see cref="Model"/>.
		/// </summary>
		/// <param name="distributeNodalLoads"><inheritdoc cref="NodalLoadsToSubdomainsDistributor"/></param>
		public void AssignNodalLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		{
			var globalNodalLoads = new Table<INode, IDofType, double>();
			foreach (var load in Loads) globalNodalLoads.TryAdd(load.Node, load.DOF, load.Amount);

			var subdomainNodalLoads = distributeNodalLoads(globalNodalLoads);
			foreach (var idSubdomainLoads in subdomainNodalLoads)
			{
				PatchesDictionary[idSubdomainLoads.Key].Forces.AddIntoThis(idSubdomainLoads.Value);
			}
		}

		/// <summary>
		/// Assigns mass acceleration loads of the time step to the <see cref="Model"/>.
		/// </summary>
		/// <param name="timeStep">An <see cref="int"/> denoting the number of the time step.</param>
		public void AssignTimeDependentNodalLoads(int timeStep, NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		{
			var globalNodalLoads = new Table<INode, IDofType, double>();
			foreach (ITimeDependentNodalLoad load in TimeDependentNodalLoads)
			{
				globalNodalLoads.TryAdd(load.Node, load.DOF, load.GetLoadAmount(timeStep));
			}

			Dictionary<int, SparseVector> subdomainNodalLoads = distributeNodalLoads(globalNodalLoads);
			foreach (var idSubdomainLoads in subdomainNodalLoads)
			{
				PatchesDictionary[idSubdomainLoads.Key].Forces.AddIntoThis(idSubdomainLoads.Value);
			}
		}

		/// <summary>
		/// Clear the <see cref="Model"/>.
		/// </summary>
		public void Clear()
		{
			Loads.Clear();
			PatchesDictionary.Clear();
			ElementsDictionary.Clear();
			ControlPointsDictionary.Clear();
			_globalDofOrdering = null;
			Constraints.Clear();
			MassAccelerationHistoryLoads.Clear();
		}

		/// <summary>
		/// Interconnects Data Structures of the <see cref="Model"/>.
		/// </summary>
		public void ConnectDataStructures()
		{
			BuildInterconnectionData();
			AssignConstraints();
			RemoveInactiveNodalLoads();
		}

		private static void AssignEdgeLoads(Patch patch)
		{
			foreach (var edge in patch.EdgesDictionary.Values)
			{
				var edgeLoadDictionary = edge.CalculateLoads();
				foreach (var dof in edgeLoadDictionary.Keys)
				{
					if (dof != -1)
					{
						patch.Forces[dof] += edgeLoadDictionary[dof];
					}
				}
			}
		}

		private static void AssignFaceLoads(Patch patch)
		{
			foreach (Face face in patch.FacesDictionary.Values)
			{
				Dictionary<int, double> faceLoadDictionary = face.CalculateLoads();
				foreach (int dof in faceLoadDictionary.Keys)
					if (dof != -1)
						patch.Forces[dof] += faceLoadDictionary[dof];
			}
		}

		private void AssignBoundaryLoads()
		{
			foreach (var patch in PatchesDictionary.Values)
			{
				AssignEdgeLoads(patch);

				AssignFaceLoads(patch);

				AssignSurfaceLoads(patch);
			}

		}

		private void AssignSurfaceLoads(Patch patch)
		{
			foreach (var surfaceLoad in SurfaceLoads)
			{
				foreach (var element in ElementsDictionary.Values.Where(e => e is ISurfaceLoadedElement))
				{
					var dofs = element.ElementType.GetElementDofTypes(element);
					var loadedDofs = new Dictionary<int, double>();
					switch (surfaceLoad)
					{
						case SurfacePressureLoad pressure:
							loadedDofs =
								(element as ISurfaceLoadedElement).CalculateSurfacePressure(element, pressure.Pressure);
							break;
						case SurfaceDistributedLoad distributedLoad:
							loadedDofs =
								(element as ISurfaceLoadedElement).CalculateSurfaceDistributedLoad(element, distributedLoad.Dof,
									distributedLoad.Magnitude);
							break;
					}

					foreach (var dof in loadedDofs.Keys)
					{
						if (dof != -1)
							patch.Forces[dof] += loadedDofs[dof];
					}
				}
			}
		}

		private void AssignConstraints()
		{
			foreach (ControlPoint controlPoint in ControlPointsDictionary.Values)
			{
				if (controlPoint.Constraints == null) continue;
				foreach (Constraint constraint in controlPoint.Constraints)
					Constraints[controlPoint, constraint.DOF] = constraint.Amount;
			}

			foreach (Patch patch in PatchesDictionary.Values) patch.ExtractConstraintsFromGlobal(Constraints);
		}

		private void BuildElementDictionaryOfEachControlPoint()
		{
			foreach (var element in ElementsDictionary.Values)
			{
				foreach (var controlPoint in element.ControlPoints)
					controlPoint.ElementsDictionary.Add(element.ID, element);
			}
		}

		private void BuildInterconnectionData()
		{
			BuildPatchOfEachElement();
			BuildElementDictionaryOfEachControlPoint();
			foreach (var controlPoint in ControlPointsDictionary.Values) controlPoint.BuildPatchesDictionary();
		}

		private void BuildPatchOfEachElement()
		{
			foreach (var patch in PatchesDictionary.Values)
			{
				foreach (var element in patch.Elements)
				{
					element.Patch = patch;
					element.Model = this;
				}
			}
		}

		private void RemoveInactiveNodalLoads()
		{
			var activeLoads = new List<Load>(Loads.Count);
			activeLoads.AddRange(from load in Loads
								 let isConstrained = Constraints.Contains(load.Node, load.DOF)
								 where !isConstrained
								 select load);
			Loads = activeLoads;
		}
	}
}
