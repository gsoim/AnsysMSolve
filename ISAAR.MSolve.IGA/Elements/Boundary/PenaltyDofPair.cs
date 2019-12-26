using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using Element = ISAAR.MSolve.IGA.Entities.Element;

namespace ISAAR.MSolve.IGA.Elements.Boundary
{
    public class PenaltyDofPair : Element, IStructuralIsogeometricElement
	{
		private const double PenaltyCoefficient = 10e8;

		public PenaltyDofPair(NodalDof firstPenaltyDof, NodalDof secondPenaltyDof, double dofDifference = 0.0)
		{
			FirstPenaltyDof = firstPenaltyDof;
			SecondPenaltyDof = secondPenaltyDof;
			DofDifference = dofDifference;
			this.ControlPointsDictionary.Add(FirstPenaltyDof.Node.ID, FirstPenaltyDof.Node);
			this.ControlPointsDictionary.Add(SecondPenaltyDof.Node.ID, SecondPenaltyDof.Node);
		}

		public NodalDof FirstPenaltyDof { get; private set; }

		public NodalDof SecondPenaltyDof { get; private set; }

		public double DofDifference { get; private set; }

		public ElementDimensions ElementDimensions => ElementDimensions.Unknown;

		public CellType CellType => CellType.Unknown;

		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		public bool MaterialModified => false;

		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			throw new NotImplementedException();
		}

		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			return new double[2];
		}

		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure)
		{
			throw new NotImplementedException();
		}

		public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure)
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements)
		{
			return new Tuple<double[], double[]>(new double[0], new double[0]);
		}

		public void ClearMaterialState()
		{
		}

		public void ClearMaterialStresses()
		{
			throw new NotImplementedException();
		}

		public IMatrix DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			var penaltyElement = (PenaltyDofPair)element;
			var dofTypes = new IDofType[2][];
			dofTypes[0] = new IDofType[] { penaltyElement.FirstPenaltyDof.DofType };
			dofTypes[1] = new IDofType[] { penaltyElement.SecondPenaltyDof.DofType };
			return dofTypes;
		}

		public IMatrix MassMatrix(IElement element) => Matrix.CreateZero(2, 2);

		public void ResetMaterialModified()
		{
			throw new NotImplementedException();
		}

		public void SaveMaterialState()
		{
			throw new NotImplementedException();
		}

		public IMatrix StiffnessMatrix(IElement element) => Matrix2by2.CreateFromArray(new double[2, 2]
			{
				{PenaltyCoefficient,-PenaltyCoefficient },
				{-PenaltyCoefficient ,PenaltyCoefficient}
			});
	}

	public class NodalDof
	{
		public NodalDof(ControlPoint node, IDofType dofType)
		{
			Node = node;
			DofType = dofType;
		}

		public ControlPoint Node { get; private set; }

		public IDofType DofType { get; private set; }
	}
}
