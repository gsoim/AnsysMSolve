using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    /// <summary>
	/// Defines a Control Point. Implements <see cref="IWeightedPoint"/>.
	/// </summary>
	public class ControlPoint:INode
	{
		/// <summary>
		/// List containing degree of freedom constraints.
		/// </summary>
		public List<Constraint> Constraints { get; } = new List<Constraint>();

		/// <summary>
		/// Dictionary containing the Elements adjacent to the <see cref="ControlPoint"/>.
		/// </summary>
		public Dictionary<int, Element> ElementsDictionary { get; } = new Dictionary<int, Element>();
        Dictionary<int, IElement> INode.ElementsDictionary => throw new NotImplementedException();
        /// <summary>
        /// Parametric coordinate Heta of the <see cref="ControlPoint"/>.
        /// </summary>
        public double Heta { get; set; }

		/// <summary>
		/// ID of the <see cref="ControlPoint"/>.
		/// </summary>
		public int ID { get; set; }

		/// <summary>
		/// Parametric coordinate Ksi of the <see cref="ControlPoint"/>.
		/// </summary>
		public double Ksi { get; set; }

		/// <summary>
		/// Dictionary that contains the patches the <see cref="ControlPoint"/> belongs to.
		/// </summary>
		public Dictionary<int, Patch> PatchesDictionary { get; } = new Dictionary<int, Patch>();

		/// <summary>
		/// Dictionary that contains the patches the <see cref="ControlPoint"/> belongs to.
		/// </summary>
		public Dictionary<int, ISubdomain> SubdomainsDictionary => throw new NotImplementedException();

		/// <summary>
		/// Weight factor of the <see cref="ControlPoint"/>.
		/// </summary>
		public double WeightFactor { get; set; }

		/// <summary>
		/// Cartesian coordinate X of the <see cref="ControlPoint"/>.
		/// </summary>
		public double X { get; set; }

		/// <summary>
		/// Cartesian coordinate Y of the <see cref="ControlPoint"/>.
		/// </summary>
		public double Y { get; set; }

		/// <summary>
		/// Cartesian coordinate Z of the <see cref="ControlPoint"/>.
		/// </summary>
		public double Z { get; set; }

		/// <summary>
		/// Parametric coordinate Zeta of the <see cref="ControlPoint"/>.
		/// </summary>
		public double Zeta { get; set; }

		/// <summary>
		/// Find the patches that the <see cref="ControlPoint"/> belongs.
		/// </summary>
		public void BuildPatchesDictionary()
		{
			foreach (var element in ElementsDictionary.Values)
			{
				if (!PatchesDictionary.ContainsKey(element.Patch.ID))
					PatchesDictionary.Add(element.Patch.ID, element.Patch);
			}
		}

		/// <summary>
		/// Clones the <see cref="ControlPoint"/> object.
		/// </summary>
		/// <returns>A copy of the <see cref="ControlPoint"/></returns>
		public ControlPoint Clone() => new ControlPoint()
		{
			ID = this.ID,
			X = X,
			Y = Y,
			Z = Z,
			Ksi = Ksi,
			Heta = Heta,
			Zeta = Zeta,
			WeightFactor = WeightFactor,
		};

		/// <summary>
		/// Compares <see cref="ControlPoint"/>s based on their IDs.
		/// </summary>
		/// <param name="other">A <see cref="ControlPoint"/> used for comparison with the current one.</param>
		/// <returns>An <see cref="int"/> denoting the result of the comparison.</returns>
		public int CompareTo(INode other) => this.ID - other.ID;

		/// <summary>
		/// Converts the <see cref="ControlPoint"/> to string.
		/// </summary>
		/// <returns> A <see cref="string"/> serialization of the <see cref="ControlPoint"/> object.</returns>
		public override string ToString()
		{
			var header = $"{ID}: ({X}, {Y}, {Z})";
			var constrains = new StringBuilder();
			foreach (var c in Constraints)
			{
				var con = new StringBuilder();
				if (c.DOF == StructuralDof.TranslationX) con.Append("X ,");
				if (c.DOF == StructuralDof.TranslationY) con.Append("Y ,");
				con.Append(c.DOF == StructuralDof.TranslationZ ? "Z ," : "?");
				constrains.Append(con);
			}

			var constraintsDescription = constrains.ToString();
			constraintsDescription = constraintsDescription.Length > 1
				? constraintsDescription.Substring(0, constraintsDescription.Length - 2)
				: constraintsDescription;

			return $"{header} - Con({constraintsDescription})";
		}
	}
}
