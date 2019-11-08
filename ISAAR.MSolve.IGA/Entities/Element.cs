using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    /// <summary>
	/// Generic element class that contains the basic functionality of an Isogeometric Element.
	/// </summary>
	public class Element : IElement
	{
		/// <summary>
		/// Returns an <see cref="IEnumerable{ControlPoint}"/> that belong to the support of the <see cref="Element"/>.
		/// </summary>
		public IEnumerable<ControlPoint> ControlPoints => ControlPointsDictionary.Values;

		/// <summary>
		/// <see cref="Dictionary{TKey,TValue}"/> where int is the <see cref="ControlPoint"/> ID.
		/// </summary>
		public Dictionary<int, ControlPoint> ControlPointsDictionary { get; } = new Dictionary<int, ControlPoint>();

		/// <summary>
		/// Defines the Element type of the generic element class.
		/// </summary>
		public IIsogeometricElement ElementType { get; set; }

		/// <summary>
		/// Returns the element type.
		/// </summary>
		IElementType IElement.ElementType => ElementType;

		/// <summary>
		/// The ID of the <see cref="Element"/>.
		/// </summary>
		public int ID { get; set; }

		/// <summary>
		/// Returns an <see cref="IEnumerable{Knot}"/> that define the boundaries of the <see cref="Element"/>.
		/// </summary>
		public IEnumerable<Knot> Knots => KnotsDictionary.Values;

		/// <summary>
		/// <see cref="Dictionary{TKey,TValue}"/> where int is the <see cref="Knot"/> ID.
		/// </summary>
		public Dictionary<int, Knot> KnotsDictionary { get; } = new Dictionary<int, Knot>();

		/// <summary>
		/// The model that the element belongs to.
		/// </summary>
		public Model Model { get; set; }

		/// <summary>
		/// Returns an <see cref="IReadOnlyList{ControlPoint}"/> that belong to the support of the <see cref="Element"/>.
		/// </summary>
		public IReadOnlyList<INode> Nodes
		{
			get
			{
				var a = new List<INode>();
				foreach (var controlPoint in ControlPointsDictionary.Values)
					a.Add(controlPoint);
				return a;
			}
		}

		/// <summary>
		/// The patch that contains the <see cref="Element"/>.
		/// </summary>
		public Patch Patch { get; set; }

		/// <summary>
		/// The patch that contains the <see cref="Element"/> using the <see cref="ISubdomain"/> interface.
		/// </summary>
		public ISubdomain Subdomain => this.Patch;

		/// <summary>
		/// Adds a <see cref="ControlPoints"/> to the <see cref="Element"/>.
		/// </summary>
		/// <param name="controlPoint">A <see cref="ControlPoint"/> object.</param>
		public void AddControlPoint(ControlPoint controlPoint)
		{
			ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
		}

		/// <summary>
		/// Adds a <see cref="IList{ControlPoint}"/> to the <see cref="Element"/>.
		/// </summary>
		/// <param name="controlPoints">An <see cref="IList{T}"/> of <see cref="ControlPoint"/> objects.</param>
		public void AddControlPoints(IList<ControlPoint> controlPoints)
		{
			foreach (ControlPoint controlPoint in controlPoints) AddControlPoint(controlPoint);
		}

		/// <summary>
		/// Adds a <see cref="Knot"/> to the <see cref="Element"/>.
		/// </summary>
		/// <param name="knot">A <see cref="Knot"/> object.</param>
		public void AddKnot(Knot knot)
		{
			KnotsDictionary.Add(knot.ID, knot);
		}

		/// <summary>
		/// Adds a <see cref="IList{Knot}"/> to the <see cref="Element"/>.
		/// </summary>
		/// <param name="knots">An <see cref="IList{T}"/> of <see cref="Knot"/> objects.</param>
		public void AddKnots(IList<Knot> knots)
		{
			foreach (Knot knot in knots) AddKnot(knot);
		}
	}
}
