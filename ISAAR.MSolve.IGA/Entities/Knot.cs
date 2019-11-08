namespace ISAAR.MSolve.IGA.Entities
{
	/// <summary>
	/// Defines a Knot.
	/// </summary>
	public class Knot
	{
		/// <summary>
		/// Parametric coordinate Heta of the <see cref="Knot"/>.
		/// </summary>
		public double Heta { get; set; }

		/// <summary>
		/// ID of the <see cref="Knot"/>.
		/// </summary>
		public int ID { get; set; }

		/// <summary>
		/// Parametric coordinate Ksi of the <see cref="Knot"/>.
		/// </summary>
		public double Ksi { get; set; }

		/// <summary>
		/// Parametric coordinate Zeta of the <see cref="Knot"/>.
		/// </summary>
		public double Zeta { get; set; }
	}
}
