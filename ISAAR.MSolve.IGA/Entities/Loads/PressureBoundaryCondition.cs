namespace MGroup.IGA.Entities.Loads
{
	/// <summary>
	/// Pressure loading boundary condition. Works as a load normal to the boundary entity applied on.
	/// </summary>
	public class PressureBoundaryCondition : LoadingCondition
	{
		/// <summary>
		/// Describes the distribution of the load. See <see cref="Loads.Value"/>.
		/// </summary>
		public double Value { get; private set; }

		/// <summary>
		/// Defines a Pressure BoundaryCondition
		/// </summary>
		/// <param name="pressureValue">A delegate that calculates the value of the Neumann BC at a given point.</param>
		public PressureBoundaryCondition(double pressureValue)
		{
			this.Value = pressureValue;
		}
	}
}
