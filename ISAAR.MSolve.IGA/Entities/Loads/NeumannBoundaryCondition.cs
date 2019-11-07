namespace MGroup.IGA.Entities.Loads
{
	using System.Collections.Generic;

	/// <summary>
	/// Delegate that calculates the load distribution that can vary depending on the position.
	/// </summary>
	/// <param name="x">Cartesian coordinate x of the load.</param>
	/// <param name="y">Cartesian coordinate y of the load.</param>
	/// <param name="z">Cartesian coordinate z of the load.</param>
	/// <returns>The of the boundary condition for the point <paramref name="x"/>,<paramref name="y"/>,<paramref name="z"/></returns>
	public delegate double[] Value(double x, double y, double z);

	/// <summary>
	/// Calculates the load depending on the the type of boundary condition and geometrical entity enforced.
	/// </summary>
	public class LoadProvider
	{
		/// <summary>
		/// Calculates Neumann load on an edge.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> that can handle a <see cref="NeumannBoundaryCondition"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="neumann">The <see cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> LoadNeumann(Element element, Edge edge, NeumannBoundaryCondition neumann) => element.ElementType.CalculateLoadingCondition(element, edge, neumann);

		/// <summary>
		/// Calculates Neumann load on a face.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> that can handle a <see cref="NeumannBoundaryCondition"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="NeumannBoundaryCondition"/> was applied to.</param>
		/// <param name="neumann">The <see cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="NeumannBoundaryCondition"/>.</returns>
		public Dictionary<int, double> LoadNeumann(Element element, Face face, NeumannBoundaryCondition neumann) => element.ElementType.CalculateLoadingCondition(element, face, neumann);

		/// <summary>
		/// Calculates Pressure load on an edge.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> that can handle a <see cref="PressureBoundaryCondition"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> LoadPressure(Element element, Edge edge, PressureBoundaryCondition pressure) => element.ElementType.CalculateLoadingCondition(element, edge, pressure);

		/// <summary>
		/// Calculates Pressure load on a face.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> that can handle a <see cref="PressureBoundaryCondition"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="NeumannBoundaryCondition"/> was applied to.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load due to the <see cref="PressureBoundaryCondition"/>.</returns>
		public Dictionary<int, double> LoadPressure(Element element, Face face, PressureBoundaryCondition pressure) => element.ElementType.CalculateLoadingCondition(element, face, pressure);
	}

	/// <summary>
	/// Neumann loading boundary condition. Works as a distributed load.
	/// </summary>
	public class NeumannBoundaryCondition : LoadingCondition
	{
		/// <summary>
		/// Defines a Neumann BoundaryCondition.
		/// </summary>
		/// <param name="neumannValue">A delegate that calculates the value of the Neumann BC at a given point.</param>
		public NeumannBoundaryCondition(Value neumannValue)
		{
			this.Value = neumannValue;
		}

		/// <summary>
		/// Describes the distribution of the load. See <see cref="Loads.Value"/>.
		/// </summary>
		public Value Value { get; private set; }
	}
}
