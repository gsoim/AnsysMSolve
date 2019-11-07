namespace MGroup.IGA.Interfaces
{
	using System.Collections.Generic;

	using MGroup.IGA.Entities;
	using MGroup.IGA.Entities.Loads;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.Discretization.Interfaces;

	/// <summary>
	/// Isogeometric element interface. Implements <see cref="IElementType"/>.
	/// </summary>
	public interface IIsogeometricElement : IElementType
	{
		/// <summary>
		/// Retries the element dimensions.
		/// </summary>
		ElementDimensions ElementDimensions { get; }

		/// <summary>
		/// Retrieves the element ID.
		/// </summary>
		int ID { get; }

		/// <summary>
		/// Calculates the knot displacements for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An isogeometric <see cref="Element"/>.</param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements);

		/// <summary>
		/// Neumann loading condition imposed on edges.
		/// </summary>
		/// <param name="element">An isogeometric <see cref="Element"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="neumann"><inheritdoc cref="NeumannBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="NeumannBoundaryCondition"/>.</returns>
		Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann);

		/// <summary>
		/// Neumann loading condition imposed on faces.
		/// </summary>
		/// <param name="element">An isogeometric <see cref="Element"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="NeumannBoundaryCondition"/> was applied to.</param>
		/// <param name="neumann">A <see cref="Dictionary{TKey,TValue}"/> where integer values denote the degree of freedom that has a value double load value due to the enforcement of the <see cref="NeumannBoundaryCondition"/>.</param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load.</returns>
		Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann);

		/// <summary>
		/// Pressure loading condition imposed on edges.
		/// </summary>
		/// <param name="element">An isogeometric <see cref="Element"/>.</param>
		/// <param name="edge">An one dimensional boundary entity. For more info see <see cref="Edge"/>.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load.</returns>
		Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure);

		/// <summary>
		/// Pressure loading condition imposed on faces.
		/// </summary>
		/// <param name="element">An isogeometric <see cref="Element"/>.</param>
		/// <param name="face">The <see cref="Face"/> that the <see cref="PressureBoundaryCondition"/> was applied to.</param>
		/// <param name="pressure"><inheritdoc cref="PressureBoundaryCondition"/></param>
		/// <returns>A <see cref="Dictionary{TKey,TValue}"/> whose keys are the numbering of the degree of freedom and values are the magnitude of the load.</returns>
		Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure);
	}
}
