using MGroup.IGA.Elements;

namespace MGroup.IGA.SupportiveClasses
{
	using System;
	using System.Collections.Generic;
	using System.Linq;

	using MGroup.IGA.Entities;
	using MGroup.LinearAlgebra.Interpolation;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// One-dimensional NURBS shape functions.
	/// </summary>
	public class Nurbs1D
	{
		/// <summary>
		/// Defines an 1D NURBS shape function for an element.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		public Nurbs1D(Element element, ControlPoint[] controlPoints)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Knots.ToArray());

			var parametricGaussPointKsi = Vector.CreateZero(element.Patch.DegreeKsi + 1);
			for (int i = 0; i < element.Patch.DegreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i].Ksi;
			}
			var bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi, parametricGaussPointKsi);
			bsplinesKsi.calculateBSPLinesAndDerivatives();

			int supportKsi = element.Patch.DegreeKsi + 1;
			int numberOfElementControlPoints = supportKsi;

			NurbsValues = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Count);
			NurbsDerivativeValuesKsi = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Count);
			for (int i = 0; i < supportKsi; i++)
			{
				double sumKsi = 0;
				double sumdKsi = 0;

				for (int j = 0; j < numberOfElementControlPoints; j++)
				{
					int indexKsi = controlPoints[j].ID;
					sumKsi += bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[j].WeightFactor;
					sumdKsi += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * controlPoints[j].WeightFactor;
				}
				for (int j = 0; j < numberOfElementControlPoints; j++)
				{
					int indexKsi = controlPoints[j].ID;
					NurbsValues[j, i] = bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[j].WeightFactor / sumKsi;
					NurbsDerivativeValuesKsi[j, i] = controlPoints[j].WeightFactor * (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsi -
						bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsi) / Math.Pow(sumKsi, 2);
				}
			}
		}

		/// <summary>
		/// Defines an 1D NURBS shape function for an edge element.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="edge">The one-dimensional boundary entity that contains the <paramref name="element"/>.</param>
		public Nurbs1D(Element element, IList<ControlPoint> controlPoints, Edge edge)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(edge.Degree, element.Knots.ToArray());

			var parametricGaussPointKsi = Vector.CreateZero(edge.Degree + 1);
			for (int i = 0; i < edge.Degree + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i].Ksi;
			}
			var bsplinesKsi = new BSPLines1D(edge.Degree, edge.KnotValueVector, parametricGaussPointKsi);
			bsplinesKsi.calculateBSPLinesAndDerivatives();

			int supportKsi = edge.Degree + 1;
			int numberOfElementControlPoints = supportKsi;

			NurbsValues = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Count);
			NurbsDerivativeValuesKsi = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Count);

			for (int i = 0; i < supportKsi; i++)
			{
				double sumKsi = 0;
				double sumdKsi = 0;

				for (int j = 0; j < numberOfElementControlPoints; j++)
				{
					int index = controlPoints[j].ID;
					sumKsi += bsplinesKsi.BSPLineValues[index, i] * controlPoints[j].WeightFactor;
					sumdKsi += bsplinesKsi.BSPLineDerivativeValues[index, i] * controlPoints[j].WeightFactor;
				}
				for (int j = 0; j < numberOfElementControlPoints; j++)
				{
					int indexKsi = controlPoints[j].ID;
					NurbsValues[j, i] = bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[j].WeightFactor / sumKsi;
					NurbsDerivativeValuesKsi[j, i] = controlPoints[j].WeightFactor * (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsi -
						bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsi) / Math.Pow(sumKsi, 2);
				}
			}
		}

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public Matrix NurbsDerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public Matrix NurbsValues { get; private set; }
	}
}
