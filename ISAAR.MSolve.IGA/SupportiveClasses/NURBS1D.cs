using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
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
		public Nurbs1D(int degree, double[] KnotValueVector, ControlPoint[] controlPoints, GaussLegendrePoint3D[] gaussPoints)
        {
            var parametricGaussPointKsi = gaussPoints.Select(gp => gp.Ksi).ToArray();
			var bsplinesKsi = new BSPLines1D(degree, KnotValueVector, parametricGaussPointKsi);
			bsplinesKsi.calculateBSPLinesAndDerivatives();

			int supportKsi = degree + 1;
			int numberOfElementControlPoints = supportKsi;

			NurbsValues = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Length);
			NurbsDerivativeValuesKsi = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Length);
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
