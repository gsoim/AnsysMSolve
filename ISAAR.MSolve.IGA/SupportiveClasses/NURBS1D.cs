using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
    /// <summary>
	/// One-dimensional NURBS shape functions.
	/// </summary>
	public class Nurbs1D:IShapeFunction1D
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

            Values = new double[numberOfElementControlPoints, gaussPoints.Length];
			DerivativeValues = new double[numberOfElementControlPoints, gaussPoints.Length];
			for (int i = 0; i < supportKsi; i++)
			{
				double sumKsi = 0;
				double sumdKsi = 0;

				for (int j = 0; j < numberOfElementControlPoints; j++)
				{
					int indexKsi = controlPoints[j].ID;
					sumKsi += bsplinesKsi.Values[indexKsi, i] * controlPoints[j].WeightFactor;
					sumdKsi += bsplinesKsi.DerivativeValues[indexKsi, i] * controlPoints[j].WeightFactor;
				}
				for (int j = 0; j < numberOfElementControlPoints; j++)
				{
					int indexKsi = controlPoints[j].ID;
					Values[j, i] = bsplinesKsi.Values[indexKsi, i] * controlPoints[j].WeightFactor / sumKsi;
					DerivativeValues[j, i] = controlPoints[j].WeightFactor * (bsplinesKsi.DerivativeValues[indexKsi, i] * sumKsi -
						bsplinesKsi.Values[indexKsi, i] * sumdKsi) / Math.Pow(sumKsi, 2);
				}
			}
		}

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValues { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] Values { get; private set; }

		public double[,] SecondDerivativeValues => throw new NotImplementedException();
	}
}
