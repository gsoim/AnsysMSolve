using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
    /// <summary>
	/// Two-dimensional NURBS shape functions.
	/// </summary>
	public class Nurbs2D
	{
		/// <summary>
		/// Defines a 2D NURBS shape function for an element given the per axis gauss point coordinates.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NURBSElement2D"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="parametricGaussPointKsi">An <see cref="IVector"/> containing Gauss points of axis Ksi.</param>
		/// <param name="parametricGaussPointHeta">An <see cref="IVector"/> containing Gauss points of axis Heta.</param>
		public Nurbs2D(int degreeKsi, double[] knotValueVectorKsi,
                       int degreeHeta, double[] knotValueVectorHeta,
                       ControlPoint[] controlPoints, double[] parametricGaussPointKsi,
                       double[] parametricGaussPointHeta)
		{
			var parametricPointsCount = parametricGaussPointKsi.Length * parametricGaussPointHeta.Length;
            var numberOfControlPointsHeta = knotValueVectorHeta.Length - degreeHeta - 1;
			BSPLines1D bsplinesKsi = new BSPLines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
				parametricGaussPointHeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();

			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;
			int numberOfElementControlPoints = (degreeKsi + 1) * (degreeHeta + 1);

			NurbsValues = new double[numberOfElementControlPoints, parametricPointsCount];
			NurbsDerivativeValuesKsi = new double[numberOfElementControlPoints, parametricPointsCount];
			NurbsDerivativeValuesHeta = new double[numberOfElementControlPoints, parametricPointsCount];
			NurbsSecondDerivativeValueKsi = new double[numberOfElementControlPoints, parametricPointsCount];
			NurbsSecondDerivativeValueHeta = new double[numberOfElementControlPoints, parametricPointsCount];
			NurbsSecondDerivativeValueKsiHeta = new double[numberOfElementControlPoints, parametricPointsCount];

			for (int i = 0; i < supportKsi; i++)
			{
				for (int j = 0; j < supportHeta; j++)
				{
					double sumKsiHeta = 0;
					double sumdKsiHeta = 0;
					double sumKsidHeta = 0;
					double sumdKsidKsi = 0;
					double sumdHetadHeta = 0;
					double sumdKsidHeta = 0;

					for (int k = 0; k < numberOfElementControlPoints; k++)
					{
						int indexKsi = controlPoints[k].ID / numberOfControlPointsHeta;
						int indexHeta = controlPoints[k].ID % numberOfControlPointsHeta;
						sumKsiHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
									  bsplinesHeta.BSPLineValues[indexHeta, j] *
									  controlPoints[k].WeightFactor;
						sumdKsiHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
									   bsplinesHeta.BSPLineValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumKsidHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
									   bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumdKsidKsi += bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] *
									   bsplinesHeta.BSPLineValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumdHetadHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
										 bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] *
										 controlPoints[k].WeightFactor;
						sumdKsidHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
										bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
										controlPoints[k].WeightFactor;
					}

					for (int k = 0; k < numberOfElementControlPoints; k++)
					{
						int indexKsi = controlPoints[k].ID / numberOfControlPointsHeta;
						int indexHeta = controlPoints[k].ID % numberOfControlPointsHeta;

						NurbsValues[k, i * supportHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] *
							bsplinesHeta.BSPLineValues[indexHeta, j] *
							controlPoints[k].WeightFactor / sumKsiHeta;

						NurbsDerivativeValuesKsi[k, i * supportHeta + j] =
							bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHeta -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHeta) / Math.Pow(sumKsiHeta, 2);

						NurbsDerivativeValuesHeta[k, i * supportHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHeta -
							 bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHeta) / Math.Pow(sumKsiHeta, 2);

						NurbsSecondDerivativeValueKsi[k, i * supportHeta + j] =
							bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHeta -
							 2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHeta /
							 Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) /
							 Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueHeta[k, i * supportHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHeta -
							 2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHeta /
							 Math.Pow(sumKsiHeta, 2) -
							 bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHeta, 2) /
							 Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueKsiHeta[k, i * supportHeta + j] =
							controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
							 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] / sumKsiHeta -
							 bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
							 bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] *
							 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
							 sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumdKsidHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
							 sumdKsiHeta * sumKsidHeta / Math.Pow(sumKsiHeta, 3));
					}
				}
			}
		}

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives per Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] NurbsDerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives per Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] NurbsDerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function mixed second derivatives per Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] NurbsSecondDerivativeValueHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] NurbsSecondDerivativeValueKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] NurbsSecondDerivativeValueKsiHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] NurbsValues { get; private set; }
	}
}
