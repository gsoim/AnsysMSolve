using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
    /// <summary>
	/// Two-dimensional NURBS shape functions.
	/// </summary>
	public class Nurbs2D:IShapeFunction2D
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

			Values = new double[numberOfElementControlPoints, parametricPointsCount];
			DerivativeValuesKsi = new double[numberOfElementControlPoints, parametricPointsCount];
			DerivativeValuesHeta = new double[numberOfElementControlPoints, parametricPointsCount];
			SecondDerivativeValuesKsi = new double[numberOfElementControlPoints, parametricPointsCount];
			SecondDerivativeValuesHeta = new double[numberOfElementControlPoints, parametricPointsCount];
			SecondDerivativeValuesKsiHeta = new double[numberOfElementControlPoints, parametricPointsCount];

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
						sumKsiHeta += bsplinesKsi.Values[indexKsi, i] *
									  bsplinesHeta.Values[indexHeta, j] *
									  controlPoints[k].WeightFactor;
						sumdKsiHeta += bsplinesKsi.DerivativeValues[indexKsi, i] *
									   bsplinesHeta.Values[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumKsidHeta += bsplinesKsi.Values[indexKsi, i] *
									   bsplinesHeta.DerivativeValues[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumdKsidKsi += bsplinesKsi.SecondDerivativeValues[indexKsi, i] *
									   bsplinesHeta.Values[indexHeta, j] *
									   controlPoints[k].WeightFactor;
						sumdHetadHeta += bsplinesKsi.Values[indexKsi, i] *
										 bsplinesHeta.SecondDerivativeValues[indexHeta, j] *
										 controlPoints[k].WeightFactor;
						sumdKsidHeta += bsplinesKsi.DerivativeValues[indexKsi, i] *
										bsplinesHeta.DerivativeValues[indexHeta, j] *
										controlPoints[k].WeightFactor;
					}

					for (int k = 0; k < numberOfElementControlPoints; k++)
					{
						int indexKsi = controlPoints[k].ID / numberOfControlPointsHeta;
						int indexHeta = controlPoints[k].ID % numberOfControlPointsHeta;

						Values[k, i * supportHeta + j] =
							bsplinesKsi.Values[indexKsi, i] *
							bsplinesHeta.Values[indexHeta, j] *
							controlPoints[k].WeightFactor / sumKsiHeta;

						DerivativeValuesKsi[k, i * supportHeta + j] =
							bsplinesHeta.Values[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.DerivativeValues[indexKsi, i] * sumKsiHeta -
							 bsplinesKsi.Values[indexKsi, i] * sumdKsiHeta) / Math.Pow(sumKsiHeta, 2);

						DerivativeValuesHeta[k, i * supportHeta + j] =
							bsplinesKsi.Values[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.DerivativeValues[indexHeta, j] * sumKsiHeta -
							 bsplinesHeta.Values[indexHeta, j] * sumKsidHeta) / Math.Pow(sumKsiHeta, 2);

						SecondDerivativeValuesKsi[k, i * supportHeta + j] =
							bsplinesHeta.Values[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.SecondDerivativeValues[indexKsi, i] / sumKsiHeta -
							 2 * bsplinesKsi.DerivativeValues[indexKsi, i] * sumdKsiHeta /
							 Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.Values[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.Values[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) /
							 Math.Pow(sumKsiHeta, 3));

						SecondDerivativeValuesHeta[k, i * supportHeta + j] =
							bsplinesKsi.Values[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.SecondDerivativeValues[indexHeta, j] / sumKsiHeta -
							 2 * bsplinesHeta.DerivativeValues[indexHeta, j] * sumKsidHeta /
							 Math.Pow(sumKsiHeta, 2) -
							 bsplinesHeta.Values[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesHeta.Values[indexHeta, j] * Math.Pow(sumKsidHeta, 2) /
							 Math.Pow(sumKsiHeta, 3));

						SecondDerivativeValuesKsiHeta[k, i * supportHeta + j] =
							controlPoints[k].WeightFactor *
							(bsplinesKsi.DerivativeValues[indexKsi, i] *
							 bsplinesHeta.DerivativeValues[indexHeta, j] / sumKsiHeta -
							 bsplinesKsi.DerivativeValues[indexKsi, i] *
							 bsplinesHeta.Values[indexHeta, j] *
							 sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.Values[indexKsi, i] *
							 bsplinesHeta.DerivativeValues[indexHeta, j] *
							 sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.Values[indexKsi, i] * bsplinesHeta.Values[indexHeta, j] *
							 sumdKsidHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.Values[indexKsi, i] * bsplinesHeta.Values[indexHeta, j] *
							 sumdKsiHeta * sumKsidHeta / Math.Pow(sumKsiHeta, 3));
					}
				}
			}
		}

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives per Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives per Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function mixed second derivatives per Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsiHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] Values { get; private set; }
	}
}
