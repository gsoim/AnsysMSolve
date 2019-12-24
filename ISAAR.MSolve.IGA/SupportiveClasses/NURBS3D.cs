using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
    public class Nurbs3D: IShapeFunction3D
	{
		/// <summary>
		/// Define 3D  NURBS shape function without needing an element definition.
		/// </summary>
		/// <param name="numberOfControlPointsKsi">Number of control points along parametric axis Ksi.</param>
		/// <param name="numberOfControlPointsHeta">Number of control points along parametric axis Heta.</param>
		/// <param name="numberOfControlPointsZeta">Number of control points along parametric axis Zeta.</param>
		/// <param name="degreeKsi">Polynomial degree of the parametric axis Ksi.</param>
		/// <param name="degreeHeta">Polynomial degree of the parametric axis Heta.</param>
		/// <param name="degreeZeta">Polynomial degree of the parametric axis Zeta.</param>
		/// <param name="knotValueVectorKsi">Knot value vector of the parametric axis Ksi.</param>
		/// <param name="knotValueVectorHeta">Knot value vector of the parametric axis Heta.</param>
		/// <param name="knotValueVectorZeta">Knot value vector of the parametric axis Zeta.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="gaussPoints">A <see cref="List{T}"/> containing the control points where the shape functions will be evaluated.</param>
		public Nurbs3D(int numberOfControlPointsKsi, int numberOfControlPointsHeta, int numberOfControlPointsZeta,
			int degreeKsi, int degreeHeta, int degreeZeta, double[] knotValueVectorKsi,
			double[] knotValueVectorHeta, double[] knotValueVectorZeta, ControlPoint[] controlPoints,
			GaussLegendrePoint3D[] gaussPoints)
		{
			var numberOfGaussPoints = gaussPoints.Length;
			var parametricGaussPointKsi = new double[degreeKsi + 1];
			for (int i = 0; i < degreeKsi + 1; i++)
				parametricGaussPointKsi[i] = gaussPoints[i * (degreeZeta + 1) * (degreeHeta + 1)].Ksi;

			var parametricGaussPointHeta = new double[degreeHeta + 1];
			for (int i = 0; i < degreeHeta + 1; i++)
				parametricGaussPointHeta[i] = gaussPoints[i * (degreeZeta + 1)].Heta;

			var parametricGaussPointZeta = new double[degreeZeta + 1];
			for (int i = 0; i < degreeZeta + 1; i++)
				parametricGaussPointZeta[i] = gaussPoints[i].Zeta;

			BSPLines1D bsplinesKsi =
				new BSPLines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
				parametricGaussPointHeta);
			BSPLines1D bsplinesZeta = new BSPLines1D(degreeZeta, knotValueVectorZeta,
				parametricGaussPointZeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();
			bsplinesZeta.calculateBSPLinesAndDerivatives();

			int supportKsi = degreeKsi + 1;
			int supportHeta = degreeHeta + 1;
			int supportZeta = degreeZeta + 1;
			int numberOfElementControlPoints = supportKsi * supportHeta * supportZeta;

			Values = new double[numberOfElementControlPoints, numberOfGaussPoints];
			DerivativeValuesKsi = new double[numberOfElementControlPoints, numberOfGaussPoints];
			DerivativeValuesHeta = new double[numberOfElementControlPoints, numberOfGaussPoints];
			DerivativeValuesZeta = new double[numberOfElementControlPoints, numberOfGaussPoints];

			for (int i = 0; i < supportKsi; i++)
			{
				for (int j = 0; j < supportHeta; j++)
				{
					for (int k = 0; k < supportZeta; k++)
					{
						double sumKsiHetaZeta = 0;
						double sumdKsiHetaZeta = 0;
						double sumKsidHetaZeta = 0;
						double sumKsiHetadZeta = 0;

						for (int m = 0; m < numberOfElementControlPoints; m++)
						{
							int indexKsi = controlPoints[m].ID /
										   (numberOfControlPointsHeta *
											numberOfControlPointsZeta);
							int indexHeta = controlPoints[m].ID %
											(numberOfControlPointsHeta *
											 numberOfControlPointsZeta) /
											numberOfControlPointsZeta;
							int indexZeta = controlPoints[m].ID %
											(numberOfControlPointsHeta *
											 numberOfControlPointsZeta) %
											numberOfControlPointsZeta;

							sumKsiHetaZeta += bsplinesKsi.Values[indexKsi, i] *
											  bsplinesHeta.Values[indexHeta, j] *
											  bsplinesZeta.Values[indexZeta, k] *
											  controlPoints[m].WeightFactor;

							sumdKsiHetaZeta += bsplinesKsi.DerivativeValues[indexKsi, i] *
											   bsplinesHeta.Values[indexHeta, j] *
											   bsplinesZeta.Values[indexZeta, k] *
											   controlPoints[m].WeightFactor;

							sumKsidHetaZeta += bsplinesKsi.Values[indexKsi, i] *
											   bsplinesHeta.DerivativeValues[indexHeta, j] *
											   bsplinesZeta.Values[indexZeta, k] *
											   controlPoints[m].WeightFactor;

							sumKsiHetadZeta += bsplinesKsi.Values[indexKsi, i] *
											   bsplinesHeta.Values[indexHeta, j] *
											   bsplinesZeta.DerivativeValues[indexZeta, k] *
											   controlPoints[m].WeightFactor;
						}

						for (int m = 0; m < numberOfElementControlPoints; m++)
						{
							int indexKsi = controlPoints[m].ID /
										   (numberOfControlPointsHeta *
											numberOfControlPointsZeta);
							int indexHeta = controlPoints[m].ID %
											(numberOfControlPointsHeta *
											 numberOfControlPointsZeta) /
											numberOfControlPointsZeta;
							int indexZeta = controlPoints[m].ID %
											(numberOfControlPointsHeta *
											 numberOfControlPointsZeta) %
											numberOfControlPointsZeta;

							Values[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								bsplinesKsi.Values[indexKsi, i] *
								bsplinesHeta.Values[indexHeta, j] *
								bsplinesZeta.Values[indexZeta, k] *
								controlPoints[m].WeightFactor / sumKsiHetaZeta;

							DerivativeValuesKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								(bsplinesKsi.DerivativeValues[indexKsi, i] * sumKsiHetaZeta -
								 bsplinesKsi.Values[indexKsi, i] * sumdKsiHetaZeta) *
								bsplinesHeta.Values[indexHeta, j] *
								bsplinesZeta.Values[indexZeta, k] *
								controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

							DerivativeValuesHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								bsplinesKsi.Values[indexKsi, i] *
								(bsplinesHeta.DerivativeValues[indexHeta, j] * sumKsiHetaZeta -
								 bsplinesHeta.Values[indexHeta, j] * sumKsidHetaZeta) *
								bsplinesZeta.Values[indexZeta, k] *
								controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

							DerivativeValuesZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
								bsplinesKsi.Values[indexKsi, i] *
								bsplinesHeta.Values[indexHeta, j] *
								(bsplinesZeta.DerivativeValues[indexZeta, k] * sumKsiHetaZeta -
								 bsplinesZeta.Values[indexZeta, k] * sumKsiHetadZeta) *
								controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);
						}
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
		/// <see cref="Matrix"/> containing NURBS shape function derivatives per Zeta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesZeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function mixed second derivatives per Heta and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesHetaZeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function mixed second derivatives per Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsiHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function mixed second derivatives per Ksi and Zeta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsiZeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Zeta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesZeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] Values { get; private set; }
	}
}
