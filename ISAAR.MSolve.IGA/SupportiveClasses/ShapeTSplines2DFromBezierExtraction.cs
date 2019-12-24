using System;
using System.Collections.Generic;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
    /// <summary>
	/// Two-dimensional T-spline shape functions from Bezier extraction.
	/// </summary>
	public class ShapeTSplines2DFromBezierExtraction:IShapeFunction2D
	{
		/// <summary>
		/// Two-dimensional T-spline shape functions from Bezier extraction for <see cref="TSplineKirchhoffLoveShellElement"/>.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="TSplineKirchhoffLoveShellElement"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		public ShapeTSplines2DFromBezierExtraction(int degreeKsi, int degreeHeta, Matrix extractionOperator, ControlPoint[] controlPoints)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta,
				new List<Knot>
				{
					new Knot(){ID=0,Ksi=-1,Heta = -1,Zeta = 0},
					new Knot(){ID=1,Ksi=-1,Heta = 1,Zeta = 0},
					new Knot(){ID=2,Ksi=1,Heta = -1,Zeta = 0},
					new Knot(){ID=3,Ksi=1,Heta = 1,Zeta = 0}
				});

			var parametricGaussPointKsi = new double[(degreeKsi + 1) * 2];
            for (int i = 0; i < degreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i * (degreeHeta + 1)].Ksi;
			}

			var parametricGaussPointHeta = new double[(degreeHeta + 1) * 2];
            for (int i = 0; i < degreeHeta + 1; i++)
			{
				parametricGaussPointHeta[i] = gaussPoints[i].Heta;
			}

			var knotValueVectorKsi = new double[(degreeKsi + 1) * 2];
			var knotValueVectorHeta = new double[(degreeHeta + 1) * 2];
            for (int i = 0; i < degreeKsi + 1; i++)
			{
				knotValueVectorKsi[i] = -1;
				knotValueVectorKsi[degreeKsi + 1 + i] = 1;
			}
			for (int i = 0; i < degreeHeta + 1; i++)
			{
				knotValueVectorHeta[i] = -1;
				knotValueVectorHeta[degreeHeta + 1 + i] = 1;
			}

			BSPLines1D bernsteinKsi = new BSPLines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bernsteinHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
			bernsteinKsi.calculateBSPLinesAndDerivatives();
			bernsteinHeta.calculateBSPLinesAndDerivatives();

			int supportKsi = degreeKsi + 1;
			int supportHeta = degreeHeta + 1;

			var bKsi = MatrixPart(supportKsi, bernsteinKsi.Values);
			var bdKsi = MatrixPart(supportKsi, bernsteinKsi.DerivativeValues);
			var bddKsi = MatrixPart(supportKsi, bernsteinKsi.SecondDerivativeValues);

			var bheta = MatrixPart(supportHeta, bernsteinHeta.Values);
			var bdheta = MatrixPart(supportHeta, bernsteinHeta.DerivativeValues);
			var bddheta = MatrixPart(supportHeta, bernsteinHeta.SecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);

			Matrix rationalTSplines = extractionOperator * bernsteinShapeFunctions;
			Matrix rationalTSplineDerivativesKsi = extractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix rationalTSplineDerivativesHeta = extractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsi = extractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix rationalTSplineSecondDerivativesHeta = extractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsiHeta = extractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

			Values = new double[controlPoints.Length, supportKsi * supportHeta];
			DerivativeValuesKsi = new double[controlPoints.Length, supportKsi * supportHeta];
			DerivativeValuesHeta = new double[controlPoints.Length, supportKsi * supportHeta];
			SecondDerivativeValuesKsi = new double[controlPoints.Length, supportKsi * supportHeta];
			SecondDerivativeValuesHeta = new double[controlPoints.Length, supportKsi * supportHeta];
			SecondDerivativeValuesKsiHeta = new double[controlPoints.Length, supportKsi * supportHeta];

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

					var index = i * supportHeta + j;

					for (int k = 0; k < controlPoints.Length; k++)
					{
						sumKsiHeta += rationalTSplines[k, index] * controlPoints[k].WeightFactor;
						sumdKsiHeta += rationalTSplineDerivativesKsi[k, index] * controlPoints[k].WeightFactor;
						sumKsidHeta += rationalTSplineDerivativesHeta[k, index] * controlPoints[k].WeightFactor;
						sumdKsidKsi += rationalTSplineSecondDerivativesKsi[k, index] * controlPoints[k].WeightFactor;
						sumdHetadHeta += rationalTSplineSecondDerivativesHeta[k, index] * controlPoints[k].WeightFactor;
						sumdKsidHeta += rationalTSplineSecondDerivativesKsiHeta[k, index] * controlPoints[k].WeightFactor;
					}

					for (int k = 0; k < controlPoints.Length; k++)
					{
						Values[k, index] = rationalTSplines[k, index] * controlPoints[k].WeightFactor / sumKsiHeta;
						DerivativeValuesKsi[k, index] = (rationalTSplineDerivativesKsi[k, index] * sumKsiHeta -
																rationalTSplines[k, index] * sumdKsiHeta) /
															   Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						DerivativeValuesHeta[k, index] = (rationalTSplineDerivativesHeta[k, index] * sumKsiHeta -
																 rationalTSplines[k, index] * sumKsidHeta) /
																Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						SecondDerivativeValuesKsi[k, index] = (rationalTSplineSecondDerivativesKsi[k, index] / sumKsiHeta -
																	  2 * rationalTSplineDerivativesKsi[k, index] * sumdKsiHeta /
																	  Math.Pow(sumKsiHeta, 2) -
																	  rationalTSplines[k, index] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
																	  2 * rationalTSplines[k, index] * Math.Pow(sumdKsiHeta, 2) /
																	  Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						SecondDerivativeValuesHeta[k, index] = (rationalTSplineSecondDerivativesHeta[k, index] / sumKsiHeta -
																	   2 * rationalTSplineDerivativesHeta[k, index] * sumKsidHeta /
																	   Math.Pow(sumKsiHeta, 2) -
																	   rationalTSplines[k, index] * sumdHetadHeta /
																	   Math.Pow(sumKsiHeta, 2) +
																	   2 * rationalTSplines[k, index] * Math.Pow(sumKsidHeta, 2) /
																	   Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						SecondDerivativeValuesKsiHeta[k, index] = (rationalTSplineSecondDerivativesKsiHeta[k, index] / sumKsiHeta -
																		  rationalTSplineDerivativesKsi[k, index] * sumKsidHeta /
																		  Math.Pow(sumKsiHeta, 2) -
																		  rationalTSplineDerivativesHeta[k, index] * sumdKsiHeta /
																		  Math.Pow(sumKsiHeta, 2) -
																		  rationalTSplines[k, index] * sumdKsidHeta /
																		  Math.Pow(sumKsiHeta, 2) +
																		  2 * rationalTSplines[k, index] * sumdKsiHeta * sumKsidHeta /
																		  Math.Pow(sumKsiHeta, 3)) *
																		 controlPoints[k].WeightFactor;
					}
				}
			}
		}
        
		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function derivatives per axis Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function derivatives per axis Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function second derivatives per axis Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function second derivatives per axis Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function mixed second derivatives per axis Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsiHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] Values { get; private set; }

		private static Matrix KroneckerProduct(Matrix A, Matrix B)
		{
			Matrix C = Matrix.CreateZero(A.NumRows * B.NumRows, A.NumColumns * B.NumColumns);
			for (int rowAIndex = 0; rowAIndex < A.NumRows; rowAIndex++)
			{
				for (int rowBIndex = 0; rowBIndex < B.NumRows; rowBIndex++)
				{
					for (int columnAIndex = 0; columnAIndex < A.NumColumns; columnAIndex++)
					{
						for (int columnBIndex = 0; columnBIndex < B.NumColumns; columnBIndex++)
						{
							C[rowAIndex * B.NumRows + rowBIndex, columnAIndex * B.NumColumns + columnBIndex] =
								A[rowAIndex, columnAIndex] * B[rowBIndex, columnBIndex];
						}
					}
				}
			}

			return C;
		}

		private static Matrix MatrixPart(int support, double[,] matrix)
		{
			var A = Matrix.CreateZero(support, support);
			for (int i = 0; i < support; i++)
			{
				for (int j = 0; j < support; j++)
				{
					A[i, j] = matrix[i, j];
				}
			}

			return A;
		}

		private static Matrix MatrixPart(int support1, int support2, double[,] matrix)
		{
			var A = Matrix.CreateZero(support1, support2);
			for (int i = 0; i < support1; i++)
			{
				for (int j = 0; j < support2; j++)
				{
					A[i, j] = matrix[i, j];
				}
			}

			return A;
		}
	}
}
