namespace MGroup.IGA.SupportiveClasses
{
	using System;
	using System.Collections.Generic;

	using MGroup.IGA.Elements;
	using MGroup.IGA.Entities;
	using MGroup.LinearAlgebra.Interpolation;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	/// <summary>
	/// Two-dimensional T-spline shape functions from Bezier extraction.
	/// </summary>
	public class ShapeTSplines2DFromBezierExtraction
	{
		/// <summary>
		/// Two-dimensional T-spline shape functions from Bezier extraction for <see cref="TSplineKirchhoffLoveShellElement"/>.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="TSplineKirchhoffLoveShellElement"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		public ShapeTSplines2DFromBezierExtraction(TSplineKirchhoffLoveShellElement element, ControlPoint[] controlPoints)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta,
				new List<Knot>
				{
					new Knot(){ID=0,Ksi=-1,Heta = -1,Zeta = 0},
					new Knot(){ID=1,Ksi=-1,Heta = 1,Zeta = 0},
					new Knot(){ID=2,Ksi=1,Heta = -1,Zeta = 0},
					new Knot(){ID=3,Ksi=1,Heta = 1,Zeta = 0}
				});

			var parametricGaussPointKsi = Vector.CreateZero(element.DegreeKsi + 1);
			for (int i = 0; i < element.DegreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i * (element.DegreeHeta + 1)].Ksi;
			}

			var parametricGaussPointHeta = Vector.CreateZero(element.DegreeHeta + 1);
			for (int i = 0; i < element.DegreeHeta + 1; i++)
			{
				parametricGaussPointHeta[i] = gaussPoints[i].Heta;
			}

			Vector knotValueVectorKsi = Vector.CreateZero((element.DegreeKsi + 1) * 2);
			Vector knotValueVectorHeta = Vector.CreateZero((element.DegreeHeta + 1) * 2);
			for (int i = 0; i < element.DegreeKsi + 1; i++)
			{
				knotValueVectorKsi[i] = -1;
				knotValueVectorKsi[element.DegreeKsi + 1 + i] = 1;
			}
			for (int i = 0; i < element.DegreeHeta + 1; i++)
			{
				knotValueVectorHeta[i] = -1;
				knotValueVectorHeta[element.DegreeHeta + 1 + i] = 1;
			}

			BSPLines1D bernsteinKsi = new BSPLines1D(element.DegreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bernsteinHeta = new BSPLines1D(element.DegreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
			bernsteinKsi.calculateBSPLinesAndDerivatives();
			bernsteinHeta.calculateBSPLinesAndDerivatives();

			int supportKsi = element.DegreeKsi + 1;
			int supportHeta = element.DegreeHeta + 1;

			var bKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);

			Matrix rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

			TSplineValues = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineDerivativeValuesKsi = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineDerivativeValuesHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsi = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsiHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);

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
						TSplineValues[k, index] = rationalTSplines[k, index] * controlPoints[k].WeightFactor / sumKsiHeta;
						TSplineDerivativeValuesKsi[k, index] = (rationalTSplineDerivativesKsi[k, index] * sumKsiHeta -
																rationalTSplines[k, index] * sumdKsiHeta) /
															   Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						TSplineDerivativeValuesHeta[k, index] = (rationalTSplineDerivativesHeta[k, index] * sumKsiHeta -
																 rationalTSplines[k, index] * sumKsidHeta) /
																Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsi[k, index] = (rationalTSplineSecondDerivativesKsi[k, index] / sumKsiHeta -
																	  2 * rationalTSplineDerivativesKsi[k, index] * sumdKsiHeta /
																	  Math.Pow(sumKsiHeta, 2) -
																	  rationalTSplines[k, index] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
																	  2 * rationalTSplines[k, index] * Math.Pow(sumdKsiHeta, 2) /
																	  Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueHeta[k, index] = (rationalTSplineSecondDerivativesHeta[k, index] / sumKsiHeta -
																	   2 * rationalTSplineDerivativesHeta[k, index] * sumKsidHeta /
																	   Math.Pow(sumKsiHeta, 2) -
																	   rationalTSplines[k, index] * sumdHetadHeta /
																	   Math.Pow(sumKsiHeta, 2) +
																	   2 * rationalTSplines[k, index] * Math.Pow(sumKsidHeta, 2) /
																	   Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsiHeta[k, index] = (rationalTSplineSecondDerivativesKsiHeta[k, index] / sumKsiHeta -
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
		/// Two-dimensional T-spline shape functions from Bezier extraction for <see cref="TSplineKirchhoffLoveShellElement"/>.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="TSplineKirchhoffLoveShellElement"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="parametricGaussPointKsi">An <see cref="IVector"/> containing Gauss points of axis Ksi.</param>
		/// <param name="parametricGaussPointHeta">An <see cref="IVector"/> containing Gauss points of axis Heta.</param>
		public ShapeTSplines2DFromBezierExtraction(TSplineKirchhoffLoveShellElement element, ControlPoint[] controlPoints, Vector parametricGaussPointKsi, Vector parametricGaussPointHeta)
		{
			Vector knotValueVectorKsi = Vector.CreateZero((element.DegreeKsi + 1) * 2);
			Vector knotValueVectorHeta = Vector.CreateZero((element.DegreeHeta + 1) * 2);
			for (int i = 0; i < element.DegreeKsi + 1; i++)
			{
				knotValueVectorKsi[i] = -1;
				knotValueVectorKsi[element.DegreeKsi + 1 + i] = 1;
			}
			for (int i = 0; i < element.DegreeHeta + 1; i++)
			{
				knotValueVectorHeta[i] = -1;
				knotValueVectorHeta[element.DegreeHeta + 1 + i] = 1;
			}

			BSPLines1D bernsteinKsi = new BSPLines1D(element.DegreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bernsteinHeta = new BSPLines1D(element.DegreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
			bernsteinKsi.calculateBSPLinesAndDerivatives();
			bernsteinHeta.calculateBSPLinesAndDerivatives();

			int supportElementKsi = element.DegreeKsi + 1;
			int supportElementHeta = element.DegreeHeta + 1;
			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;

			var bKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);

			Matrix rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

			TSplineValues = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineDerivativeValuesKsi = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineDerivativeValuesHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsi = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsiHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);

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
						TSplineValues[k, index] = rationalTSplines[k, index] * controlPoints[k].WeightFactor / sumKsiHeta;
						TSplineDerivativeValuesKsi[k, index] = (rationalTSplineDerivativesKsi[k, index] * sumKsiHeta -
																rationalTSplines[k, index] * sumdKsiHeta) /
															   Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						TSplineDerivativeValuesHeta[k, index] = (rationalTSplineDerivativesHeta[k, index] * sumKsiHeta -
																 rationalTSplines[k, index] * sumKsidHeta) /
																Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsi[k, index] = (rationalTSplineSecondDerivativesKsi[k, index] / sumKsiHeta -
																	  2 * rationalTSplineDerivativesKsi[k, index] * sumdKsiHeta /
																	  Math.Pow(sumKsiHeta, 2) -
																	  rationalTSplines[k, index] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
																	  2 * rationalTSplines[k, index] * Math.Pow(sumdKsiHeta, 2) /
																	  Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueHeta[k, index] = (rationalTSplineSecondDerivativesHeta[k, index] / sumKsiHeta -
																	   2 * rationalTSplineDerivativesHeta[k, index] * sumKsidHeta /
																	   Math.Pow(sumKsiHeta, 2) -
																	   rationalTSplines[k, index] * sumdHetadHeta /
																	   Math.Pow(sumKsiHeta, 2) +
																	   2 * rationalTSplines[k, index] * Math.Pow(sumKsidHeta, 2) /
																	   Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsiHeta[k, index] = (rationalTSplineSecondDerivativesKsiHeta[k, index] / sumKsiHeta -
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
		/// Two-dimensional T-spline shape functions from Bezier extraction for <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="parametricGaussPointKsi">An <see cref="IVector"/> containing Gauss points of axis Ksi.</param>
		/// <param name="parametricGaussPointHeta">An <see cref="IVector"/> containing Gauss points of axis Heta.</param>
		public ShapeTSplines2DFromBezierExtraction(TSplineKirchhoffLoveShellElementMaterial element, ControlPoint[] controlPoints, Vector parametricGaussPointKsi, Vector parametricGaussPointHeta)
		{
			Vector knotValueVectorKsi = Vector.CreateZero((element.DegreeKsi + 1) * 2);
			Vector knotValueVectorHeta = Vector.CreateZero((element.DegreeHeta + 1) * 2);
			for (int i = 0; i < element.DegreeKsi + 1; i++)
			{
				knotValueVectorKsi[i] = -1;
				knotValueVectorKsi[element.DegreeKsi + 1 + i] = 1;
			}
			for (int i = 0; i < element.DegreeHeta + 1; i++)
			{
				knotValueVectorHeta[i] = -1;
				knotValueVectorHeta[element.DegreeHeta + 1 + i] = 1;
			}

			BSPLines1D bernsteinKsi = new BSPLines1D(element.DegreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bernsteinHeta = new BSPLines1D(element.DegreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
			bernsteinKsi.calculateBSPLinesAndDerivatives();
			bernsteinHeta.calculateBSPLinesAndDerivatives();

			int supportElementKsi = element.DegreeKsi + 1;
			int supportElementHeta = element.DegreeHeta + 1;
			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;

			var bKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);

			Matrix rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

			TSplineValues = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineDerivativeValuesKsi = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineDerivativeValuesHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsi = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsiHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);

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
						TSplineValues[k, index] = rationalTSplines[k, index] * controlPoints[k].WeightFactor / sumKsiHeta;
						TSplineDerivativeValuesKsi[k, index] = (rationalTSplineDerivativesKsi[k, index] * sumKsiHeta -
																rationalTSplines[k, index] * sumdKsiHeta) /
															   Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						TSplineDerivativeValuesHeta[k, index] = (rationalTSplineDerivativesHeta[k, index] * sumKsiHeta -
																 rationalTSplines[k, index] * sumKsidHeta) /
																Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsi[k, index] = (rationalTSplineSecondDerivativesKsi[k, index] / sumKsiHeta -
																	  2 * rationalTSplineDerivativesKsi[k, index] * sumdKsiHeta /
																	  Math.Pow(sumKsiHeta, 2) -
																	  rationalTSplines[k, index] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
																	  2 * rationalTSplines[k, index] * Math.Pow(sumdKsiHeta, 2) /
																	  Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueHeta[k, index] = (rationalTSplineSecondDerivativesHeta[k, index] / sumKsiHeta -
																	   2 * rationalTSplineDerivativesHeta[k, index] * sumKsidHeta /
																	   Math.Pow(sumKsiHeta, 2) -
																	   rationalTSplines[k, index] * sumdHetadHeta /
																	   Math.Pow(sumKsiHeta, 2) +
																	   2 * rationalTSplines[k, index] * Math.Pow(sumKsidHeta, 2) /
																	   Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsiHeta[k, index] = (rationalTSplineSecondDerivativesKsiHeta[k, index] / sumKsiHeta -
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
		/// Two-dimensional T-spline shape functions from Bezier extraction for <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="TSplineKirchhoffLoveShellElementMaterial"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		public ShapeTSplines2DFromBezierExtraction(TSplineKirchhoffLoveShellElementMaterial element, ControlPoint[] controlPoints)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta,
				new List<Knot>
				{
					new Knot(){ID=0,Ksi=-1,Heta = -1,Zeta = 0},
					new Knot(){ID=1,Ksi=-1,Heta = 1,Zeta = 0},
					new Knot(){ID=2,Ksi=1,Heta = -1,Zeta = 0},
					new Knot(){ID=3,Ksi=1,Heta = 1,Zeta = 0}
				});

			var parametricGaussPointKsi = Vector.CreateZero(element.DegreeKsi + 1);
			for (int i = 0; i < element.DegreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i * (element.DegreeHeta + 1)].Ksi;
			}

			var parametricGaussPointHeta = Vector.CreateZero(element.DegreeHeta + 1);
			for (int i = 0; i < element.DegreeHeta + 1; i++)
			{
				parametricGaussPointHeta[i] = gaussPoints[i].Heta;
			}

			Vector knotValueVectorKsi = Vector.CreateZero((element.DegreeKsi + 1) * 2);
			Vector knotValueVectorHeta = Vector.CreateZero((element.DegreeHeta + 1) * 2);
			for (int i = 0; i < element.DegreeKsi + 1; i++)
			{
				knotValueVectorKsi[i] = -1;
				knotValueVectorKsi[element.DegreeKsi + 1 + i] = 1;
			}
			for (int i = 0; i < element.DegreeHeta + 1; i++)
			{
				knotValueVectorHeta[i] = -1;
				knotValueVectorHeta[element.DegreeHeta + 1 + i] = 1;
			}

			BSPLines1D bernsteinKsi = new BSPLines1D(element.DegreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bernsteinHeta = new BSPLines1D(element.DegreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
			bernsteinKsi.calculateBSPLinesAndDerivatives();
			bernsteinHeta.calculateBSPLinesAndDerivatives();

			int supportKsi = element.DegreeKsi + 1;
			int supportHeta = element.DegreeHeta + 1;

			var bKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);

			Matrix rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

			TSplineValues = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineDerivativeValuesKsi = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineDerivativeValuesHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsi = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsiHeta = Matrix.CreateZero(controlPoints.Length, supportKsi * supportHeta);

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
						TSplineValues[k, index] = rationalTSplines[k, index] * controlPoints[k].WeightFactor / sumKsiHeta;
						TSplineDerivativeValuesKsi[k, index] = (rationalTSplineDerivativesKsi[k, index] * sumKsiHeta -
																rationalTSplines[k, index] * sumdKsiHeta) /
															   Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						TSplineDerivativeValuesHeta[k, index] = (rationalTSplineDerivativesHeta[k, index] * sumKsiHeta -
																 rationalTSplines[k, index] * sumKsidHeta) /
																Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsi[k, index] = (rationalTSplineSecondDerivativesKsi[k, index] / sumKsiHeta -
																	  2 * rationalTSplineDerivativesKsi[k, index] * sumdKsiHeta /
																	  Math.Pow(sumKsiHeta, 2) -
																	  rationalTSplines[k, index] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
																	  2 * rationalTSplines[k, index] * Math.Pow(sumdKsiHeta, 2) /
																	  Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueHeta[k, index] = (rationalTSplineSecondDerivativesHeta[k, index] / sumKsiHeta -
																	   2 * rationalTSplineDerivativesHeta[k, index] * sumKsidHeta /
																	   Math.Pow(sumKsiHeta, 2) -
																	   rationalTSplines[k, index] * sumdHetadHeta /
																	   Math.Pow(sumKsiHeta, 2) +
																	   2 * rationalTSplines[k, index] * Math.Pow(sumKsidHeta, 2) /
																	   Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsiHeta[k, index] = (rationalTSplineSecondDerivativesKsiHeta[k, index] / sumKsiHeta -
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
		public Matrix TSplineDerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function derivatives per axis Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public Matrix TSplineDerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function second derivatives per axis Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public Matrix TSplineSecondDerivativesValueHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function second derivatives per axis Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public Matrix TSplineSecondDerivativesValueKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function mixed second derivatives per axis Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public Matrix TSplineSecondDerivativesValueKsiHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public Matrix TSplineValues { get; private set; }

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
