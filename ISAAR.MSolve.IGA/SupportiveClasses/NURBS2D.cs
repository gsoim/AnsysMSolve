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
	using MGroup.MSolve.Geometry.Coordinates;

	/// <summary>
	/// Two-dimensional NURBS shape functions.
	/// </summary>
	public class Nurbs2D
	{

		/// <summary>
		/// Defines a 2D NURBS shape functions for a collocation point.
		/// Calculates NURBS for only Control Points affected by the collocation point.
		/// </summary>
		/// <param name="degreeKsi">Polynomial degree of the parametric axis Ksi.</param>
		/// <param name="degreeHeta">Polynomial degree of the parametric axis Heta.</param>
		/// <param name="knotValueVectorKsi">Knot value vector of the parametric axis Ksi.</param>
		/// <param name="knotValueVectorHeta">Knot value vector of the parametric axis Heta.</param>
		/// <param name="collocationPoint">A <see cref="NaturalPoin"/> for which the shape functions will be evaluated.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		public Nurbs2D(int degreeKsi, int degreeHeta, Vector knotValueVectorKsi, Vector knotValueVectorHeta,
			NaturalPoint collocationPoint, IList<ControlPoint> controlPoints)
		{
			var numberOfControlPointsHeta = knotValueVectorHeta.Length - degreeHeta - 1;

			BSPLines1D bsplinesKsi = new BSPLines1D(degreeKsi, knotValueVectorKsi,
				Vector.CreateFromArray(new double[] { collocationPoint.Xi }));
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
				Vector.CreateFromArray(new double[] { collocationPoint.Eta }));
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();

			int numberOfGPKsi = 1;
			int numberOfGPHeta = 1;
			int numberOfElementControlPoints = (degreeKsi + 1) * (degreeHeta + 1);

			NurbsValues = new double[numberOfElementControlPoints, 1];
			NurbsDerivativeValuesKsi = new double[numberOfElementControlPoints, 1];
			NurbsDerivativeValuesHeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueKsi = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueHeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueKsiHeta = new double[numberOfElementControlPoints, 1];

			for (int i = 0; i < numberOfGPKsi; i++)
			{
				for (int j = 0; j < numberOfGPHeta; j++)
				{
					double sumKsiHeta = 0;
					double sumdKsiHeta = 0;
					double sumKsidHeta = 0;
					double sumdKsidKsi = 0;
					double sumdHetadHeta = 0;
					double sumdKsidHeta = 0;

					for (int k = 0; k < numberOfElementControlPoints; k++)
					{
						// Why type casting is needed.?

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

						NurbsValues[k, i * numberOfGPHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] *
							bsplinesHeta.BSPLineValues[indexHeta, j] *
							controlPoints[k].WeightFactor / sumKsiHeta;

						NurbsDerivativeValuesKsi[k, i * numberOfGPHeta + j] =
							bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHeta -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHeta) / Math.Pow(sumKsiHeta, 2);

						NurbsDerivativeValuesHeta[k, i * numberOfGPHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHeta -
							 bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHeta) / Math.Pow(sumKsiHeta, 2);

						NurbsSecondDerivativeValueKsi[k, i * numberOfGPHeta + j] =
							bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHeta -
							 2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHeta /
							 Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) /
							 Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueHeta[k, i * numberOfGPHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHeta -
							 2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHeta /
							 Math.Pow(sumKsiHeta, 2) -
							 bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHeta, 2) /
							 Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueKsiHeta[k, i * numberOfGPHeta + j] =
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
		/// Defines a 2D NURBS shape functions for a collocation point.
		/// Calculates NURBS for all Control Points.
		/// </summary>
		/// <param name="degreeKsi">Polynomial degree of the parametric axis Ksi.</param>
		/// <param name="degreeHeta">Polynomial degree of the parametric axis Heta.</param>
		/// <param name="knotValueVectorKsi">Knot value vector of the parametric axis Ksi.</param>
		/// <param name="knotValueVectorHeta">Knot value vector of the parametric axis Heta.</param>
		/// <param name="collocationPoint">A <see cref="NaturalPoin"/> for which the shape functions will be evaluated.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="calculateAllFunctions"></param>
		public Nurbs2D(int degreeKsi, int degreeHeta, Vector knotValueVectorKsi, Vector knotValueVectorHeta,
			NaturalPoint collocationPoint, IList<ControlPoint> controlPoints, bool calculateAllFunctions)
		{
			var numberOfControlPointsHeta = knotValueVectorHeta.Length - degreeHeta - 1;

			BSPLines1D bsplinesKsi = new BSPLines1D(degreeKsi, knotValueVectorKsi,
				Vector.CreateFromArray(new double[] { collocationPoint.Xi }));
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
				Vector.CreateFromArray(new double[] { collocationPoint.Eta }));
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();

			int numberOfGPKsi = 1;
			int numberOfGPHeta = 1;
			int numberOfElementControlPoints = controlPoints.Count;

			NurbsValues = new double[numberOfElementControlPoints, 1];
			NurbsDerivativeValuesKsi = new double[numberOfElementControlPoints, 1];
			NurbsDerivativeValuesHeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueKsi = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueHeta = new double[numberOfElementControlPoints, 1];
			NurbsSecondDerivativeValueKsiHeta = new double[numberOfElementControlPoints, 1];

			for (int i = 0; i < numberOfGPKsi; i++)
			{
				for (int j = 0; j < numberOfGPHeta; j++)
				{
					double sumKsiHeta = 0;
					double sumdKsiHeta = 0;
					double sumKsidHeta = 0;
					double sumdKsidKsi = 0;
					double sumdHetadHeta = 0;
					double sumdKsidHeta = 0;

					for (int k = 0; k < numberOfElementControlPoints; k++)
					{
						// Why type casting is needed.?

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

						NurbsValues[k, i * numberOfGPHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] *
							bsplinesHeta.BSPLineValues[indexHeta, j] *
							controlPoints[k].WeightFactor / sumKsiHeta;

						NurbsDerivativeValuesKsi[k, i * numberOfGPHeta + j] =
							bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHeta -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHeta) / Math.Pow(sumKsiHeta, 2);

						NurbsDerivativeValuesHeta[k, i * numberOfGPHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHeta -
							 bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHeta) / Math.Pow(sumKsiHeta, 2);

						NurbsSecondDerivativeValueKsi[k, i * numberOfGPHeta + j] =
							bsplinesHeta.BSPLineValues[indexHeta, j] * controlPoints[k].WeightFactor *
							(bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHeta -
							 2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHeta /
							 Math.Pow(sumKsiHeta, 2) -
							 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) /
							 Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueHeta[k, i * numberOfGPHeta + j] =
							bsplinesKsi.BSPLineValues[indexKsi, i] * controlPoints[k].WeightFactor *
							(bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHeta -
							 2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHeta /
							 Math.Pow(sumKsiHeta, 2) -
							 bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
							 2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHeta, 2) /
							 Math.Pow(sumKsiHeta, 3));

						NurbsSecondDerivativeValueKsiHeta[k, i * numberOfGPHeta + j] =
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
		/// Defines a 2D NURBS shape function for an element given the control points.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		public Nurbs2D(Element element, ControlPoint[] controlPoints)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi,
				element.Patch.DegreeHeta, element.Knots.ToArray());

			var parametricGaussPointKsi = Vector.CreateZero(element.Patch.DegreeKsi + 1);
			for (int i = 0; i < element.Patch.DegreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i * (element.Patch.DegreeHeta + 1)].Ksi;
			}

			var parametricGaussPointHeta = Vector.CreateZero(element.Patch.DegreeHeta + 1);
			for (int i = 0; i < element.Patch.DegreeHeta + 1; i++)
			{
				parametricGaussPointHeta[i] = gaussPoints[i].Heta;
			}

			BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi,
				parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(element.Patch.DegreeHeta, element.Patch.KnotValueVectorHeta,
				parametricGaussPointHeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();

			int supportKsi = element.Patch.DegreeKsi + 1;
			int supportHeta = element.Patch.DegreeHeta + 1;
			int numberOfElementControlPoints = supportKsi * supportHeta;

			NurbsValues = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsDerivativeValuesKsi = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsDerivativeValuesHeta = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsSecondDerivativeValueKsi = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsSecondDerivativeValueHeta = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsSecondDerivativeValueKsiHeta = new double[numberOfElementControlPoints, gaussPoints.Count];

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
						int indexKsi = controlPoints[k].ID / element.Patch.NumberOfControlPointsHeta;
						int indexHeta = controlPoints[k].ID % element.Patch.NumberOfControlPointsHeta;
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
						int indexKsi = controlPoints[k].ID / element.Patch.NumberOfControlPointsHeta;
						int indexHeta = controlPoints[k].ID % element.Patch.NumberOfControlPointsHeta;

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
		/// Defines a 2D NURBS shape function for an element given the per axis gauss point coordinates.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="parametricGaussPointKsi">An <see cref="IVector"/> containing Gauss points of axis Ksi.</param>
		/// <param name="parametricGaussPointHeta">An <see cref="IVector"/> containing Gauss points of axis Heta.</param>
		public Nurbs2D(Element element, ControlPoint[] controlPoints, IVector parametricGaussPointKsi,
			IVector parametricGaussPointHeta)
		{
			var parametricPointsCount = parametricGaussPointKsi.Length * parametricGaussPointHeta.Length;

			BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi,
				parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(element.Patch.DegreeHeta, element.Patch.KnotValueVectorHeta,
				parametricGaussPointHeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();

			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;
			int numberOfElementControlPoints = (element.Patch.DegreeKsi + 1) * (element.Patch.DegreeHeta + 1);

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
						int indexKsi = controlPoints[k].ID / element.Patch.NumberOfControlPointsHeta;
						int indexHeta = controlPoints[k].ID % element.Patch.NumberOfControlPointsHeta;
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
						int indexKsi = controlPoints[k].ID / element.Patch.NumberOfControlPointsHeta;
						int indexHeta = controlPoints[k].ID % element.Patch.NumberOfControlPointsHeta;

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
		/// Defines n 2D NURBS shape function for a face element.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement2D"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="face">The two-dimensional boundary entities where the <paramref name="element"/> shape functions will be evaluated.</param>
		public Nurbs2D(Element element, ControlPoint[] controlPoints, Face face)
		{
			var degreeKsi = face.Degrees[0];
			var degreeHeta = face.Degrees[1];
			var knotValueVectorKsi = face.KnotValueVectors[0];
			var knotValueVectorHeta = face.KnotValueVectors[1];
			var numberOfControlPointsHeta = knotValueVectorHeta.Length - degreeHeta - 1;

			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints =
				gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, element.Knots.ToArray());

			var parametricGaussPointKsi = Vector.CreateZero(degreeKsi + 1);
			for (int i = 0; i < degreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i * (degreeHeta + 1)].Ksi;
			}

			var parametricGaussPointHeta = Vector.CreateZero(degreeHeta + 1);
			for (int i = 0; i < degreeHeta + 1; i++)
			{
				parametricGaussPointHeta[i] = gaussPoints[i].Heta;
			}

			BSPLines1D bsplinesKsi = new BSPLines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
			bsplinesKsi.calculateBSPLinesAndDerivatives();
			bsplinesHeta.calculateBSPLinesAndDerivatives();

			int supportKsi = degreeKsi + 1;
			int supportHeta = degreeHeta + 1;
			int numberOfElementControlPoints = supportKsi * supportHeta;

			NurbsValues = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsDerivativeValuesKsi = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsDerivativeValuesHeta = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsSecondDerivativeValueKsi = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsSecondDerivativeValueHeta = new double[numberOfElementControlPoints, gaussPoints.Count];
			NurbsSecondDerivativeValueKsiHeta = new double[numberOfElementControlPoints, gaussPoints.Count];

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
						int indexKsi =
							face.ControlPointsDictionary.First(cp => cp.Value.ID == controlPoints[k].ID).Key /
							numberOfControlPointsHeta;
						int indexHeta =
							face.ControlPointsDictionary.First(cp => cp.Value.ID == controlPoints[k].ID).Key %
							numberOfControlPointsHeta;
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
						int indexKsi =
							face.ControlPointsDictionary.First(cp => cp.Value.ID == controlPoints[k].ID).Key /
							numberOfControlPointsHeta;
						int indexHeta =
							face.ControlPointsDictionary.First(cp => cp.Value.ID == controlPoints[k].ID).Key %
							numberOfControlPointsHeta;

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
