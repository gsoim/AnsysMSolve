using System.Collections.Generic;
using ISAAR.MSolve.IGA.Entities;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
    /// <summary>
	/// Gauss Quadrature rules.
	/// </summary>
	public class GaussQuadrature
	{
		#region Coordinates and Weights

		//https://pomax.github.io/bezierinfo/legendre-gauss.html
		private static readonly double[][] coordinate = new double[][]
		{
			new double[]{ 0.0 },
			new double[]{ -0.5773502691896257, 0.5773502691896257 },
			new double[]{ -0.7745966692414836, 0.000000000000000, 0.7745966692414836 },
			new double[]{ -0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526 },
			new double[]{ -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640 },
			new double[]{ -0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645,0.9324695142031521 },
			new double[]{ -0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972,0.7415311855993945, 0.9491079123427585 },
			new double[]{ -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 }
		};

		private static readonly double[][] weight = new double[][]
		{
			new double[]{ 2.0 },
			new double[]{ 1.0, 1.0 },
			new double[]{ 0.555555555555556, 0.888888888888889, 0.555555555555556 },
			new double[]{ 0.34785484513745385, 0.6521451548625461, 0.6521451548625461, 0.34785484513745385 },
			new double[]{ 0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 },
			new double[]{ 0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704 },
			new double[]{ 0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697 },
			new double[]{ 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 }
		};

		#endregion Coordinates and Weights

		/// <summary>
		/// Create full gauss quadrature rule for a two-dimensional isogeometric element.
		/// </summary>
		/// <param name="degreeKsi">Parametric degree Ksi.</param>
		/// <param name="degreeHeta">Parametric degree Heta.</param>
		/// <param name="knotsOfElement">Knots of that will define the integration boundaries.</param>
		/// <returns>An <see cref="IList{T}"/> with the genertaed <see cref="GaussLegendrePoint3D"/>.</returns>
		public IList<GaussLegendrePoint3D> CalculateElementGaussPoints(int degreeKsi, int degreeHeta, IList<Knot> knotsOfElement)
		{
			IList<GaussLegendrePoint3D> gaussPointsofElement = new List<GaussLegendrePoint3D>();
			int numberOfGPKsi = degreeKsi + 1;
			int numberOfGPHeta = degreeHeta + 1;
			double[] coordinatesKsi = new double[numberOfGPKsi];
			double[] weightKsi = new double[numberOfGPKsi];
			double[] coordinatesHeta = new double[numberOfGPHeta];
			double[] weightHeta = new double[numberOfGPHeta];

			for (int i = 0; i < numberOfGPKsi; i++)
			{
				coordinatesKsi[i] = 0.5 * (knotsOfElement[0].Ksi + knotsOfElement[2].Ksi
					+ (knotsOfElement[2].Ksi - knotsOfElement[0].Ksi) * coordinate[degreeKsi][i]);
				weightKsi[i] = 0.5 * ((knotsOfElement[2].Ksi - knotsOfElement[0].Ksi) * weight[degreeKsi][i]);
			}

			for (int i = 0; i < numberOfGPHeta; i++)
			{
				coordinatesHeta[i] = 0.5 * (knotsOfElement[0].Heta + knotsOfElement[1].Heta
					+ (knotsOfElement[1].Heta - knotsOfElement[0].Heta) * coordinate[degreeHeta][i]);
				weightHeta[i] = 0.5 * ((knotsOfElement[1].Heta - knotsOfElement[0].Heta) * weight[degreeHeta][i]);
			}

			for (int i = 0; i < numberOfGPKsi; i++)
			{
				for (int j = 0; j < numberOfGPHeta; j++)
				{
					gaussPointsofElement.Add(new GaussLegendrePoint3D(coordinatesKsi[i], coordinatesHeta[j], 0, weightKsi[i] * weightHeta[j]));
				}
			}
			return gaussPointsofElement;
		}

		internal IList<GaussLegendrePoint3D> CalculateElementGaussPoints(int degreeKsi, IList<Knot> knotsOfElement)
		{
			IList<GaussLegendrePoint3D> gaussPointsofElement = new List<GaussLegendrePoint3D>();
			int numberOfGPKsi = degreeKsi + 1;
			double[] coordinatesKsi = new double[numberOfGPKsi];
			double[] weightKsi = new double[numberOfGPKsi];

			for (int i = 0; i < numberOfGPKsi; i++)
			{
				coordinatesKsi[i] = 0.5 * (knotsOfElement[0].Ksi + knotsOfElement[1].Ksi
					+ (knotsOfElement[1].Ksi - knotsOfElement[0].Ksi) * coordinate[degreeKsi][i]);
				weightKsi[i] = 0.5 * ((knotsOfElement[1].Ksi - knotsOfElement[0].Ksi) * weight[degreeKsi][i]);
			}

			for (int i = 0; i < numberOfGPKsi; i++)
			{
				gaussPointsofElement.Add(new GaussLegendrePoint3D(coordinatesKsi[i], 0, 0, weightKsi[i]));
			}
			return gaussPointsofElement;
		}

		/// <summary>
		/// Create full gauss quadrature rule for a three-dimensional isogeometric element.
		/// </summary>
		/// <param name="degreeKsi">Parametric degree Ksi.</param>
		/// <param name="degreeHeta">Parametric degree Heta.</param>
		/// /// <param name="degreeZeta">Parametric degree Zeta.</param>
		/// <param name="knotsOfElement">An <see cref="IList{T}"/> with the genertaed <see cref="GaussLegendrePoint3D"/>.</param>
		/// <returns>An <see cref="IList{T}"/> with the genertaed <see cref="GaussLegendrePoint3D"/>.</returns>
		public GaussLegendrePoint3D[] CalculateElementGaussPoints(int degreeKsi, int degreeHeta, int degreeZeta, IList<Knot> knotsOfElement)
		{
			int numberOfGPKsi = degreeKsi + 1;
			int numberOfGPHeta = degreeHeta + 1;
			int numberOfGPZeta = degreeZeta + 1;
			var gaussPointsofElement = new GaussLegendrePoint3D[numberOfGPKsi * numberOfGPHeta * numberOfGPZeta];
			double[] coordinatesKsi = new double[numberOfGPKsi];
			double[] weightKsi = new double[numberOfGPKsi];
			double[] coordinatesHeta = new double[numberOfGPHeta];
			double[] weightHeta = new double[numberOfGPHeta];
			double[] coordinatesZeta = new double[numberOfGPZeta];
			double[] weightZeta = new double[numberOfGPZeta];

			for (int i = 0; i < numberOfGPKsi; i++)
			{
				coordinatesKsi[i] = 0.5 * (knotsOfElement[0].Ksi + knotsOfElement[4].Ksi
					+ (knotsOfElement[4].Ksi - knotsOfElement[0].Ksi) * coordinate[degreeKsi][i]);
				weightKsi[i] = 0.5 * ((knotsOfElement[4].Ksi - knotsOfElement[0].Ksi) * weight[degreeKsi][i]);
			}

			for (int i = 0; i < numberOfGPHeta; i++)
			{
				coordinatesHeta[i] = 0.5 * (knotsOfElement[0].Heta + knotsOfElement[2].Heta
					+ (knotsOfElement[2].Heta - knotsOfElement[0].Heta) * coordinate[degreeHeta][i]);
				weightHeta[i] = 0.5 * ((knotsOfElement[2].Heta - knotsOfElement[0].Heta) * weight[degreeHeta][i]);
			}

			for (int i = 0; i < numberOfGPZeta; i++)
			{
				coordinatesZeta[i] = 0.5 * (knotsOfElement[0].Zeta + knotsOfElement[1].Zeta
					+ (knotsOfElement[1].Zeta - knotsOfElement[0].Zeta) * coordinate[degreeZeta][i]);
				weightZeta[i] = 0.5 * ((knotsOfElement[1].Zeta - knotsOfElement[0].Zeta) * weight[degreeZeta][i]);
			}

			var index = 0;
			for (int i = 0; i < numberOfGPKsi; i++)
			{
				for (int j = 0; j < numberOfGPHeta; j++)
				{
					for (int k = 0; k < numberOfGPZeta; k++)
					{
						gaussPointsofElement[index++] = new GaussLegendrePoint3D(coordinatesKsi[i], coordinatesHeta[j], coordinatesZeta[k], weightKsi[i] * weightHeta[j] * weightZeta[k]);
					}
				}
			}
			return gaussPointsofElement;
		}
	}

	/// <summary>
	/// General three dimensional Gauss Point.
	/// </summary>
	public class GaussLegendrePoint3D
	{
		/// <summary>
		/// Parametric coordinate Ksi of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double Ksi { get; private set; }

		/// <summary>
		/// Parametric coordinate Heta of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double Heta { get; private set; }

		/// <summary>
		/// Parametric coordinate Zeta of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double Zeta { get; private set; }

		/// <summary>
		/// Weight factor of the <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		public double WeightFactor { get; private set; }

		/// <summary>
		/// Defines a <see cref="GaussLegendrePoint3D"/>.
		/// </summary>
		/// <param name="ksi">Parametric coordinate Ksi.</param>
		/// <param name="heta">Parametric coordinate Heta.</param>
		/// <param name="zeta">Parametric coordinate Zeta.</param>
		/// <param name="weightFactor">Weight factor.</param>
		public GaussLegendrePoint3D(double ksi, double heta, double zeta, double weightFactor)
		{
			this.Ksi = ksi;
			this.Heta = heta;
			this.Zeta = zeta;
			this.WeightFactor = weightFactor;
		}
	}
}
