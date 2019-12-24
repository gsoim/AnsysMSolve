using System;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
    public class BSPLines1D: IShapeFunction1D
    {
        public int Degree { get; private set; }

        public double[] KnotValueVector { get; private set; }

        public double[] ParametricCoordinates { get; private set; }

        public double[,] Values { get; private set; }

        public double[,] DerivativeValues { get; private set; }

		public double [,] SecondDerivativeValues { get; private set; }

        public BSPLines1D(int degree, double[] knotValueVector, double[] parametricCoordinates)
        {
            if (degree <= 0)
            {
                throw new ArgumentException();
            }else if(knotValueVector == null){ throw new ArgumentNullException(); }
            this.Degree = degree;
            this.KnotValueVector = knotValueVector;
            this.ParametricCoordinates = parametricCoordinates;
        }

        public void calculateBSPLinesAndDerivatives()
        {
            int order = Degree + 1;
            int numberOfControlPoints = KnotValueVector.Length - order;
            int numberOfGaussPoints = ParametricCoordinates.Length;
            Values = new double[numberOfControlPoints + Degree, numberOfGaussPoints];
            DerivativeValues = new double[numberOfControlPoints + Degree, numberOfGaussPoints];
			SecondDerivativeValues = new double[numberOfControlPoints + Degree, numberOfGaussPoints];

            for (int i = 0; i < numberOfGaussPoints; i++)
                for (int j = 0; j < numberOfControlPoints+Degree; j++)
	                if (KnotValueVector[j]<=ParametricCoordinates[i]&& ParametricCoordinates[i] <= KnotValueVector[j + 1])
		                Values[j, i] = 1;
	                else
		                Values[j, i] = 0;

            for (int i = 1; i <= Degree; i++)
            {
                for (int j = 0; j < numberOfControlPoints+Degree-i; j++)
                {
                    for (int k = 0; k < numberOfGaussPoints; k++)
                    {
                        double additive1 = 0;
                        double additive2 = 0;
                        double additiveDerivative1 = 0;
                        double additiveDerivative2 = 0;
	                    double additiveSecondDerivative1 = 0;
	                    double additiveSecondDerivative2 = 0;
						double denominator1 = KnotValueVector[j + i] - KnotValueVector[j];
                        double denominator2 = KnotValueVector[j + i + 1] - KnotValueVector[j + 1];
                        if (denominator1 != 0)
                        {
                            additive1 = (ParametricCoordinates[k] - KnotValueVector[j]) 
                                / denominator1 * Values[j, k];
                            additiveDerivative1 = ((ParametricCoordinates[k] - KnotValueVector[j])
                                * DerivativeValues[j, k] + Values[j, k]) / denominator1;
	                        additiveSecondDerivative1 = Degree / denominator1 * DerivativeValues[j, k];
                        }
                        if (denominator2 != 0)
                        {
                            additive2 = (KnotValueVector[j + i + 1] - ParametricCoordinates[k])
                                / denominator2 * Values[j + 1, k];
                            additiveDerivative2 = ((KnotValueVector[j + i + 1] - ParametricCoordinates[k])
                                * DerivativeValues[j + 1, k] - Values[j + 1, k]) / denominator2;
	                        additiveSecondDerivative2 = -Degree / denominator2 * DerivativeValues[j + 1, k];
                        }
                        Values[j, k] = additive1 + additive2;
                        DerivativeValues[j, k] = additiveDerivative1 + additiveDerivative2;
	                    SecondDerivativeValues[j, k] = additiveSecondDerivative1 + additiveSecondDerivative2;
                    }
                }
            }
        }
    }
}
