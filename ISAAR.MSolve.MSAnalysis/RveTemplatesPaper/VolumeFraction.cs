using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.MSAnalysis.RveTemplatesPaper
{
    public static class VolumeFraction
    {
        public static (double volumeFraction, double weightFraction) CalculatePercentage()
        {
            var a = 0.241;
            var cntDensity = 1.8;
            var matrixDensity = 1.4;
            var alpha = cntDensity / matrixDensity;
            var cntThickness = 0.34;

            //CNT Geometry
            var mi = 8.0;
            var ni = 8.0;
            var numberOfCNTs = 790;
            var cntLength = 98.0;

            var cntDiameter = (a / Math.PI) * Math.Sqrt(ni * ni + ni * mi + mi * mi);
            var cntRadius = cntDiameter / 2.0;
            var cntOuterRadius = cntRadius + (cntThickness / 2.0);
            var cntInnerRadius = cntRadius - (cntThickness / 2.0);


            // Matrix geometry
            var rveLength = 100.0;
            var rveHeight = 100.0;
            var rveWidth = 100.0;

            var outerCntVolume = Math.PI * (cntOuterRadius * cntOuterRadius) * cntLength;
            var innerCntVolume = Math.PI * (cntInnerRadius * cntInnerRadius) * cntLength;
            var cntVolume = outerCntVolume - innerCntVolume;

            var rveVolume = rveLength * rveWidth * rveHeight;
            var volumeFraction = ((numberOfCNTs * cntVolume) / (rveVolume - numberOfCNTs * outerCntVolume));
            var weightFraction = alpha * volumeFraction / (alpha * volumeFraction + 1.0);

            return (volumeFraction, weightFraction);
        }
    }
}
