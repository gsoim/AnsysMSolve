using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests.FEMpartB
{
    public class LinearCntRveTest
    {
        [Fact]
        public static void RunCntRveModel()
        {
            // Rve matrix material parameters
            var solidMaterial = new ElasticMaterial3D()
            {
                YoungModulus = 4.00, //GPa
                PoissonRatio = 0.40,
            };

            // cnt paramaters
            // define mechanical properties
            double youngModulus = 1.0; // 5490; // 
            double shearModulus = 1.0; // 871; // 
            double poissonRatio = (youngModulus / (2 * shearModulus)) - 1; //2.15; // 0.034;
            double area = 1776.65;  // CNT(20,20)-LinearEBE-TBT-L = 10nm
            double inertiaY = 1058.55;
            double inertiaZ = 1058.55;
            double torsionalInertia = 496.38;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;
            // CNT beams material and properties
            var beamMaterial = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };
            var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

        }
    }
}
