using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Materials;
using MathNet.Numerics.Data.Matlab;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
    public class ShellSectionMaterialTests
    {
        private const double Tolerance = 1e-9;

        private ShellElasticSectionMaterial2D Material => new ShellElasticSectionMaterial2D()
        {
            YoungModulus = 1200000,
            PoissonRatio = 0.0,
            TangentVectorV1 = new double[] { 1.000000000000000000000000000000000000, 0.000000000000000053342746886286800000, 0.000000000000000000000000000000000000 },
            TangentVectorV2 = new double[] { 3.90312782094781000000000000000000E-18, 9.99999999999999000000000000000000E-01, 0.00000000000000000000000000000000E+00, },
            NormalVectorV3 = new double[] { 0, 0, 1 },
            Thickness = 0.1
        };
        
        [Fact]
        public void ConstitutiveMatricesTests()
        {
            var MembraneConstitutiveMatrix = Material.MembraneConstitutiveMatrix;
            var BendingConstitutiveMatrix = Material.BendingConstitutiveMatrix;

            var expectedConstitutiveMembrane = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "MembraneConstitutiveMatrix");
            var expectedConstitutiveBending = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "BendingConstitutiveMatrix");

            for (int i = 0; i < MembraneConstitutiveMatrix.NumRows; i++)
            {
                for (int j = 0; j < MembraneConstitutiveMatrix.NumColumns; j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedConstitutiveMembrane[i, j], MembraneConstitutiveMatrix[i, j], Tolerance));
                    Assert.True(Utilities.AreValuesEqual(expectedConstitutiveBending[i, j], BendingConstitutiveMatrix[i, j], Tolerance));
                }
            }
        }
        
        [Fact]
        public void StressesTests()
        {
            var membraneStrains = new double[]
            {
                -2.72674706425224000000000E-07,
                0.00000000000000000000000E+00,
                2.81820360153431000000000E-18
            };

            var bendingStrains = new double[]
            {
                0.00426183634012108000000000000000000000,
                0.00000000000000000000756698343313268000,
                -0.00000000000000013630925450307800000000,
            };
            var material = Material;
            material.UpdateMaterial(membraneStrains, bendingStrains);

            var membraneForces = material.MembraneForces;
            var bendingMoments= material.Moments;

            var expectedMembraneForces = new double[]
            {
                -0.03272096477102780000000000000000000000000000,
                -0.00000000000000000000000000001936098854493600,
                0.00000000000016909408923230800000000000000000
            };

            var expectedBendingMoments = new double[]
            {
                -0.4261836340121090000000000000000000000000,
                -0.0000000000000000007566983433140510000000,
                0.0000000000000068398599800688700000000000
            };

            for (int i = 0; i < 3; i++)
            {
                Assert.True(Utilities.AreValuesEqual(expectedMembraneForces[i], membraneForces[i], Tolerance));
                Assert.True(Utilities.AreValuesEqual(expectedBendingMoments[i], bendingMoments[i], Tolerance));
            }


        }
    }
}
