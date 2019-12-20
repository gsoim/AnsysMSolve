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
    public static class LinearRvesCompositeMaterialTet
    {
        [Fact]
        public static void CheckShellScaleTransitionsAndMicrostructureDense()
        {
            

            double E_disp = 3.5e9; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var outterMaterial = new ElasticMaterial3DtotalStrain() { PoissonRatio = 0.4, YoungModulus = 3.5e9 }; //pascal
            var innerMaterial = new ElasticMaterial3DtotalStrain() { PoissonRatio = 0.4, YoungModulus = 3.5e9 }; //pascal

            var Vec1 = Vector.CreateFromArray(new double[3] { 1, 0, 0 });
            var Vec2 = Vector.CreateFromArray(new double[3] { 0.5, 2, 0 });
            var strain = new double[3] { 0.01, 0, 0 };

            
            var material3 = new ShellElasticMaterial2Dtransformationb() { YoungModulus = E_disp, PoissonRatio = ni_disp, TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } };
            var Matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }
            material3.UpdateMaterial(strain);
            double[] stressesCheck3 = new double[3] { material3.Stresses[0], material3.Stresses[1], material3.Stresses[2] };
            double[,] consCheck1 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }


            IdegenerateRVEbuilder homogeneousRveBuilder1 = new CompositeMaterialModeluilderTet2(outterMaterial,innerMaterial,100,100,100);
            var material4 = new MicrostructureShell2D(homogeneousRveBuilder1, 
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1)
            {
                TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] },
                TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] }
            };
            var Matrix2 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix2[i1, i2] = material4.ConstitutiveMatrix[i1, i2]; } }
            material4.UpdateMaterial(strain);
            double[] stressesCheck4 = new double[3] { material4.Stresses[0], material4.Stresses[1], material4.Stresses[2] };
            double[,] consCheck2 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck2[i1, i2] = material4.ConstitutiveMatrix[i1, i2]; } }

            //-------------Check 2 steps savestate etc---------------
            material4.SaveState();
            material4.UpdateMaterial(new double[3] { 2 * strain[0], 2 * strain[1], 2 * strain[2] });
            double[] stressesCheck5 = new double[3] { material4.Stresses[0], material4.Stresses[1], material4.Stresses[2] };


            Assert.True(AreDisplacementsSame(stressesCheck3, stressesCheck4));
            Assert.True(AreDisplacementsSame(new double[3] { 2 * stressesCheck3[0], 2 * stressesCheck3[1], 2 * stressesCheck3[2] },
                                                                            stressesCheck5));
            Assert.True(BondSlipTest.AreDisplacementsSame(Matrix1.CopyToArray2D(), consCheck1));
            Assert.True(BondSlipTest.AreDisplacementsSame(Matrix1.CopyToArray2D(), material4.ConstitutiveMatrix.CopytoArray2D()));
        }

        [Fact]
        public static void CheckShellScaleTransitionsAndMicrostructure()
        {


            double E_disp = 3.5e9; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var outterMaterial = new ElasticMaterial3DtotalStrain() { PoissonRatio = 0.4, YoungModulus = 3.5e9 }; //pascal
            var innerMaterial = new ElasticMaterial3DtotalStrain() { PoissonRatio = 0.4, YoungModulus = 3.5e9 }; //pascal

            var Vec1 = Vector.CreateFromArray(new double[3] { 1, 0, 0 });
            var Vec2 = Vector.CreateFromArray(new double[3] { 0.5, 2, 0 });
            var strain = new double[3] { 0.01, 0, 0 };


            var material3 = new ShellElasticMaterial2Dtransformationb() { YoungModulus = E_disp, PoissonRatio = ni_disp, TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } };
            var Matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }
            material3.UpdateMaterial(strain);
            double[] stressesCheck3 = new double[3] { material3.Stresses[0], material3.Stresses[1], material3.Stresses[2] };
            double[,] consCheck1 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }


            IdegenerateRVEbuilder homogeneousRveBuilder1 = new CompositeMaterialModeluilderTet(outterMaterial, innerMaterial, 100, 100, 100);
            var material4 = new MicrostructureShell2D(homogeneousRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1)
            {
                TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] },
                TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] }
            };
            var Matrix2 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix2[i1, i2] = material4.ConstitutiveMatrix[i1, i2]; } }
            material4.UpdateMaterial(strain);
            double[] stressesCheck4 = new double[3] { material4.Stresses[0], material4.Stresses[1], material4.Stresses[2] };
            double[,] consCheck2 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck2[i1, i2] = material4.ConstitutiveMatrix[i1, i2]; } }

            //-------------Check 2 steps savestate etc---------------
            material4.SaveState();
            material4.UpdateMaterial(new double[3] { 2 * strain[0], 2 * strain[1], 2 * strain[2] });
            double[] stressesCheck5 = new double[3] { material4.Stresses[0], material4.Stresses[1], material4.Stresses[2] };


            Assert.True(AreDisplacementsSame(stressesCheck3, stressesCheck4));
            Assert.True(AreDisplacementsSame(new double[3] { 2 * stressesCheck3[0], 2 * stressesCheck3[1], 2 * stressesCheck3[2] },
                                                                            stressesCheck5));
            Assert.True(BondSlipTest.AreDisplacementsSame(Matrix1.CopyToArray2D(), consCheck1));
            Assert.True(BondSlipTest.AreDisplacementsSame(Matrix1.CopyToArray2D(), material4.ConstitutiveMatrix.CopytoArray2D()));
        }

        [Fact]
        public static void CheckSharedShellMaterialClasses()
        {


            double E_disp = 3.5e9; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var outterMaterial = new ElasticMaterial3DtotalStrain() { PoissonRatio = 0.4, YoungModulus = 3.5e9 }; //pascal
            var innerMaterial = new ElasticMaterial3DtotalStrain() { PoissonRatio = 0.4, YoungModulus = 3.5e9 }; //pascal

            var Vec1 = Vector.CreateFromArray(new double[3] { 1, 0, 0 });
            var Vec2 = Vector.CreateFromArray(new double[3] { 0.5, 2, 0 });
            var strain = new double[3] { 0.01, 0, 0 };


            var material3 = new ShellElasticMaterial2Dtransformationb() { YoungModulus = E_disp, PoissonRatio = ni_disp, TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } };
            var Matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }
            material3.UpdateMaterial(strain);
            double[] stressesCheck3 = new double[3] { material3.Stresses[0], material3.Stresses[1], material3.Stresses[2] };
            double[,] consCheck1 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }


            IdegenerateRVEbuilder homogeneousRveBuilder1 = new CompositeMaterialModeluilderTet(outterMaterial, innerMaterial, 100, 100, 100);
            var material4 = new Shell2dRVEMaterialHostConst(1,1,1,homogeneousRveBuilder1,
                model => (new SkylineSolver.Builder()).BuildSolver(model))
            {
                TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] },
                TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] }
            };
            //var Matrix2 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix2[i1, i2] = material4.ConstitutiveMatrix[i1, i2]; } }
            //material4.UpdateMaterial(strain);
            //double[] stressesCheck4 = new double[3] { material4.Stresses[0], material4.Stresses[1], material4.Stresses[2] };

            var clonedMAterial = material4.Clone();
            clonedMAterial.TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] };
            clonedMAterial.TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] };
            double[,] consCheck2 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck2[i1, i2] = clonedMAterial.ConstitutiveMatrix[i1, i2]; } }

            var clonedMAterialB = material4.Clone();
            clonedMAterialB.TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] };
            clonedMAterialB.TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] };
            var Matrix2 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix2[i1, i2] = clonedMAterialB.ConstitutiveMatrix[i1, i2]; } }

            var clonedMAterial3 = material4.Clone();
            clonedMAterial3.TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] };
            clonedMAterial3.TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] };
            var Matrixc = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrixc[i1, i2] = clonedMAterial3.ConstitutiveMatrix[i1, i2]; } }


            Assert.True(BondSlipTest.AreDisplacementsSame(Matrix2.CopyToArray2D(), consCheck2));
            Assert.True(BondSlipTest.AreDisplacementsSame(Matrixc.CopyToArray2D(), consCheck2));
            //Assert.True(BondSlipTest.AreDisplacementsSame(Matrix1.CopyToArray2D(), material4.ConstitutiveMatrix.CopytoArray2D()));
        }

        public static bool AreDisplacementsSame(double[] expectedValues,
            double[] computedValues)
        {
            var comparer = new ValueComparer(1E-12);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                if (!comparer.AreEqual(expectedValues[i1], computedValues[i1]))
                {
                    return false;
                }
            }
            return true;
        }


    }
}
