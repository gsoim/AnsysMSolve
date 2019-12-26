
using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Solvers;


namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class Shell2dRVESectionMaterialHostConst : IShellSectionMaterial
    {
        public double[] NormalVectorV3 { get; set; }
        public double[] TangentVectorV1 { get; set; }
        public double[] TangentVectorV2 { get; set; }
        private Matrix consMat;

        private int rveDatabaseSize;
        private int gpsPerRve;
        private int gpCounter;
        private IdegenerateRVEbuilder rveBuilderToClone;

        //IShellMaterial coreMaterial;
        double[,] coreConstitutive;

        private Dictionary<int, IShellMaterial> coreMaterials;
        private Dictionary<int, double[,]> MaterialDatabase;

        private Func<Model, ISolver> rveSolver = model => (new Solvers.Direct.SkylineSolver.Builder()).BuildSolver(model);

        public Shell2dRVESectionMaterialHostConst(int rveDatabaseSize, int gpsPerRve, int gpCounter,
            IdegenerateRVEbuilder rVEbuilder,Func<Model,ISolver> rveSolver, double thickness)
        {
            this.rveDatabaseSize = rveDatabaseSize;
            this.gpsPerRve = gpsPerRve;
            this.gpCounter = gpCounter;
            coreMaterials = new Dictionary<int, IShellMaterial>();
            this.rveBuilderToClone = rVEbuilder;
            CreateConstitutiveTensorsDatabase();
            this.rveSolver = rveSolver;
            Thickness = thickness;
        }

        private void CreateConstitutiveTensorsDatabase()
        {
            MaterialDatabase = new Dictionary<int, double[,]>(rveDatabaseSize);

            for (int rve_id=1; rve_id<rveDatabaseSize+1;rve_id++)
            {
                var newMaterial = new
                MicrostructureShell2D((IdegenerateRVEbuilder)rveBuilderToClone.Clone(rve_id),rveSolver, true, 1);
                newMaterial.TangentVectorV1 = new double[3] { 1, 0, 0 };
                newMaterial.TangentVectorV2 = new double[3] { 0, 1, 0 };
                newMaterial.NormalVectorV3 = new double[3] { 0, 0, 1 };
                var newCOnstitutive = new double[3, 3];
                for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { newCOnstitutive[i1, i2] = newMaterial.ConstitutiveMatrix[i1, i2]; } }
                MaterialDatabase.Add(rve_id, newCOnstitutive);
            }

        }

        //private Shell2dRVEMaterialHostConst(IShellMaterial coreMaterial)
        //{
        //    this.coreMaterial = coreMaterial;
        //}

        private Shell2dRVESectionMaterialHostConst(double[,] coreCons)
        {
            this.coreConstitutive = coreCons;
        }

        public IShellSectionMaterial Clone()
        {
            gpCounter += 1;

            //if ((gpCounter/ (rveDatabaseSize*gpsPerRve)==1)&&(gpCounter% (rveDatabaseSize * gpsPerRve)==1))
            //{
            //    gpCounter += -(rveDatabaseSize * gpsPerRve);
            //}

            int rve_id = 1;// ((gpCounter-1) / gpsPerRve) + 1;      //TODO: xrhsh tou rand()     

            var sharedRVEmaterial = new Shell2dRVESectionMaterialHostConst(MaterialDatabase[rve_id]);
            return sharedRVEmaterial;
        }

        object ICloneable.Clone()
        {
            return this.Clone();
        }

        public IMatrixView ConstitutiveMatrix
        {
            get
            {
                if (consMat == null) BuildConstitutive();
                return consMat;
            }
        }

        private void BuildConstitutive()
        {
            var transformationMatrix=this.CalculateTransformationMatrix(Vector.CreateFromArray(TangentVectorV1),Vector.CreateFromArray(TangentVectorV2));
            var coreConstCopy = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { coreConstCopy[i1, i2] = coreConstitutive[i1, i2]; } }
            consMat = (transformationMatrix.Transpose() * ( Matrix.CreateFromArray(coreConstCopy)) * transformationMatrix);
        }

        private Matrix CalculateTransformationMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
        {
            var auxMatrix1 = Matrix.CreateZero(2, 2);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            (Matrix inverse, double det) = auxMatrix1.InvertAndDeterminant(); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)

            //Contravariant base vectors
            double[][] G_i = new double[2][];
            for (int i1 = 0; i1 < 2; i1++)
            {
                G_i[i1] = new double[3];
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i1][i2] = inverse[i1, 0] * surfaceBasisVector1[i2] + inverse[i1, 1] * surfaceBasisVector2[i2];
                }
            }

            //Normalised covariant base vectors
            double[][] Ei = new double[2][];// to trito den xreiazetai

            //Ei[0] = new double[3]; Array.Copy(surfaceBasisVector1.Data, Ei[0], 3);
            Ei[0] = surfaceBasisVector1.CopyToArray();
            double G1_norm = surfaceBasisVector1.Norm2();
            for (int i1 = 0; i1 < 3; i1++) { Ei[0][i1] = Ei[0][i1] / G1_norm; }

            double G2_dot_E1 = 0;
            for (int i1 = 0; i1 < 3; i1++) { G2_dot_E1 += surfaceBasisVector2[i1] * Ei[0][i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = G2_dot_E1 * Ei[0][i1]; }

            Ei[1] = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = surfaceBasisVector2[i1] - projection[i1]; }
            double norm1 = Vector.CreateFromArray(Ei[1]).Norm2();
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = Ei[1][i1] / norm1; }

            double[,] EiDOTG_j = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    EiDOTG_j[i1, i2] = Vector.CreateFromArray(Ei[i1]).DotProduct(Vector.CreateFromArray(G_i[i2]));
                }
            }

            var transformationMatrix = Matrix.CreateFromArray(new double[3, 3] { {EiDOTG_j[0,0]*EiDOTG_j[0,0],EiDOTG_j[0,1]*EiDOTG_j[0,1],EiDOTG_j[0,0]*EiDOTG_j[0,1]  },
                 {EiDOTG_j[1,0]*EiDOTG_j[1,0],EiDOTG_j[1,1]*EiDOTG_j[1,1],EiDOTG_j[1,0]*EiDOTG_j[1,1]  },
                {2*EiDOTG_j[1,0]*EiDOTG_j[0,0],2*EiDOTG_j[1,1]*EiDOTG_j[0,1],EiDOTG_j[1,0]*EiDOTG_j[0,1]+EiDOTG_j[1,1]*EiDOTG_j[0,0]   } });

            return transformationMatrix;
        }

        private double[] trueStressVec;
        public double[] Stresses => throw new NotImplementedException();

        public int ID => throw new NotImplementedException();

        public bool Modified => throw new NotImplementedException();

        public double[] Coordinates { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public double YoungModulus => throw new NotSupportedException();

        public double PoissonRatio => throw new NotSupportedException();

        public double[] MembraneForces => throw new NotImplementedException();

        public double[] Moments => throw new NotImplementedException();

        public double Thickness { get; set; }

        public Matrix MembraneConstitutiveMatrix =>
            ConstitutiveMatrix.CopyToFullMatrix();

        public Matrix BendingConstitutiveMatrix =>
            ConstitutiveMatrix.CopyToFullMatrix();

        public Matrix CouplingConstitutiveMatrix => Matrix.CreateZero(3, 3);

        public void UpdateMaterial(double[] GLvec)
        {
            throw new NotSupportedException();
        }

        public void ResetModified()
        {
            throw new NotSupportedException();
        }

        public void SaveState()
        {
            throw new NotSupportedException();
        }

        public void ClearState()
        {
            throw new NotSupportedException();
        }

        public void ClearStresses()
        {
            throw new NotSupportedException();
        }
        
        public void UpdateMaterial(double[] membraneStrains, double[] bendingStrains)
        {
            throw new NotImplementedException();
        }
    }
}
