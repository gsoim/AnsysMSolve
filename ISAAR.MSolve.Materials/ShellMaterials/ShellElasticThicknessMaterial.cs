using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.Materials.ShellMaterials
{
    public class ShellElasticThicknessMaterial:IShellSectionMaterial
    {

        public ShellElasticThicknessMaterial(int integrationOrder, double thickness,
            double youngModulus, double poissonRatio)
        {
            Thickness = thickness;
            _integrationOrder = integrationOrder;
            YoungModulus = youngModulus;
            PoissonRatio = poissonRatio;
            _thicknessIntegrationPoints = GaussLegendre1D.GetQuadratureOfSpan(_integrationOrder,
                -Thickness / 2.0, Thickness / 2.0);

            foreach (var gaussPoint in _thicknessIntegrationPoints.IntegrationPoints)
            {
                _thicknessIntegrationPointMaterials.Add(gaussPoint, new ShellElasticMaterial2Dtransformationb()
                {
                    YoungModulus = youngModulus,
                    PoissonRatio = PoissonRatio
                });
            }

            MembraneConstitutiveMatrix = Matrix.CreateZero(3, 3);
            BendingConstitutiveMatrix = Matrix.CreateZero(3, 3);
            CouplingConstitutiveMatrix = Matrix.CreateZero(3, 3);
        }

        private Dictionary<GaussPoint, IShellMaterial> _thicknessIntegrationPointMaterials =
            new Dictionary<GaussPoint, IShellMaterial>();
        private Matrix _couplingConstitutiveMatrix;
        private Matrix _bendingConstitutiveMatrix;
        private Matrix _membraneConstitutiveMatrix;

        private double[] _normalVector3;
        public double[] NormalVectorV3
        {
            get=>_normalVector3;
            set
            {
                _normalVector3 = value;
                foreach (var material in _thicknessIntegrationPointMaterials.Values)
                {
                    material.NormalVectorV3 = value;
                }
            }
        }
        private double[] _tangentVector1;
        public double[] TangentVectorV1
        {
            get=>_tangentVector1;
            set
            {
                _tangentVector1 = value;
                foreach (var material in _thicknessIntegrationPointMaterials.Values)
                {
                    material.TangentVectorV1 = _tangentVector1;
                }
            }
        }

        private double[] _tangentVector2;
        public double[] TangentVectorV2
        {
            get => _tangentVector2;
            set
            {
                _tangentVector2 = value;
                foreach (var material in _thicknessIntegrationPointMaterials.Values)
                {
                    material.TangentVectorV2 = value;
                }
            }
        }

        public double[] MembraneForces { get; }= new double[3];
        public double[] Moments { get; }= new double[3];
        public double Thickness { get; set; }

        private readonly int _integrationOrder;
        private readonly GaussLegendre1D _thicknessIntegrationPoints;

        public Matrix MembraneConstitutiveMatrix { get; }
        public Matrix BendingConstitutiveMatrix { get; }
        public Matrix CouplingConstitutiveMatrix { get; }
        public void UpdateMaterial(double[] membraneStrains, double[] bendingStrains)
        {
            foreach (var keyValuePair in _thicknessIntegrationPointMaterials)
            {
                var thicknessPoint = keyValuePair.Key;
                var material = keyValuePair.Value;
                var gpStrain = new double[bendingStrains.Length];
                var z = thicknessPoint.Xi;
                for (var i = 0; i < bendingStrains.Length; i++)
                {
                    gpStrain[i] += membraneStrains[i] + bendingStrains[i] * z;
                }

                material.UpdateMaterial(gpStrain);
            }

            MembraneConstitutiveMatrix.Clear();
            BendingConstitutiveMatrix.Clear();
            CouplingConstitutiveMatrix.Clear();

            foreach (var keyValuePair in _thicknessIntegrationPointMaterials)
            {
                var thicknessPoint = keyValuePair.Key;
                var material = keyValuePair.Value;
                var constitutiveMatrixM = material.ConstitutiveMatrix;
                double tempc = 0;
                double w = thicknessPoint.Weight;
                double z = thicknessPoint.Xi;
                for (int i = 0; i < 3; i++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        tempc = constitutiveMatrixM[i, k];
                        MembraneConstitutiveMatrix[i, k] += tempc * w;
                        CouplingConstitutiveMatrix[i, k] += tempc * w * z;
                        BendingConstitutiveMatrix[i, k] += tempc * w * z * z;
                    }
                }
            }

            MembraneForces.Clear();
            Moments.Clear();

            foreach (var keyValuePair in _thicknessIntegrationPointMaterials)
            {
                var thicknessPoint = keyValuePair.Key;
                var material = keyValuePair.Value;

                var w = thicknessPoint.Weight;
                var z = thicknessPoint.Xi;
                MembraneForces[0] += material.Stresses[0] * w;
                MembraneForces[1] += material.Stresses[1] * w;
                MembraneForces[2] += material.Stresses[2] * w;

                Moments[0] -= material.Stresses[0] * w * z;
                Moments[1] -= material.Stresses[1] * w * z;
                Moments[2] -= material.Stresses[2] * w * z;
            }
        }

        object ICloneable.Clone()
        {
            return Clone();
        }

        public int ID { get; }
        public bool Modified { get; private set; }
        public void ResetModified()
        {
            Modified = false;
        }

        public double[] Coordinates { get; set; }
        public double YoungModulus { get; }
        public double PoissonRatio { get; }
        public void SaveState()
        {
        }

        public void ClearState()
        {
        }

        public void ClearStresses()
        {
        }

        public IShellSectionMaterial Clone() =>
            new ShellElasticThicknessMaterial(_integrationOrder, Thickness, YoungModulus, PoissonRatio);
    }
}
