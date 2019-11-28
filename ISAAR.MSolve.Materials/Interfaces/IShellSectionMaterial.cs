using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Materials.Interfaces
{
	public interface IShellSectionMaterial:IFiniteElementMaterial
	{
        double[] NormalVectorV3 { get; set; }
        double[] TangentVectorV1 { get; set; }
        double[] TangentVectorV2 { get; set; }
        new IShellSectionMaterial Clone();
		double[] MembraneForces { get; }
		double[] Moments { get; }
		Matrix MembraneConstitutiveMatrix { get; }
		Matrix BendingConstitutiveMatrix { get; }
		Matrix CouplingConstitutiveMatrix { get; }
		void UpdateMaterial(double[] membraneStrains, double[] bendingStrains);
	}
}