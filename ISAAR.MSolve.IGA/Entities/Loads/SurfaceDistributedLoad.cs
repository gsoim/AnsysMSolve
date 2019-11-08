using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Interfaces;

namespace ISAAR.MSolve.IGA.Entities.Loads
{
	public class SurfaceDistributedLoad : ISurfaceLoad
	{
		public SurfaceDistributedLoad(double magnitude, IDofType loadedDof)
		{
			Magnitude = magnitude;
			Dof = loadedDof;
		}

		public double Magnitude { get; private set; }

		public IDofType Dof { get; private set; }
	}
}
