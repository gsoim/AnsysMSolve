using ISAAR.MSolve.IGA.Interfaces;

namespace ISAAR.MSolve.IGA.Entities.Loads
{
	public class SurfacePressureLoad : ISurfaceLoad
	{
		public SurfacePressureLoad(double pressure) => Pressure = pressure;

		public double Pressure { get; private set; }
	}
}
