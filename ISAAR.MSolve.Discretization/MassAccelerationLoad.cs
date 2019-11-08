using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.Discretization
{
    public class MassAccelerationLoad
    {
        public IDofType DOF { get; set; }
        public double Amount { get; set; }
    }
}
