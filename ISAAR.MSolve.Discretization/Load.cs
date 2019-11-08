using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Many boundary (including this) should depend on INode and IElement, in order to work for all discretization methods
namespace ISAAR.MSolve.Discretization
{
    public class Load
    {
        public INode Node { get; set; }
        public IDofType DOF { get; set; }
        public double Amount { get; set; }
    }
}
