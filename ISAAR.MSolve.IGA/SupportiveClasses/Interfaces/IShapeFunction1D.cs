namespace ISAAR.MSolve.IGA.SupportiveClasses.Interfaces
{
    public interface IShapeFunction1D
    {
        double[,] Values { get; }
        double[,] DerivativeValues { get; }
        double[,] SecondDerivativeValues { get; }
    }
}