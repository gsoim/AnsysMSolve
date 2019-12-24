namespace ISAAR.MSolve.IGA.SupportiveClasses.Interfaces
{
    public interface IShapeFunction3D
    {
        double[,] Values { get; }

        double[,] DerivativeValuesKsi { get; }
        double[,] DerivativeValuesHeta { get; }
        double[,] DerivativeValuesZeta { get; }

        double[,] SecondDerivativeValuesKsi { get; }
        double[,] SecondDerivativeValuesHeta { get; }
        double[,] SecondDerivativeValuesZeta { get; }

        double[,] SecondDerivativeValuesKsiHeta { get; }
        double[,] SecondDerivativeValuesHetaZeta { get; }
        double[,] SecondDerivativeValuesKsiZeta { get; }
    }
}