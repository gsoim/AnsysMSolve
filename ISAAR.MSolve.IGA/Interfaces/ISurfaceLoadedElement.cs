using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Entities;

namespace ISAAR.MSolve.IGA.Interfaces
{
	public interface ISurfaceLoadedElement
	{
		Dictionary<int, double> CalculateSurfacePressure(Element element, double pressureMagnitude);

		Dictionary<int, double> CalculateSurfaceDistributedLoad(Element element, IDofType loadedDof, double loadMagnitude);
	}
}
