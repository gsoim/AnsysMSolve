﻿using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.Discretization.Interfaces
{
	public interface IStructuralAsymmetricModel:IStructuralModel_v2
	{
		IGlobalFreeDofOrdering GlobalRowDofOrdering { get; set; }
		IGlobalFreeDofOrdering GlobalColDofOrdering { get; set; }

		IReadOnlyList<IAsymmetricSubdomain> Subdomains { get; }
	}
}