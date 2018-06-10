﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Integration.Points;

namespace ISAAR.MSolve.FEM.Integration.Quadratures
{
    /// <summary>
    /// Collection of integration points that are generated by a traditional quadrature rule, independent of the 
    /// element type. Thus these points are cached as static fields of each enum class, so that accessing them is fast 
    /// and storing them is done only once for all elements. All integration points are with respect to a natural 
    /// coordinate system.
    /// </summary>
    public interface IQuadrature1D
    {
        /// <summary>
        /// The order of the integration points is strictly defined for each quadrature.
        /// </summary>
        IReadOnlyList<GaussPoint1D> IntegrationPoints { get; }
    }
}
