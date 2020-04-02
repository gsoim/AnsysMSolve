using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using Ansys.ACT.Automation.Mechanical;
using Ansys.ACT.Automation.Mechanical.BoundaryConditions;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.DataModel;
using Ansys.ACT.Interfaces.Geometry;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces.Mesh;
using Ansys.ACT.Interfaces.UserObject;
using Ansys.EngineeringData.Material;
using Ansys.Mechanical.DataModel.Enums;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Iterative;
using Newtonsoft.Json;
using Microsoft.Scripting.Hosting;
using Model = ISAAR.MSolve.FEM.Entities.Model;

namespace AnsysMSolve
{
	public class MSolveStatic:IExtension
	{
		private readonly IMechanicalExtAPI _api;
		private IMechanicalUserSolver _solver;
		private Dictionary<ElementTypeEnum, Dictionary<int, int[]>> ElementTypesNodesDictionary =
			new Dictionary<ElementTypeEnum, Dictionary<int, int[]>>()
			{
				{
					ElementTypeEnum.kTri3,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1}},
						{1, new int[] {0, 2}},
						{2, new int[] {1, 0}}
					}
				},
				{
					ElementTypeEnum.kTri6,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1}},
						{1, new int[] {0, 2}},
						{2, new int[] {1, 0}},
						{3, new int[] {0, 1}},
						{4, new int[] {1, 2}},
						{5, new int[] {2, 3}}
					}
				},
				{
					ElementTypeEnum.kQuad4,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1}},
						{1, new int[] {0, 2}},
						{2, new int[] {1, 3}},
						{3, new int[] {2, 0}},
					}
				},
				{
					ElementTypeEnum.kQuad8,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1}},
						{1, new int[] {0, 2}},
						{2, new int[] {1, 3}},
						{3, new int[] {2, 0}},
						{4, new int[] {0, 1}},
						{5, new int[] {1, 2}},
						{6, new int[] {2, 3}},
						{7, new int[] {3, 0}},
					}
				},
				{
					ElementTypeEnum.kTet4,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1, 3}},
						{1, new int[] {0, 2, 3}},
						{2, new int[] {1, 0, 3}},
						{3, new int[] {0, 1, 2}}
					}
				},
				{
					ElementTypeEnum.kTet10,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1, 3}},
						{1, new int[] {0, 2, 3}},
						{2, new int[] {1, 0, 3}},
						{3, new int[] {0, 1, 2}},
						{4, new int[] {0, 1}},
						{5, new int[] {1, 2}},
						{6, new int[] {2, 0}},
						{7, new int[] {0, 3}},
						{8, new int[] {1, 3}},
						{9, new int[] {2, 3}},
					}
				},
				{
					ElementTypeEnum.kPyramid5,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1, 4}},
						{1, new int[] {0, 2, 4}},
						{2, new int[] {1, 3, 4}},
						{3, new int[] {2, 0, 4}},
						{4, new int[] {0, 1, 2, 3}},
					}
				},
				{
					ElementTypeEnum.kPyramid13,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1, 4}},
						{1, new int[] {0, 2, 4}},
						{2, new int[] {1, 3, 4}},
						{3, new int[] {2, 0, 4}},
						{4, new int[] {0, 1, 2, 3}},
						{5, new int[] {0, 1}},
						{6, new int[] {1, 2}},
						{7, new int[] {2, 3}},
						{8, new int[] {3, 0}},
						{9, new int[] {0, 4}},
						{10, new int[] {1, 4}},
						{11, new int[] {2, 4}},
						{12, new int[] {3, 4}},
					}
				},
				{
					ElementTypeEnum.kWedge6,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1, 3}},
						{1, new int[] {0, 2, 4}},
						{2, new int[] {1, 0, 5}},
						{3, new int[] {5, 4, 0}},
						{4, new int[] {3, 5, 1}},
						{5, new int[] {4, 3, 2}},
					}
				},
				{
					ElementTypeEnum.kWedge15,
					new Dictionary<int, int[]>
					{
						{0, new int[] {2, 1, 3}},
						{1, new int[] {0, 2, 4}},
						{2, new int[] {1, 0, 5}},
						{3, new int[] {5, 4, 0}},
						{4, new int[] {3, 5, 1}},
						{5, new int[] {4, 3, 2}},
						{6, new int[] {0, 1}},
						{7, new int[] {1, 2}},
						{8, new int[] {2, 0}},
						{9, new int[] {3, 4}},
						{10, new int[] {4, 5}},
						{11, new int[] {5, 3}},
						{12, new int[] {0, 3}},
						{13, new int[] {1, 4}},
						{14, new int[] {2, 5}},
					}
				},
				{
					ElementTypeEnum.kHex8,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1, 4}},
						{1, new int[] {0, 2, 5}},
						{2, new int[] {1, 3, 6}},
						{3, new int[] {2, 0, 7}},
						{4, new int[] {5, 0, 7}},
						{5, new int[] {1, 4, 6}},
						{6, new int[] {5, 7, 2}},
						{7, new int[] {6, 3, 4}}
					}
				},
				{
					ElementTypeEnum.kHex20,
					new Dictionary<int, int[]>
					{
						{0, new int[] {3, 1, 4}},
						{1, new int[] {0, 2, 5}},
						{2, new int[] {1, 3, 6}},
						{3, new int[] {2, 0, 7}},
						{4, new int[] {5, 0, 7}},
						{5, new int[] {1, 4, 6}},
						{6, new int[] {5, 7, 2}},
						{7, new int[] {6, 3, 4}},
						{8, new int[] {0, 1}},
						{9, new int[] {1, 2}},
						{10, new int[] {2, 3}},
						{11, new int[] {3, 0}},
						{12, new int[] {4, 5}},
						{13, new int[] {5, 6}},
						{14, new int[] {6, 7}},
						{15, new int[] {7, 4}},
						{16, new int[] {0, 4}},
						{17, new int[] {1, 5}},
						{18, new int[] {2, 6}},
						{19, new int[] {3, 7}},
					}
				},
			};

		private Dictionary<ElementTypeEnum, CellType> AnsysMSolveElementDictionary =
			new Dictionary<ElementTypeEnum, CellType>
			{
				{ElementTypeEnum.kHex8, CellType.Hexa8},
				{ElementTypeEnum.kHex20, CellType.Hexa20},
				{ElementTypeEnum.kTet4, CellType.Tet4},
				{ElementTypeEnum.kTet10, CellType.Tet10},
				{ElementTypeEnum.kWedge6, CellType.Wedge6},
				{ElementTypeEnum.kWedge15, CellType.Wedge15},
				{ElementTypeEnum.kPyramid5, CellType.Pyra5},
				{ElementTypeEnum.kPyramid13, CellType.Pyra13},
			};

		public MSolveStatic(IExtAPI api, IUserSolver solver)
		{
			_api = api as IMechanicalExtAPI;
            _solver = (IMechanicalUserSolver)solver;
        }

        private const string guid = "EB5D5C22-DF17-4755-AFAE-813C733F4791";

		public virtual bool onsolve(IUserSolver userSolver, Func<int, string, bool> progressFunction)
		{
			
			var solver = userSolver as IMechanicalUserSolver;
			try
			{
				var model = new Model();
				model.SubdomainsDictionary.Add(0, new Subdomain(0));
				solver.Analysis.MeshData.Nodes.ToList()
					.ForEach(n => model.NodesDictionary.Add(n.Id, new Node(n.Id, n.X, n.Y, n.Z)));
				var mat = (MaterialClass)(_api.DataModel.GeoData.Assemblies[0].Parts[0].Bodies[0] as IGeoBody).Material;
				
                var materialProperties = mat.MaterialProperties;
                var elasticity=materialProperties.Get(12);
                var materialDatas = elasticity.MaterialPropertyDatas.Get(0);
                var youngModulus=materialDatas.Variables.Get(0).MaxValue;
                var poissonRatio=materialDatas.Variables.Get(1).MaxValue;

                var elasticMaterial3D = new ElasticMaterial3D()
                {
                    YoungModulus = youngModulus,
                    PoissonRatio = poissonRatio
                };
				var factory = new ContinuumElement3DFactory(elasticMaterial3D, null);

				foreach (var ansysElement in solver.Analysis.MeshData.Elements)
				{
					var elementNodes = new List<Node>();
					ansysElement.NodeIds.ToList().ForEach(id => elementNodes.Add(model.NodesDictionary[id]));
					var element = factory.CreateElement(AnsysMSolveElementDictionary[ansysElement.Type], elementNodes);
					var elementWrapper = new Element() { ID = ansysElement.Id, ElementType = element };
					foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
					model.ElementsDictionary.Add(ansysElement.Id, elementWrapper);
					model.SubdomainsDictionary[0].Elements.Add(elementWrapper);
				}

				AnsysUtilities.CalculateDisplacements(_api, solver, model);
				AnsysUtilities.CalculateForces(_api, solver, model);
				AnsysUtilities.CalculatePressure(_api, _solver, model);

                //AnsysUtilities.ImposeFixedSupport(_api,solver,model);

                //var analysis = ((Analysis) solver.Analysis);
                //var analysisChildren = analysis.Children;
                //var force=analysisChildren.FirstOrDefault(c => c is Force);
                //var props=_api.Application.InvokeUIThread(() => force.VisibleProperties);
                //var dataObjects = analysis.DataObjects.List;
                //var loads=solver.Analysis.GetLoadObjects("AnsysMSolve");
                //var loads1=solver.Analysis.GetLoadObjects(solver.Extension);
                
                //var loads1 = solver.Analysis.GetLoadObjects("ansysmsolve");
                //var forces = dataObjects.Where(o => o.Name == "Force");
                //var support = dataObjects.Where(o => o.Name == "FixedSupport");
                //var analysisChildren = ((Analysis) solver.Analysis).Children;
				             

                var solverBuilder = new PcgSolver.Builder();
				ISolver msolveSolver = solverBuilder.BuildSolver(model);
				var provider = new ProblemStructural(model, msolveSolver);
				var childAnalyzer = new LinearAnalyzer(model, msolveSolver, provider);
				var parentAnalyzer = new StaticAnalyzer(model, msolveSolver, provider, childAnalyzer);
				parentAnalyzer.Initialize();
				parentAnalyzer.Solve();
			}
			catch (Exception e)
			{
				return false;
			}

			return true;
		}

		
		private List<TempLoad> _loads = new List<TempLoad>();

		public ObjectTypeEnum ObjectType { get; }
		public ObjectLocationEnum ObjectLocation { get; }
		public IDataModelObject Parent { get; }
		public ObjectLocationEnum[] SupportedLocations { get; }
		public void NotifyChange()
		{
			throw new NotImplementedException();
		}

		public int Id { get; }
		public string Name { get; }
		public int Version { get; }
		public string State { get; }
		public string Caption { get; }
		public string Icon { get; }
		public IExtension Extension { get; }
		public object Controller { get; }
		public IAttributeCollection Attributes { get; }
		public IPropertyCollection Properties { get; }
		public IEnumerable<ISimProperty> AllProperties { get; }
		public IEnumerable<IUserObject> Children { get; }
		public string AnalysisType { get; }
		public string PhysicsType { get; }

        public int MinorVersion => 1;

		public Microsoft.Scripting.Hosting.ScriptEngine ScriptEngine => throw new NotImplementedException();

		public Microsoft.Scripting.Hosting.ScriptScope ScriptScope => throw new NotImplementedException();

        public string UniqueId => "EB5D5C22-DF17-4755-AFAE-813C733F4791";

        public string InstallDir => Directory.GetCurrentDirectory();
    }

	public class TempLoad
	{
		public double Amount { get; set; }
		public IDofType DofType { get; set; }
		public int nodeId { get; set; }
	}


}
