﻿using System;
using System.Collections.Generic;
using System.Linq;
using Ansys.ACT.Automation.Mechanical;
using Ansys.ACT.Automation.Mechanical.BoundaryConditions;
using Ansys.ACT.Interfaces.Common;
using Ansys.ACT.Interfaces.Geometry;
using Ansys.ACT.Interfaces.Mechanical;
using Ansys.ACT.Interfaces.Mesh;
using Ansys.ACT.Mechanical.Fields;
using Ansys.Core.Units;
using Ansys.EngineeringData.Material;
using Ansys.Mechanical.DataModel.Enums;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using Model = ISAAR.MSolve.FEM.Entities.Model;

namespace AnsysMSolve
{
	public static class AnsysUtilities
	{
		private static readonly Dictionary<ElementTypeEnum, CellType> _ansysMSolveElementDictionary =
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

		private static readonly Dictionary<ElementTypeEnum, int[]> _ansysMSolveElementLocalCoordinates =
			new Dictionary<ElementTypeEnum, int[]>()
			{
				{ElementTypeEnum.kHex8, new[] {4,5,1,0,7,6,2,3}},
				{ElementTypeEnum.kHex20, new int[20] {4,5,1,0,7,6,2,3,12,16,15,17,13,8,9,11,14,19,18,10}},
				{ElementTypeEnum.kTet4, new int[4]{0,3,1,2} },
				{ElementTypeEnum.kTet10, new int[10]{0,3,1,2,7,8,4,6,5,9} },
				{ElementTypeEnum.kWedge6, new int[6]{5,4,3,2,1,0} },
				{ElementTypeEnum.kWedge15, new int[15]{5,4,3,2,1,0,10,11,14,9,13,12,7,8,6} }
			};


		public static ElasticMaterial3D GetAnsysMaterial(IMechanicalExtAPI _api)
		{
			var mat = (_api.DataModel.GeoData.Assemblies[0].Parts[0].Bodies[0] as IGeoBody).Material as MaterialClass;
			var dataDictionary = GetMaterialPropertyByName(mat, "Elasticity");
			var material = new ElasticMaterial3D()
			{
				YoungModulus = dataDictionary["Young's Modulus"].Value,
				PoissonRatio = dataDictionary["Poisson's Ratio"].Value
			};
			return material;
		}

		private static List<string> GetListMaterialProperties(MaterialClass material)
		{
			var propertyNames = new List<string>();
			for (int i = 0; i < material.MaterialProperties.Count; i++)
			{
				propertyNames.Add(material.MaterialProperties.Get(i).TypeName);
			}

			return propertyNames;
		}

		private static Dictionary<string, (string Unit, double Value)> GetMaterialPropertyByName(MaterialClass material, string propertyName)
		{
			var properties = new List<IMaterialProperty>();
			for (var i = 0; i < material.MaterialProperties.Count; i++)
			{
				properties.Add(material.MaterialProperties.Get(i));
			}

			var property = properties.Single(p => p.TypeName == propertyName);

			var variablesDictionary = new Dictionary<string, (string unit, double value)>();
			for (int i = 0; i < property.MaterialPropertyDatas.Get(0).Variables.Count; i++)
			{
				var variable = property.MaterialPropertyDatas.Get(0).Variables.Get(i);
				variablesDictionary.Add(variable.TypeName, (variable.Unit, variable.DatumColl.Get(0).Value));
			}

			return variablesDictionary;
		}

		private static List<Node> RenumberNodesFromAnsysToMSolve(IReadOnlyList<Node> ansysNodes, ElementTypeEnum elementType)
		{
			var connectivity = _ansysMSolveElementLocalCoordinates[elementType];
			return ansysNodes.Select((t, i) => ansysNodes[connectivity[i]]).ToList();
		}


		public static void ImposeFixedSupport(IMechanicalExtAPI _api,IMechanicalUserSolver solver, Model model)
		{
			var ansysFixedSupports =
				_api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
						.Where(c => c.GetType() == typeof(FixedSupport)).ToList()) as List<DataModelObject>;

			foreach (var ansysFixedSupport in ansysFixedSupports)
			{
				var fixedSupport = ansysFixedSupport as FixedSupport;
				var fixedLocation = _api.Application.InvokeUIThread(() => fixedSupport.Location) as ISelectionInfo;
				var fixedSurfaceId = fixedLocation.Ids[0];
				var fixedNodes = solver.Analysis.MeshData.MeshRegionById(fixedSurfaceId).Nodes;

				foreach (var node in fixedNodes)
				{
					model.NodesDictionary[node.Id].Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationX});
                    model.NodesDictionary[node.Id].Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationY});
                    model.NodesDictionary[node.Id].Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationZ});
				}
			}
		}

		public static void GenerateElements(ElasticMaterial3D material, DynamicMaterial dynamicMaterial, IMechanicalUserSolver solver, Model model)
		{
			var factory = new ContinuumElement3DFactory(material, dynamicMaterial);
			foreach (var ansysElement in solver.Analysis.MeshData.Elements)
			{
				var ansysNodes = new List<Node>();
				ansysElement.NodeIds.ToList().ForEach(id => ansysNodes.Add(model.NodesDictionary[id]));
				var msolveNodes = RenumberNodesFromAnsysToMSolve(ansysNodes, ansysElement.Type);
				var element = factory.CreateElement(_ansysMSolveElementDictionary[ansysElement.Type], msolveNodes);
				var elementWrapper = new Element() { ID = ansysElement.Id, ElementType = element };
				foreach (var node in element.Nodes) elementWrapper.AddNode(node);
				model.ElementsDictionary.Add(ansysElement.Id, elementWrapper);
				model.SubdomainsDictionary[0].Elements.Add(elementWrapper);
			}
		}
		
		public static void CalculateForces(IMechanicalExtAPI _api,IMechanicalUserSolver solver, Model model)
		{
			var forces =
				_api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
					.Where(c => c.GetType() == typeof(Force)).ToList()) as List<DataModelObject>;

			foreach (var ansysForce in forces)
			{
				var force = ansysForce as Force;
				var isParsed = Enum.TryParse<LoadDefineBy>(_api.Application.InvokeUIThread(() => force.DefineBy).ToString(), out var defineBy);
				if (defineBy == LoadDefineBy.Vector)
					CalculateVectorForce(_api,solver,model, force);
				else
					CalculateComponentsForce(_api,solver,model, force);
			}
		}

		public static void CalculateDisplacements(IMechanicalExtAPI _api, IMechanicalUserSolver solver, Model model)
		{
			var displacements =
				_api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
					.Where(c => c.GetType() == typeof(Displacement)).ToList()) as List<DataModelObject>;

			foreach (var ansysDisplacements in displacements)
			{
				var displacement = ansysDisplacements as Displacement;
				CalculateComponentsDisplacement(_api, solver, model, displacement);
			}
		}

		public static void CalculateComponentsDisplacement(IMechanicalExtAPI _api, IMechanicalUserSolver solver, Model model, Displacement displacement)
		{
			Quantity xValue = null;
			Quantity yValue = null;
			Quantity zValue = null;
			var displacementLocation = _api.Application.InvokeUIThread(() => displacement.Location) as ISelectionInfo;
			var xFree = (VariableDefinitionType)_api.Application.InvokeUIThread(() => displacement.XComponent.Output.DefinitionType) == VariableDefinitionType.Free;
			var yFree = (VariableDefinitionType)_api.Application.InvokeUIThread(() => displacement.YComponent.Output.DefinitionType) == VariableDefinitionType.Free;
			var zFree = (VariableDefinitionType)_api.Application.InvokeUIThread(() => displacement.ZComponent.Output.DefinitionType) == VariableDefinitionType.Free;
			if (!xFree)
				xValue = _api.Application.InvokeUIThread(() => displacement.XComponent.Output.DiscreteValues[1]) as Quantity;
			if (!yFree)
				yValue = _api.Application.InvokeUIThread(() => displacement.YComponent.Output.DiscreteValues[1]) as Quantity;
			if (!zFree)
				zValue = _api.Application.InvokeUIThread(() => displacement.ZComponent.Output.DiscreteValues[1]) as Quantity;
			var displacementSurfaceId = displacementLocation.Ids[0];
			var displacementNodes = solver.Analysis.MeshData.MeshRegionById(displacementSurfaceId).Nodes;

			foreach (var node in displacementNodes)
			{
				if (!xFree)
				{
					model.NodesDictionary[node.Id].Constraints.Add(new Constraint()
					{
						Amount = xValue.Value,
						DOF = StructuralDof.TranslationX
					});
				}
				if (!yFree)
				{
					model.NodesDictionary[node.Id].Constraints.Add(new Constraint()
					{
						Amount = yValue.Value,
						DOF = StructuralDof.TranslationY
					});
				}
				if (!zFree)
				{
					model.NodesDictionary[node.Id].Constraints.Add(new Constraint()
					{
						Amount = zValue.Value,
						DOF = StructuralDof.TranslationZ
					});
				}
			}
		}
		
		public static void CalculatePressure(IMechanicalExtAPI _api,IMechanicalUserSolver solver, Model model)
        {
            var pressures =
                _api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
                    .Where(c => c.GetType() == typeof(Pressure)).ToList()) as List<DataModelObject>;

            foreach (var ansysPressures in pressures)
            {
                var pressure = ansysPressures as Pressure;
                var isParsed = Enum.TryParse<LoadDefineBy>(_api.Application.InvokeUIThread(() => pressure.DefineBy).ToString(), out var defineBy);

                var forceLocation = _api.Application.InvokeUIThread(() => pressure.Location) as ISelectionInfo;
                var forceSurfaceId = forceLocation.Ids[0];
                var forceNodes = solver.Analysis.MeshData.MeshRegionById(forceSurfaceId).Nodes;
                var a = solver.Analysis.MeshData.MeshRegionById(forceSurfaceId);

                a.Mesh.GetQuad4ExteriorFaces(out int[] quad4Elements);

                solver.Analysis.MeshData.ElementById(quad4Elements[0]);
                var element = solver.Analysis.MeshData.ElementById(quad4Elements[0]);
				
                //if (defineBy == LoadDefineBy.Vector)
                //    CalculateVectorForce(_api,solver,model, force);
                //else
                //    CalculateComponentsForce(_api,solver,model, force);
            }
        }

		//public static void CalculateAcceleration(IMechanicalExtAPI _api, IMechanicalUserSolver solver, Model model)
		//{
		//	var accelerations = _api.Application.InvokeUIThread(() => _api.DataModel.Project.Model.Analyses[0].Children
		//		.Where(c => c.GetType() == typeof(Acceleration)).ToList()) as List<DataModelObject>;

		//	foreach (var ansysAcceleration in accelerations)
		//	{
		//		var acceleration = ansysAcceleration as Acceleration;
		//		var accelerationLocation =
		//			_api.Application.InvokeUIThread(() => acceleration.Location) as ISelectionInfo;

		//		var xValues = _api.Application.InvokeUIThread(() => acceleration.XComponent.Output.DiscreteValues) as List<Quantity>;
		//		var yValues = _api.Application.InvokeUIThread(() => acceleration.YComponent.Output.DiscreteValues) as List<Quantity>;
		//		var zValues = _api.Application.InvokeUIThread(() => acceleration.ZComponent.Output.DiscreteValues) as List<Quantity>;

		//		model.MassAccelerationHistoryLoads.Add(
		//			new MassAccelerationHistoryLoad(xValues.Select(v => v.Value).ToList()) {DOF =  StructuralDof.TranslationX});
		//		model.MassAccelerationHistoryLoads.Add(
		//			new MassAccelerationHistoryLoad(yValues.Select(v => v.Value).ToList()) { DOF =  StructuralDof.TranslationY });
		//		model.MassAccelerationHistoryLoads.Add(
		//			new MassAccelerationHistoryLoad(zValues.Select(v => v.Value).ToList()) { DOF =  StructuralDof.TranslationZ });
		//	}
		//}


		public static void CalculateComponentsForce(IMechanicalExtAPI _api,IMechanicalUserSolver solver,Model model, Force force)
		{
			var forceLocation = _api.Application.InvokeUIThread(() => force.Location) as ISelectionInfo;
			var xValue = _api.Application.InvokeUIThread(() => force.XComponent.Output.DiscreteValues[1]) as Quantity;
			var yValue = _api.Application.InvokeUIThread(() => force.YComponent.Output.DiscreteValues[1]) as Quantity;
			var zValue = _api.Application.InvokeUIThread(() => force.ZComponent.Output.DiscreteValues[1]) as Quantity;
			var forceSurfaceId = forceLocation.Ids[0];
			var forceNodes = solver.Analysis.MeshData.MeshRegionById(forceSurfaceId).Nodes;

			foreach (var node in forceNodes)
			{
				model.Loads.Add(new Load()
				{
					Amount = xValue.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF =StructuralDof.TranslationX
				});
				model.Loads.Add(new Load()
				{
					Amount = yValue.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationY
				});
				model.Loads.Add(new Load()
				{
					Amount = zValue.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationZ
				});
			}
		}

		private static void CalculateVectorForce(IMechanicalExtAPI _api,IMechanicalUserSolver solver, Model model, Force force)
		{
			var forceLocation = _api.Application.InvokeUIThread(() => force.Location) as ISelectionInfo;
			var xValue = _api.Application.InvokeUIThread(() => force.XComponent.Output.DiscreteValues[0]) as Quantity;
			var yValue = _api.Application.InvokeUIThread(() => force.YComponent.Output.DiscreteValues[0]) as Quantity;
			var zValue = _api.Application.InvokeUIThread(() => force.ZComponent.Output.DiscreteValues[0]) as Quantity;
			var magnitude = _api.Application.InvokeUIThread(() => force.Magnitude.Output.DiscreteValues[1]) as Quantity;
			var forceSurfaceId = forceLocation.Ids[0];
			var forceNodes = solver.Analysis.MeshData.MeshRegionById(forceSurfaceId).Nodes;
            var a = solver.Analysis.MeshData.MeshRegionById(forceSurfaceId);
			
			foreach (var node in forceNodes)
			{
				model.Loads.Add(new Load()
				{
					Amount = xValue.Value * magnitude.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationX
				});
				model.Loads.Add(new Load()
				{
					Amount = yValue.Value * magnitude.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationY
				});
				model.Loads.Add(new Load()
				{
					Amount = zValue.Value * magnitude.Value / forceNodes.Count,
					Node = model.NodesDictionary[node.Id],
					DOF = StructuralDof.TranslationZ
				});
			}
		}
	}
}
