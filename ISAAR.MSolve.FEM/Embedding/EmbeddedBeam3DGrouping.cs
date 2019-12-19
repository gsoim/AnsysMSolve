using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Embedding
{
    public class EmbeddedBeam3DGrouping
    {
        private readonly Model model;
        private readonly IEnumerable<Element> hostGroup;
        private readonly IEnumerable<Element> embeddedGroup;
        private readonly bool hasEmbeddedRotations = false;
        private readonly int skip;
        public IEnumerable<Element> HostGroup { get { return hostGroup; } }
        public IEnumerable<Element> EmbeddedGroup { get { return embeddedGroup; } }

        public static EmbeddedBeam3DGrouping CreateFullyBonded(Model  model, IEnumerable<Element> hostGroup, IEnumerable<Element> embeddedGroup, bool hasEmbeddedRotations)
        {
            return new EmbeddedBeam3DGrouping(model, hostGroup, embeddedGroup, hasEmbeddedRotations, 0);
        }

        public static EmbeddedBeam3DGrouping CreateCohesive(Model  model, IEnumerable<Element> hostGroup, IEnumerable<Element> embeddedGroup, bool hasEmbeddedRotations)
        {
            return new EmbeddedBeam3DGrouping(model, hostGroup, embeddedGroup, hasEmbeddedRotations, 2);
        }

        private EmbeddedBeam3DGrouping(Model  model, IEnumerable<Element> hostGroup, IEnumerable<Element> embeddedGroup, bool hasEmbeddedRotations, int skip)
        {
            this.model = model;
            this.hostGroup = hostGroup;
            this.embeddedGroup = embeddedGroup;
            this.hasEmbeddedRotations = hasEmbeddedRotations;
            this.skip = skip;
            hostGroup.Select(e => e.ElementType).Distinct().ToList().ForEach(et =>
            {
                if (!(et is IEmbeddedHostElement))
                    throw new ArgumentException("EmbeddedGrouping: One or more elements of host group does NOT implement IEmbeddedHostElement.");
            });
            embeddedGroup.Select(e => e.ElementType).Distinct().ToList().ForEach(et =>
            {
                if (!(et is IEmbeddedElement))
                    throw new ArgumentException("EmbeddedGrouping: One or more elements of embedded group does NOT implement IEmbeddedElement.");
            });
            UpdateNodesBelongingToEmbeddedElements();
        }
        

        private void UpdateNodesBelongingToEmbeddedElements()
        {
            IEmbeddedDOFInHostTransformationVector transformer;
            if (hasEmbeddedRotations)
                transformer = new Hexa8TranslationAndRotationTransformationVector();
            else
                transformer = new Hexa8LAndNLTranslationTransformationVector();

            foreach (var embeddedElement in embeddedGroup)
            {
                var elType = (IEmbeddedElement)embeddedElement.ElementType;
                foreach (var node in embeddedElement.Nodes.Skip(skip))
                {
                    var embeddedNodes = hostGroup
                        .Select(e => ((IEmbeddedHostElement)e.ElementType).BuildHostElementEmbeddedNode(e, node, transformer))
                        .Where(e => e != null);
                    foreach (var embeddedNode in embeddedNodes)
                    {
                        if (elType.EmbeddedNodes.Count(x => x.Node == embeddedNode.Node) == 0)
                            elType.EmbeddedNodes.Add(embeddedNode);

                        // Update embedded node information for elements that are not inside the embedded group but contain an embedded node.
                        foreach (var element in model.Elements.Except(embeddedGroup))
                            if (element.ElementType is IEmbeddedElement && element.Nodes.Contains(embeddedNode.Node))
                            {
                                var currentElementType = (IEmbeddedElement)element.ElementType;
                                if (!currentElementType.EmbeddedNodes.Contains(embeddedNode))
                                {
                                    currentElementType.EmbeddedNodes.Add(embeddedNode);
                                    element.ElementType.DofEnumerator = new BeamElementEmbedder(model, element, transformer);
                                }
                            }
                    }
                }
                embeddedElement.ElementType.DofEnumerator = new BeamElementEmbedder(model, embeddedElement, transformer);
            }
        }
    }
}
