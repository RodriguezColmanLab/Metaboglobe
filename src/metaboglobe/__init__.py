"""Metaboglobe, for plotting metabolic pathways from KEGG, and coloring pathways by up- and downregulation."""

import importlib

__all__ = ["math", "plotting", "compass_data", "kegg_pathway", "kegg_species", "tabular_data"]

def __getattr__(name):
    """Allows calling things like `metaboglobe.kegg_species.KeggReactionId("rn:R00001")` after just
    `import metaboglobe`, so without specifying the submodule in the import statement."""
    if name in __all__:
        module = importlib.import_module("." + name, __name__)

        # Cache, so that subsequent imports will be faster
        globals()[name] = module
        return module

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
