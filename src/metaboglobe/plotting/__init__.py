"""Contains code for working with vectors and curves."""
from typing import Literal

import matplotlib
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap

from ..kegg_pathway import KeggMap

MPLColor = tuple[float, float, float] \
           | str \
           | tuple[float, float, float, float] \
           | tuple[tuple[float, float, float] | str, float]

MPLJoinStyle = Literal["miter", "round", "bevel"]
MPLCapStyle = Literal["butt", "projecting", "round"]

class PlotStyle:

    facecolor: MPLColor = "#eeeeee"
    hide_ticks_and_spines: bool = True
    plot_double_arrows: bool = True

    flux_cmap: Colormap = matplotlib.colormaps.get_cmap("coolwarm")
    flux_vmin: float = 0
    flux_vmax: float = 1
    flux_linewidth: float = 2
    flux_arrowsize: float = 5
    flux_joinstyle: MPLJoinStyle = "miter"
    flux_capstyle: MPLCapStyle = "butt"
    flux_nan_color: MPLColor = "#888888"  # Used if no flux is available
    flux_nan_linewidth: float = 1  # Used if no flux is available

    compound_nan_color: MPLColor = "#ffffff"
    compound_edgecolor: MPLColor = "#000000"
    compound_linewidth: float = 0.75
    compound_radius: float = 5

    maplink_edgecolor: MPLColor = "#bbbbbb"
    maplink_linewidth: float = 3
    maplink_linestyle: str = ":"

    enzyme_textcolor: MPLColor = "#000000"
    enzyme_linewidth: float = 0.75
    enzyme_facecolor: MPLColor = "#ffffff"
    enzyme_edgecolor: MPLColor = "#dddddd"
    enzyme_padding: float = 0.25
    enzyme_rounding: bool = True

    plot_entries_without_reactions: bool = False

    def __init__(self, **kwargs):
        """Initializes the plot style with the given parameters. Raises ValueError if any parameter is unknown."""
        self.update(**kwargs)

    def update(self, **kwargs):
        """Updates the plot style with the given parameters. Raises ValueError if any parameter is unknown."""
        for key, value in kwargs.items():
            if not hasattr(self, key):
                raise ValueError(f"Invalid plotting parameter: {key}")

            setattr(self, key, value)


def plot_kegg(ax: Axes, kegg_map: KeggMap, plot_style: PlotStyle | None = None, **kwargs) -> ScalarMappable:
    """Plots the given KEGG map on the given axes, and returns a ScalarMappable that can be used for the colorbar (if any)."""
    from . import _kegg_plotting
    if plot_style is None:
        plot_style = PlotStyle()
    plot_style.update(**kwargs)
    return _kegg_plotting.plot_kegg(ax, kegg_map, plot_style=plot_style)
