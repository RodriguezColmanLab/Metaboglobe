"""Contains code for working with vectors and curves."""
from typing import Literal

import matplotlib
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    flux_label: str | None = None
    flux_vmin: float = 0
    flux_vmax: float = 1
    flux_arrowwidth: float = 1.5
    flux_edgecolor: MPLColor = "#555555"
    flux_arrowsize: float = 5
    flux_edgewidth: float = 0.5  # Width of the edge around the arrow, set to 0 to disable
    flux_joinstyle: MPLJoinStyle = "miter"
    flux_capstyle: MPLCapStyle = "round"
    flux_separation: float = 6  # Used if `self.plot_double_arrows` is True, to separate the two arrows of a double arrow
    flux_nan_color: MPLColor = "#888888"  # Used if no flux is available
    flux_nan_arrowwidth: float = 1  # Used if no flux is available
    flux_nan_separation: float = 4  # Used instead of `self.flux_separation` if no flux is available
    flux_nan_edgecolor: MPLColor = "#000000"  # Used instead of `self.flux_edgecolor` if no flux is available
    flux_nan_edgewidth: float = 0  # Used instead of `self.flux_edgewidth` if no flux is available

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


def plot_kegg_figure(kegg_map: KeggMap, plot_style: PlotStyle | None = None, **kwargs) -> Figure:
    """Plots the given KEGG map on the given axes, and returns a Figure object. To have a colorbar, pass a `flux_label`
    in `plot_style` or as a keyword argument."""
    if plot_style is None:
        plot_style = PlotStyle()
    plot_style.update(**kwargs)

    figsize = kegg_map.figsize()
    figure = Figure(figsize=figsize)
    ax = figure.gca()

    mappable = plot_kegg(ax, kegg_map, plot_style=plot_style)
    if plot_style.flux_label is not None:
        figure.colorbar(mappable, ax=ax, orientation="vertical", label=plot_style.flux_label, shrink=0.2)

    figure.tight_layout()
    return figure
