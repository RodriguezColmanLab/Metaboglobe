from typing import NamedTuple

import matplotlib
import numpy
from numpy import ndarray
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize

from metaboglobe._util import MPLColor, point_direction
from metaboglobe.kegg_pathway import KeggMap, KeggRelation, RelationType


class PlottingParameters:

    facecolor: MPLColor = "#eeeeee"
    hide_ticks_and_spines: bool = True

    flux_cmap: Colormap = matplotlib.colormaps.get_cmap("coolwarm")
    flux_vmin: float = 0
    flux_vmax: float = 1
    flux_nan_color: MPLColor = "#888888"
    flux_linewidth: float = 2

    plot_entries_without_relations: bool = False

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            # We just check whether the attribute exists, we don't type check (would be nice to have though)
            if not hasattr(self, key):
                raise ValueError(f"Invalid plotting parameter: {key}")

            setattr(self, key, value)


def _adjust_limits(ax: Axes, kegg_map: KeggMap):
    """Adjusts the xlim and ylim of the plot to fit the map."""

    # Calculate necessary space for plot
    xmax = 0
    ymax = 0
    for entry in kegg_map.entries:
        xmax = max(entry.x, xmax)
        ymax = max(entry.y, ymax)

    if xmax == 0 or ymax == 0:
        return  # Invalid limits

    # Give some extra space
    xsize = xmax
    ysize = ymax
    xmax += xsize / 10
    ymax += ysize / 10

    # Adjust the limits
    ax.set_xlim(0, xmax)
    ax.set_ylim(ymax, 0)
    ax.set_aspect("equal")


def plot_double_arrows(ax: Axes, kegg_map: KeggMap, **kwargs) -> ScalarMappable:
    """Draws the KEGG map, with double arrows for revisble/two-way-irrervisble reactions. Returns a mappable for use
    in figure.colorbar(...).

    Any kwargs parameters are passed on to PlottingParameters, see that class for available parameters.
    """
    plotting_parameters = PlottingParameters(**kwargs)

    # Set up plot
    ax.set_facecolor(plotting_parameters.facecolor)
    _adjust_limits(ax, kegg_map)
    if plotting_parameters.hide_ticks_and_spines:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

    # Draw entries
    for entry in kegg_map.entries:
        if plotting_parameters.plot_entries_without_relations or kegg_map.has_relations(entry):
            entry.entry_type.draw_entry(ax, entry)

    # Draw relations
    for relation in kegg_map.relations:
        _draw_relation(ax, kegg_map, relation, plotting_parameters)

    return ScalarMappable(Normalize(vmin=plotting_parameters.flux_vmin, vmax=plotting_parameters.flux_vmax), plotting_parameters.flux_cmap)


def _draw_relation(ax: Axes, kegg_map: KeggMap, relation: KeggRelation, plotting_parameters: PlottingParameters):
    if relation.relation_type == RelationType.ECREL:
        return  # Not interested in subsequent enzymes
    elif relation.relation_type == RelationType.MAPLINK:
        _draw_maplink(ax, kegg_map, relation)
    elif relation.relation_type.is_reaction():
        _draw_reaction(ax, kegg_map, relation, plotting_parameters)


def _draw_reaction(ax: Axes, kegg_map: KeggMap, relation: KeggRelation, plotting_parameters: PlottingParameters):
    from_entry = kegg_map.entry(relation.substrate_id)
    to_entry = kegg_map.entry(relation.product_id)
    x1 = from_entry.x
    y1 = from_entry.y
    x2 = to_entry.x
    y2 = to_entry.y

    # Make arrow shorter
    # Horizontal arrow correction
    if abs(x1 - x2) > 50:  # Only apply if the arrow is long enough, otherwise it looks weird
        if x1 > x2:
            x1 -= from_entry.width / 2
            x2 += to_entry.width / 2
        else:
            x1 += from_entry.width / 2
            x2 -= to_entry.width / 2

    # Vertical arrow correction
    if abs(y1 - y2) > 20:  # Only apply if the arrow is long enough, otherwise it looks weird
        if y1 > y2:
            y1 -= from_entry.height / 2
            y2 += to_entry.height / 2
        else:
            y1 += from_entry.height / 2
            y2 -= to_entry.height / 2

    vmin = plotting_parameters.flux_vmin
    vspread = plotting_parameters.flux_vmax - plotting_parameters.flux_vmin
    if relation.relation_type in {RelationType.REACTION_REVERSIBLE, RelationType.REACTION_TWO_IRREVERSIBLE}:
        # Need to draw two arrows
        direction = point_direction(x1, y1, x2, y2)
        orthogonal_direction = direction.orthogonal()

        arrow_distance = 5 if relation.relation_type == RelationType.REACTION_REVERSIBLE else 12

        x1_forward = x1 + orthogonal_direction.dx() * arrow_distance
        x2_forward = x2 + orthogonal_direction.dx() * arrow_distance
        y1_forward = y1 + orthogonal_direction.dy() * arrow_distance
        y2_forward = y2 + orthogonal_direction.dy() * arrow_distance

        x1_backward = x2 - orthogonal_direction.dx() * arrow_distance
        x2_backward = x1 - orthogonal_direction.dx() * arrow_distance
        y1_backward = y2 - orthogonal_direction.dy() * arrow_distance
        y2_backward = y1 - orthogonal_direction.dy() * arrow_distance

        forward_value = (kegg_map.forward_value(relation) - vmin) / vspread
        _draw_arrow(ax, x1_forward, y1_forward, x2_forward, y2_forward, plotting_parameters.flux_cmap, forward_value,
                    nan_color=plotting_parameters.flux_nan_color, linewidth=plotting_parameters.flux_linewidth)

        backward_value = (kegg_map.backward_value(relation) - vmin) / vspread
        _draw_arrow(ax, x1_backward, y1_backward, x2_backward, y2_backward, plotting_parameters.flux_cmap, backward_value,
                    nan_color=plotting_parameters.flux_nan_color, linewidth=plotting_parameters.flux_linewidth)
    else:
        forward_value = (kegg_map.forward_value(relation) - vmin) / vspread
        _draw_arrow(ax, x1, y1, x2, y2, plotting_parameters.flux_cmap, forward_value,
                    nan_color=plotting_parameters.flux_nan_color, linewidth=plotting_parameters.flux_linewidth)

def _draw_arrow(ax: Axes, x1: float, y1: float, x2: float, y2: float, cmap: Colormap, value: float, nan_color: MPLColor, linewidth: float):
    if numpy.isnan(value):
        linewidth /= 2
        color = nan_color
    else:
        color = cmap(value)

    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(
        arrowstyle="->", color=color, linewidth=linewidth))


def _draw_maplink(ax: Axes, kegg_map: KeggMap, relation: KeggRelation):
    from_entry = kegg_map.entry(relation.substrate_id)
    to_entry = kegg_map.entry(relation.product_id)
    x1 = from_entry.x
    y1 = from_entry.y
    x2 = to_entry.x
    y2 = to_entry.y

    # Make arrow shorter
    # Horizontal arrow correction
    if abs(x1 - x2) > 50:  # Only apply if the arrow is long enough, otherwise it looks weird
        if x1 > x2:
            x1 -= from_entry.width / 2
            x2 += to_entry.width / 2
        else:
            x1 += from_entry.width / 2
            x2 -= to_entry.width / 2

    # Vertical arrow correction
    if abs(y1 - y2) > 20:  # Only apply if the arrow is long enough, otherwise it looks weird
        if y1 > y2:
            y1 -= from_entry.height / 2
            y2 += to_entry.height / 2
        else:
            y1 += from_entry.height / 2
            y2 -= to_entry.height / 2

    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(
        arrowstyle="->", color="#888888", linestyle="--", linewidth=1))