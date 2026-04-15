import matplotlib
import numpy
from matplotlib.axes import Axes
from matplotlib.colors import Colormap

from metaboglobe._util import MPLColor, point_direction
from metaboglobe.kegg_pathway import KeggMap, KeggRelation, RelationType


def plot_double_arrows(ax: Axes, kegg_map: KeggMap, *, facecolor: MPLColor = "#eeeeee",
                       hide_ticks_and_spines: bool = True, reaction_cmap: Colormap | None = None,
                       reaction_nan_color: MPLColor = "#888888", reaction_linewidth: float = 2,
                       plot_entries_without_relations: bool = False):
    """Draws a plot. Reactions """

    if reaction_cmap is None:
        reaction_cmap = matplotlib.colormaps.get_cmap("coolwarm")

    # Set up plot
    ax.set_facecolor(facecolor)
    ax.set_xlim(0, 1000)
    ax.set_ylim(1000, 0)
    if hide_ticks_and_spines:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

    # Draw entries
    for entry in kegg_map.entries:
        if plot_entries_without_relations or kegg_map.has_relations(entry):
            entry.entry_type.draw_entry(ax, entry)

    # Draw relations
    for relation in kegg_map.relations:
        _draw_relation(ax, kegg_map, relation, cmap=reaction_cmap, reaction_linewidth=reaction_linewidth,
                       reaction_nan_color=reaction_nan_color)


def _draw_relation(ax: Axes, kegg_map: KeggMap, relation: KeggRelation, *, cmap: Colormap, reaction_linewidth: float, reaction_nan_color: MPLColor):
    if relation.relation_type == RelationType.ECREL:
        return  # Not interested in subsequent enzymes
    elif relation.relation_type == RelationType.MAPLINK:
        _draw_maplink(ax, kegg_map, relation)
    elif relation.relation_type.is_reaction():
        _draw_reaction(ax, kegg_map, relation, cmap=cmap, linewidth=reaction_linewidth, nan_color=reaction_nan_color)


def _draw_reaction(ax: Axes, kegg_map: KeggMap, relation: KeggRelation, cmap: Colormap, nan_color: MPLColor, linewidth: float):
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

    if relation.relation_type in {RelationType.REACTION_REVERSIBLE, RelationType.REACTION_TWO_IRREVERSIBLE}:
        # Need to draw two arrows
        direction = point_direction(x1, y1, x2, y2)
        orthogonal_direction = direction.orthogonal()

        arrow_distance = 5 if relation.relation_type == RelationType.REACTION_REVERSIBLE else 8

        x1_forward = x1 + orthogonal_direction.dx() * arrow_distance
        x2_forward = x2 + orthogonal_direction.dx() * arrow_distance
        y1_forward = y1 + orthogonal_direction.dy() * arrow_distance
        y2_forward = y2 + orthogonal_direction.dy() * arrow_distance

        x1_backward = x2 - orthogonal_direction.dx() * arrow_distance
        x2_backward = x1 - orthogonal_direction.dx() * arrow_distance
        y1_backward = y2 - orthogonal_direction.dy() * arrow_distance
        y2_backward = y1 - orthogonal_direction.dy() * arrow_distance

        forward_value = kegg_map.forward_value(relation)
        _draw_arrow(ax, x1_forward, y1_forward, x2_forward, y2_forward, cmap, forward_value, nan_color=nan_color, linewidth=linewidth)

        backward_value = kegg_map.backward_value(relation)
        _draw_arrow(ax, x1_backward, y1_backward, x2_backward, y2_backward, cmap, backward_value, nan_color=nan_color, linewidth=linewidth)
    else:
        forward_value = kegg_map.forward_value(relation)
        _draw_arrow(ax, x1, y1, x2, y2, cmap, forward_value, nan_color=nan_color, linewidth=linewidth)

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