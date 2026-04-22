import matplotlib
import numpy
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize
from matplotlib.patches import FancyArrowPatch, ArrowStyle

import metaboglobe
from metaboglobe.math.curve_2d import Curve2
from metaboglobe.math.vector_2d import Vector2
from metaboglobe._util import MPLColor, optimize_for_display, wrap_text
from metaboglobe.kegg_pathway import KeggMap, KeggRelation, KeggReaction, RelationType, ReactionType, EntryType, \
    KeggEntry


class PlottingParameters:

    facecolor: MPLColor = "#eeeeee"
    hide_ticks_and_spines: bool = True

    flux_cmap: Colormap = matplotlib.colormaps.get_cmap("coolwarm")
    flux_vmin: float = 0
    flux_vmax: float = 1
    flux_nan_color: MPLColor = "#888888"
    flux_linewidth: float = 2

    compound_nan_color: MPLColor = "#ffffff"
    compound_edgecolor: MPLColor = "#000000"
    compound_linewidth: float = 0.75
    compound_radius: float = 5

    enzyme_textcolor: MPLColor = "#000000"
    enzyme_linewidth: float = 0.75
    enzyme_facecolor: MPLColor = "#ffffff"
    enzyme_edgecolor: MPLColor = "#dddddd"
    enzyme_padding: float = 0.5
    enzyme_rounding: bool = True

    plot_entries_without_reactions: bool = False

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
        if plotting_parameters.plot_entries_without_reactions or kegg_map.has_relations_or_reactions(entry):
            _draw_entry(ax, entry, plotting_parameters)

    # Draw relations and reactions
    for relation in kegg_map.relations:
        if relation.relation_type == RelationType.MAPLINK:
            _draw_maplink(ax, kegg_map, relation)
    for reaction in kegg_map.reactions:
        _draw_reaction(ax, kegg_map, reaction, plotting_parameters)

    return ScalarMappable(Normalize(vmin=plotting_parameters.flux_vmin, vmax=plotting_parameters.flux_vmax), plotting_parameters.flux_cmap)


def _draw_entry(ax: Axes, entry: KeggEntry, plotting_parameters: PlottingParameters):
    self = entry.entry_type
    if self == EntryType.TITLE:
        ax.set_title(entry.name)
    elif self == EntryType.MAP:
        display_name = wrap_text(entry.name, int(entry.width // 6))  # Wrap text to fit within the entry box, assuming an average character width of 6 pixels
        rect = matplotlib.patches.Rectangle((entry.x - entry.width / 2, entry.y - entry.height / 2),
                                            entry.width, entry.height, fill=True, color="#c1bcb8")
        ax.add_patch(rect)
        ax.text(entry.x, entry.y, display_name, ha="center", va="center", fontsize=6, zorder=10)
    elif self == EntryType.COMPOUND:
        display_name = metaboglobe.kegg_pathway.get_display_name(entry.name)
        ax.text(entry.x, entry.y, display_name, ha="center", va="center", fontsize=6, zorder=10)
        ax.add_patch(matplotlib.patches.Circle((entry.x, entry.y), plotting_parameters.compound_radius,
                                               facecolor=plotting_parameters.compound_nan_color,
                                               edgecolor=plotting_parameters.compound_edgecolor,
                                               linewidth=plotting_parameters.compound_linewidth))
        return
    elif self == EntryType.ORTHOLOG:
        return  # Not interested in orthologs, they just clutter the figure
    elif self == EntryType.GENE:  # Genes
        if "," in entry.name:
            display_name = entry.name.split(",")[0] + ", ..."
        else:
            display_name = entry.name
        bbox = {"facecolor": plotting_parameters.enzyme_facecolor, "edgecolor": plotting_parameters.enzyme_edgecolor,
                "linewidth": plotting_parameters.enzyme_linewidth}
        if plotting_parameters.enzyme_rounding:
            bbox["boxstyle"] = f"round,pad={plotting_parameters.enzyme_padding}"
        else:
            bbox["pad"] = plotting_parameters.enzyme_padding
        ax.text(entry.x, entry.y, display_name, ha="center", va="center", fontsize=6, zorder=10,
                color=plotting_parameters.enzyme_textcolor, bbox=bbox)
    else:
        raise ValueError(f"Unknown entry type: {self}")


def _snap_to_box(entry: Vector2, box_min: Vector2, box_max: Vector2, *, margin_px: float = 2) -> Vector2:
    """If a point is within margin_px from the edge of the box, it is moved to the edge of the box."""
    if abs(entry.x - box_min.x) < margin_px:
        entry = entry.with_x(box_min.x)
    elif abs(entry.x - box_max.x) < margin_px:
        entry = entry.with_x(box_max.x)
    if abs(entry.y - box_min.y) < margin_px:
        entry = entry.with_y(box_min.y)
    elif abs(entry.y - box_max.y) < margin_px:
        entry = entry.with_y(box_max.y)
    return entry


def _draw_reaction(ax: Axes, kegg_map: KeggMap, reaction: KeggReaction, plotting_parameters: PlottingParameters):
    from_entry = kegg_map.entry(reaction.substrate_id).xy()
    to_entry = kegg_map.entry(reaction.product_id).xy()
    enzyme_entry = kegg_map.entry(reaction.gene_id).xy()

    # We take a box around the three objects in our path (the line connecting the three should follow this box)
    box_min = Vector2.min_x_y(from_entry, enzyme_entry, to_entry)
    box_max = Vector2.max_x_y(from_entry, enzyme_entry, to_entry)

    # If any of the points are very close to the box, just move them there
    # (KEGG isn't very precise with aligning things)
    from_entry = _snap_to_box(from_entry, box_min, box_max)
    to_entry = _snap_to_box(to_entry, box_min, box_max)
    enzyme_entry = _snap_to_box(enzyme_entry, box_min, box_max)

    # Check if the enzyme is on the side of a box (necessary for determining the corners in our path)
    enzyme_on_min_or_max_x = enzyme_entry.x == box_min.x or enzyme_entry.x == box_max.x
    enzyme_on_min_or_max_y = enzyme_entry.y == box_min.y or enzyme_entry.y == box_max.y

    # Start drawing the path
    curve = Curve2(from_entry)

    # Make the path to the enzyme
    if enzyme_entry.x == from_entry.x or enzyme_entry.y == from_entry.y:
        # We can draw a straight line to the enzyme
        curve.append_line_to(enzyme_entry)
    else:
        # Need to make a corner
        if enzyme_on_min_or_max_x:  # So enzyme on left or right side of box
            curve.append_rounded_corner_to(corner=enzyme_entry.with_y(from_entry.y), end=enzyme_entry)
        elif enzyme_on_min_or_max_y:  # So enzyme on top or bottom side of box
            curve.append_rounded_corner_to(corner=enzyme_entry.with_x(from_entry.x), end=enzyme_entry)
        else: # Some weird diagonal situation, start with straight line, then rounded line
            curve.append_line_then_curve_to(enzyme_entry)

    # Make the path to the product
    if enzyme_entry.x == to_entry.x or enzyme_entry.y == to_entry.y:
        # We can draw a straight line to the product
        curve.append_line_to(to_entry)
    else:
        # Need to make a corner
        if enzyme_on_min_or_max_x:  # So enzyme on left or right side of box
            curve.append_rounded_corner_to(corner=to_entry.with_x(enzyme_entry.x), end=to_entry)
        elif enzyme_on_min_or_max_y:  # So enzyme on top or bottom side of box
            curve.append_rounded_corner_to(corner=to_entry.with_y(enzyme_entry.y), end=to_entry)
        else: # Some weird diagonal situation, start with rounded line, then straight line
            curve.append_curve_then_line_to(to_entry)

    #curve.shorten(length_start=plotting_parameters.compound_radius, length_end=plotting_parameters.compound_radius)
    style_name = "<|-|>" if reaction.reaction_type == ReactionType.REVERSIBLE else "-|>"
    ax.add_patch(FancyArrowPatch(path=curve.to_path(), arrowstyle=ArrowStyle(style_name, head_length=4, head_width=2)))


    vmin = plotting_parameters.flux_vmin
    vspread = plotting_parameters.flux_vmax - plotting_parameters.flux_vmin
    forward_value = (kegg_map.forward_value(reaction) - vmin) / vspread

def _draw_arrow(ax: Axes, x1: float, y1: float, x2: float, y2: float, cmap: Colormap, value: float, nan_color: MPLColor, linewidth: float, *, arrowstyle: str = "->") -> None:
    if numpy.isnan(value):
        linewidth /= 2
        color = nan_color
    else:
        color = cmap(value)

    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(
        arrowstyle=arrowstyle, color=color, linewidth=linewidth))


def _should_draw_maplink(entry_type_1: EntryType, entry_type_2: EntryType) -> bool:
    return (entry_type_1 == EntryType.COMPOUND and entry_type_2 == EntryType.MAP) \
        or (entry_type_1 == EntryType.MAP and entry_type_2 == EntryType.COMPOUND)


def _draw_maplink(ax: Axes, kegg_map: KeggMap, relation: KeggRelation):
    # We don't draw map links from/to genes, to avoid clutter - so only towards compounds
    from_entry = kegg_map.entry(relation.from_id)
    to_entry = kegg_map.entry(relation.to_id)
    compound = kegg_map.entry(relation.compound_id)

    if _should_draw_maplink(from_entry.entry_type, to_entry.entry_type):
        ax.annotate("", xy=(to_entry.x, to_entry.y), xytext=(from_entry.x, from_entry.y), arrowprops=dict(
            arrowstyle="-|>", color="#ff0000", linestyle="dotted", linewidth=1), zorder=20)

    if relation.from_id != relation.compound_id and _should_draw_maplink(from_entry.entry_type, compound.entry_type):
        ax.annotate("", xy=(compound.x, compound.y), xytext=(from_entry.x, from_entry.y), arrowprops=dict(
            arrowstyle="-|>", color="#00ff00", linestyle="dotted", linewidth=1), zorder=20)

    if relation.to_id != relation.compound_id and _should_draw_maplink(to_entry.entry_type, compound.entry_type):
        ax.annotate("", xy=(to_entry.x, to_entry.y), xytext=(compound.x, compound.y), arrowprops=dict(
            arrowstyle="-|>", color="#0000ff", linestyle="dotted", linewidth=1), zorder=20)
