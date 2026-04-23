import matplotlib
import numpy
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize
from matplotlib.patches import FancyArrowPatch, ArrowStyle, PathPatch
from matplotlib.path import Path

import metaboglobe
from metaboglobe.plotting._collision_map import CollisionMap, TextWithAnchor
from metaboglobe.plotting._curve_2d import Curve2
from metaboglobe.plotting._vector_2d import Vector2
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

    maplink_edgecolor: MPLColor = "#bbbbbb"
    maplink_linewidth: float = 3
    maplink_linestyle: str = ":"

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
            _draw_maplinks(ax, kegg_map, relation, plotting_parameters)
    for reaction in kegg_map.reactions:
        _draw_reaction(ax, kegg_map, reaction, plotting_parameters)

    # Draw compound names
    _draw_compound_names(ax, kegg_map, plotting_parameters)

    return ScalarMappable(Normalize(vmin=plotting_parameters.flux_vmin, vmax=plotting_parameters.flux_vmax), plotting_parameters.flux_cmap)


def _draw_compound_names(ax: Axes, kegg_map: KeggMap, plotting_parameters: PlottingParameters):
    collision_map = CollisionMap(ax)
    for artist in ax.get_children():
        collision_map.add_artist(artist)

    texts = list()
    for entry in kegg_map.entries:
        if entry.entry_type == EntryType.COMPOUND:
            if plotting_parameters.plot_entries_without_reactions or kegg_map.has_relations_or_reactions(entry):
                display_name = metaboglobe.kegg_pathway.get_display_name(entry.name)
                texts.append(TextWithAnchor(text=display_name, x=entry.x, y=entry.y))
    collision_map.fit_text(ax, texts, fontsize=6, zorder=10)


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


def _to_vector(entry: KeggEntry) -> Vector2:
    return Vector2(entry.x, entry.y)


def _draw_reaction(ax: Axes, kegg_map: KeggMap, reaction: KeggReaction, plotting_parameters: PlottingParameters):
    from_entry = _to_vector(kegg_map.entry(reaction.substrate_id))
    to_entry = _to_vector(kegg_map.entry(reaction.product_id))
    enzyme_entry = _to_vector(kegg_map.entry(reaction.gene_id))

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


def _draw_maplink(ax: Axes, entry1: KeggEntry, entry2: KeggEntry, plotting_parameters: PlottingParameters) -> None:
    if not( (entry1.entry_type == EntryType.COMPOUND and entry2.entry_type == EntryType.MAP)
        or (entry1.entry_type == EntryType.MAP and entry2.entry_type == EntryType.COMPOUND)):
        # Only draw maplinks between maps and compounds, not any of the other specified ones
        return

    map_entry = entry1
    compound_entry = entry2
    if map_entry.entry_type == EntryType.COMPOUND:
        # Switch them around so that names match
        compound_entry, map_entry = compound_entry, map_entry

    map_min_x = map_entry.x - map_entry.width / 2
    map_max_x = map_entry.x + map_entry.width / 2
    map_min_y = map_entry.y - map_entry.height / 2
    map_max_y = map_entry.y + map_entry.height / 2

    if map_min_x <= compound_entry.x <= map_max_x:
        # Easy, we can draw a straight vertical line
        to_y = map_entry.y - map_entry.height / 2 if compound_entry.y < map_entry.y else map_entry.y + map_entry.height / 2
        patch = PathPatch(Path([(compound_entry.x, compound_entry.y), (compound_entry.x, to_y)]),
                          lw=plotting_parameters.maplink_linewidth, edgecolor=plotting_parameters.maplink_edgecolor,
                          linestyle=plotting_parameters.maplink_linestyle, zorder=-5)
        ax.add_patch(patch)
        return
    if map_min_y <= compound_entry.y <= map_max_y:
        # Also easy, we can draw a straight horizontal line
        to_x = map_entry.x - map_entry.width / 2 if compound_entry.x < map_entry.x else map_entry.x + map_entry.width / 2
        patch = PathPatch(Path([(compound_entry.x, compound_entry.y), (to_x, compound_entry.y)]),
                          lw=plotting_parameters.maplink_linewidth, edgecolor=plotting_parameters.maplink_edgecolor,
                          linestyle=plotting_parameters.maplink_linestyle, zorder=-5)
        ax.add_patch(patch)
        return

    # Define all possible paths (always first attachment to pathway rectangle, then the corner)
    possible_paths = [[Vector2(map_min_x, map_entry.y), Vector2(compound_entry.x, map_entry.y)],
                      [Vector2(map_max_x, map_entry.y), Vector2(compound_entry.x, map_entry.y)],
                      [Vector2(map_entry.x, map_min_y), Vector2(map_entry.x, compound_entry.y)],
                      [Vector2(map_entry.x, map_max_y), Vector2(map_entry.x, compound_entry.y)]]
    shortest_distance_squared = float("inf")
    shortest_distance_path = None
    for attachment, path_corner in possible_paths:
        distance_squared = (attachment.x - compound_entry.x) ** 2 + (attachment.y - compound_entry.y) ** 2
        if distance_squared < shortest_distance_squared:
            shortest_distance_squared = distance_squared
            shortest_distance_path = attachment, path_corner

    curve = Curve2(Vector2(compound_entry.x, compound_entry.y))
    curve.append_cut_corner_to(corner=shortest_distance_path[1], end=shortest_distance_path[0])
    patch = PathPatch(curve.to_path(), lw=plotting_parameters.maplink_linewidth, edgecolor=plotting_parameters.maplink_edgecolor,
                      linestyle=plotting_parameters.maplink_linestyle, zorder=-5, fill=False)
    ax.add_patch(patch)


def _draw_maplinks(ax: Axes, kegg_map: KeggMap, relation: KeggRelation, plotting_parameters: PlottingParameters) -> None:
    # We don't draw map links from/to genes, to avoid clutter - so only towards compounds
    from_entry = kegg_map.entry(relation.from_id)
    to_entry = kegg_map.entry(relation.to_id)
    compound = kegg_map.entry(relation.compound_id)

    _draw_maplink(ax, from_entry, to_entry, plotting_parameters)
    _draw_maplink(ax, from_entry, compound, plotting_parameters)
    _draw_maplink(ax, to_entry, compound, plotting_parameters)
