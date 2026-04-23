import matplotlib
import numpy
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize
from matplotlib.patches import FancyArrowPatch, ArrowStyle, PathPatch

import metaboglobe
from metaboglobe._util import wrap_text
from metaboglobe.kegg_pathway import KeggMap, KeggRelation, KeggReaction, RelationType, ReactionType, EntryType, \
    KeggEntry
from metaboglobe.plotting import PlotStyle, MPLColor
from metaboglobe.plotting._collision_map import CollisionMap, TextWithAnchor
from metaboglobe.plotting._curve_2d import Curve2
from metaboglobe.plotting._vector_2d import Vector2


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


def plot_kegg(ax: Axes, kegg_map: KeggMap, plot_style: PlotStyle) -> ScalarMappable:
    """Draws the KEGG map, with double arrows for revisble/two-way-irrervisble reactions. Returns a mappable for use
    in figure.colorbar(...).
    """

    # Set up plot
    ax.set_facecolor(plot_style.facecolor)
    _adjust_limits(ax, kegg_map)
    if plot_style.hide_ticks_and_spines:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

    # Draw entries
    for entry in kegg_map.entries:
        if plot_style.plot_entries_without_reactions or kegg_map.has_relations_or_reactions(entry):
            _draw_entry(ax, entry, plot_style)

    # Draw relations and reactions
    for relation in kegg_map.relations:
        if relation.relation_type == RelationType.MAPLINK:
            _draw_maplinks(ax, kegg_map, relation, plot_style)
    for reaction in kegg_map.reactions:
        _draw_reaction(ax, kegg_map, reaction, plot_style)

    # Draw compound names
    _draw_compound_names(ax, kegg_map, plot_style)

    return ScalarMappable(Normalize(vmin=plot_style.flux_vmin, vmax=plot_style.flux_vmax), plot_style.flux_cmap)


def _draw_compound_names(ax: Axes, kegg_map: KeggMap, plot_style: PlotStyle):
    collision_map = CollisionMap(ax)
    for artist in ax.get_children():
        collision_map.add_artist(artist)

    texts = list()
    for entry in kegg_map.entries:
        if entry.entry_type == EntryType.COMPOUND:
            if plot_style.plot_entries_without_reactions or kegg_map.has_relations_or_reactions(entry):
                display_name = metaboglobe.kegg_pathway.get_display_name(entry.name)
                texts.append(TextWithAnchor(text=display_name, x=entry.x, y=entry.y))
    collision_map.fit_text(ax, texts, fontsize=6, zorder=10)


def _draw_entry(ax: Axes, entry: KeggEntry, plot_style: PlotStyle):
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
        ax.add_patch(matplotlib.patches.Circle((entry.x, entry.y), plot_style.compound_radius,
                                               facecolor=plot_style.compound_nan_color,
                                               edgecolor=plot_style.compound_edgecolor,
                                               linewidth=plot_style.compound_linewidth))
        return
    elif self == EntryType.ORTHOLOG:
        return  # Not interested in orthologs, they just clutter the figure
    elif self == EntryType.GENE:  # Genes
        if "," in entry.name:
            display_name = entry.name.split(",")[0] + ", ..."
        else:
            display_name = entry.name
        bbox = {"facecolor": plot_style.enzyme_facecolor, "edgecolor": plot_style.enzyme_edgecolor,
                "linewidth": plot_style.enzyme_linewidth}
        if plot_style.enzyme_rounding:
            bbox["boxstyle"] = f"round,pad={plot_style.enzyme_padding}"
        else:
            bbox["pad"] = plot_style.enzyme_padding
        ax.text(entry.x, entry.y, display_name, ha="center", va="center", fontsize=6, zorder=10,
                color=plot_style.enzyme_textcolor, bbox=bbox)
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


def _draw_reaction(ax: Axes, kegg_map: KeggMap, reaction: KeggReaction, plot_style: PlotStyle):
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

    vmin = plot_style.flux_vmin
    vspread = plot_style.flux_vmax - plot_style.flux_vmin

    forward_value = kegg_map.forward_value(reaction)
    forward_color = plot_style.flux_nan_color if numpy.isnan(forward_value) else plot_style.flux_cmap((forward_value - vmin) / vspread)
    forward_width = plot_style.flux_nan_linewidth if numpy.isnan(forward_value) else plot_style.flux_linewidth

    curve_forward, curve_backward = curve.split()
    ax.add_patch(FancyArrowPatch(path=curve_forward.to_path(), arrowstyle=ArrowStyle("-|>",
                 head_length=plot_style.flux_arrowsize, head_width=plot_style.flux_arrowsize / 2),
                 color=forward_color, linewidth=forward_width, joinstyle=plot_style.flux_joinstyle, capstyle=plot_style.flux_capstyle))

    if reaction.reaction_type == ReactionType.REVERSIBLE:
        # Also draw backwards arrow
        backward_value = kegg_map.backward_value(reaction)
        backward_color = plot_style.flux_nan_color if numpy.isnan(backward_value) else plot_style.flux_cmap((backward_value - vmin) / vspread)
        backward_width = plot_style.flux_nan_linewidth if numpy.isnan(backward_value) else plot_style.flux_linewidth
        ax.add_patch(FancyArrowPatch(path=curve_backward.to_path(), arrowstyle=ArrowStyle("-|>",
                     head_length=plot_style.flux_arrowsize, head_width=plot_style.flux_arrowsize / 2),
                     color=backward_color, linewidth=backward_width, joinstyle=plot_style.flux_joinstyle, capstyle=plot_style.flux_capstyle))

def _draw_arrow(ax: Axes, x1: float, y1: float, x2: float, y2: float, cmap: Colormap, value: float, nan_color: MPLColor, linewidth: float, *, arrowstyle: str = "->") -> None:
    if numpy.isnan(value):
        linewidth /= 2
        color = nan_color
    else:
        color = cmap(value)

    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(
        arrowstyle=arrowstyle, color=color, linewidth=linewidth))


def _draw_maplink(ax: Axes, entry1: KeggEntry, entry2: KeggEntry, plot_style: PlotStyle) -> None:
    """Draws a link between a reference to another pathway, and a compound. Does nothing if other entries are provided."""
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

    curve = Curve2(Vector2(compound_entry.x, compound_entry.y))

    if map_min_x <= compound_entry.x <= map_max_x:
        # Easy, we can draw a straight vertical line
        to_y = map_entry.y - map_entry.height / 2 if compound_entry.y < map_entry.y else map_entry.y + map_entry.height / 2
        curve.append_line_to(Vector2(compound_entry.x, to_y))
    elif map_min_y <= compound_entry.y <= map_max_y:
        # Also easy, we can draw a straight horizontal line
        to_x = map_entry.x - map_entry.width / 2 if compound_entry.x < map_entry.x else map_entry.x + map_entry.width / 2
        curve.append_line_to(Vector2(to_x, compound_entry.y))
    else:
        # Define all possible paths (always first attachment to pathway rectangle, then the corner)
        possible_paths = [[Vector2(map_min_x, map_entry.y), Vector2(compound_entry.x, map_entry.y)],
                          [Vector2(map_max_x, map_entry.y), Vector2(compound_entry.x, map_entry.y)],
                          [Vector2(map_entry.x, map_min_y), Vector2(map_entry.x, compound_entry.y)],
                          [Vector2(map_entry.x, map_max_y), Vector2(map_entry.x, compound_entry.y)]]

        # Find the shortest path
        shortest_distance_squared = float("inf")
        shortest_distance_path = None
        for attachment, path_corner in possible_paths:
            distance_squared = (attachment.x - compound_entry.x) ** 2 + (attachment.y - compound_entry.y) ** 2
            if distance_squared < shortest_distance_squared:
                shortest_distance_squared = distance_squared
                shortest_distance_path = attachment, path_corner

        # Build that path
        curve.append_cut_corner_to(corner=shortest_distance_path[1], end=shortest_distance_path[0])
    patch = PathPatch(curve.to_path(), lw=plot_style.maplink_linewidth, edgecolor=plot_style.maplink_edgecolor,
                      linestyle=plot_style.maplink_linestyle, zorder=-5, fill=False)
    ax.add_patch(patch)


def _draw_maplinks(ax: Axes, kegg_map: KeggMap, relation: KeggRelation, plot_style: PlotStyle) -> None:
    # We don't draw map links from/to genes, to avoid clutter - so only towards compounds
    from_entry = kegg_map.entry(relation.from_id)
    to_entry = kegg_map.entry(relation.to_id)
    compound = kegg_map.entry(relation.compound_id)

    _draw_maplink(ax, from_entry, to_entry, plot_style)
    _draw_maplink(ax, from_entry, compound, plot_style)
    _draw_maplink(ax, to_entry, compound, plot_style)
