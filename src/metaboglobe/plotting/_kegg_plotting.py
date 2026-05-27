import matplotlib
import numpy
from matplotlib import patheffects
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch, ArrowStyle, PathPatch

import metaboglobe
from metaboglobe._util import wrap_text
from metaboglobe.kegg_pathway import KeggMap, KeggReactionArrow, ReactionType, KeggCompoundInMap, \
    KeggReactionEnzymeInMap, KeggOtherPathwayRelationInMap, KeggPathwayReferenceInMap
from metaboglobe.plotting import PlotStyle
from metaboglobe.math.box_2d import Box2
from metaboglobe.plotting._collision_map import CollisionMap, TextWithAnchor
from metaboglobe.plotting._curve_2d import Curve2
from metaboglobe.math.vector_2d import Vector2


def _adjust_limits(ax: Axes, kegg_map: KeggMap):
    """Adjusts the xlim and ylim of the plot to fit the map."""
    box = kegg_map.box()

    # Adjust the limits
    ax.set_xlim(box.min.x, box.max.x)
    ax.set_ylim(box.max.y, box.min.y)  # Reversed, so that low y is on top
    ax.set_aspect("equal")


def plot_kegg(ax: Axes, kegg_map: KeggMap, plot_style: PlotStyle) -> ScalarMappable:
    """Draws the KEGG map, with double arrows for revisble/two-way-irrervisble reactions. Returns a mappable for use
    in figure.colorbar(...).
    """

    # Set up plot
    ax.set_title(kegg_map.title)
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
    for compound in kegg_map.compounds:
        if plot_style.plot_entries_without_reactions or kegg_map.has_relations_or_reactions(compound.compound_id):
            _draw_compound(ax, kegg_map, compound, plot_style)
    for pathway_reference in kegg_map.other_pathways:
        if plot_style.plot_entries_without_reactions or kegg_map.has_relations_or_reactions(pathway_reference.pathway_id):
            _draw_pathway_reference(ax, pathway_reference, plot_style)
    for enzyme in kegg_map.reaction_enzymes:
        if plot_style.plot_entries_without_reactions or kegg_map.has_relations_or_reactions(enzyme.reaction_id):
            _draw_enzyme(ax, enzyme, plot_style)

    # Draw relations and reactions
    for relation in kegg_map.other_pathway_relations:
        _draw_maplink(ax, relation, plot_style)
    for reaction in kegg_map.reaction_arrows:
        _draw_reaction(ax, kegg_map, reaction, plot_style)

    # Draw compound names
    _draw_compound_names(ax, kegg_map, plot_style)

    return plot_style.flux_mappable()


def _draw_compound_names(ax: Axes, kegg_map: KeggMap, plot_style: PlotStyle):
    collision_map = CollisionMap(ax)
    for artist in ax.get_children():
        collision_map.add_artist(artist)

    texts = list()
    for entry in kegg_map.compounds:
        if plot_style.plot_entries_without_reactions or kegg_map.has_relations_or_reactions(entry.compound_id):
            display_name = metaboglobe.kegg_pathway.get_display_name(entry.compound_id)
            texts.append(TextWithAnchor(text=display_name, x=entry.x, y=entry.y))
    collision_map.fit_text(ax, texts, fontsize=6, zorder=10)


def _draw_compound(ax: Axes, kegg_map: KeggMap, entry: KeggCompoundInMap, plot_style: PlotStyle):
    value = kegg_map.get_compound_score(entry.compound_id)
    vmin = plot_style.compound_vmin
    vspread = plot_style.compound_vmax - plot_style.compound_vmin
    color = plot_style.compound_nan_color if numpy.isnan(value) else plot_style.flux_cmap((value - vmin) / vspread)
    ax.add_patch(matplotlib.patches.Circle((entry.x, entry.y), plot_style.compound_radius,
                                           facecolor=color,
                                           edgecolor=plot_style.compound_edgecolor,
                                           linewidth=plot_style.compound_linewidth))


def _draw_pathway_reference(ax: Axes, entry: KeggPathwayReferenceInMap, plot_style: PlotStyle):
    display_name = wrap_text(entry.name,
                             int(entry.width // 6))  # Wrap text to fit within the entry box, assuming an average character width of 6 pixels
    rect = matplotlib.patches.Rectangle((entry.x - entry.width / 2, entry.y - entry.height / 2),
                                        entry.width, entry.height, fill=True,
                                        linewidth=plot_style.pathway_reference_linewidth,
                                        edgecolor=plot_style.pathway_reference_edgecolor,
                                        facecolor=plot_style.pathway_reference_facecolor, zorder=-10)
    ax.add_patch(rect)
    ax.text(entry.x, entry.y, display_name, ha="center", va="center", fontsize=6, zorder=10)


def _draw_enzyme(ax: Axes, entry: KeggReactionEnzymeInMap, plot_style: PlotStyle):
    bbox = {"facecolor": plot_style.enzyme_facecolor, "edgecolor": plot_style.enzyme_edgecolor,
            "linewidth": plot_style.enzyme_linewidth}
    if plot_style.enzyme_rounding:
        bbox["boxstyle"] = f"round,pad={plot_style.enzyme_padding}"
    else:
        bbox["pad"] = plot_style.enzyme_padding
    ax.text(entry.x, entry.y, entry.display_name, ha="center", va="center", fontsize=6, zorder=10,
            color=plot_style.enzyme_textcolor, bbox=bbox)


def _snap_to_box(entry: Vector2, box: Box2, *, margin_px: float = 2) -> Vector2:
    """If a point is within margin_px from the edge of the box, it is moved to the edge of the box."""
    if abs(entry.x - box.min.x) < margin_px:
        entry = entry.with_x(box.min.x)
    elif abs(entry.x - box.max.x) < margin_px:
        entry = entry.with_x(box.max.x)
    if abs(entry.y - box.min.y) < margin_px:
        entry = entry.with_y(box.min.y)
    elif abs(entry.y - box.max.y) < margin_px:
        entry = entry.with_y(box.max.y)
    return entry


def _to_vector(entry: KeggCompoundInMap | KeggReactionEnzymeInMap) -> Vector2:
    return Vector2(entry.x, entry.y)


def _draw_reaction(ax: Axes, kegg_map: KeggMap, reaction: KeggReactionArrow, plot_style: PlotStyle):
    """Draws the reaction arrow for a single reaction."""

    curve_forward = _find_reaction_curve(reaction, plot_style)

    vmin = plot_style.flux_vmin
    vspread = plot_style.flux_vmax - plot_style.flux_vmin

    forward_value = kegg_map.get_reaction_score(reaction.reaction_id.forwards())
    forward_color = plot_style.flux_nan_color if numpy.isnan(forward_value) else plot_style.flux_cmap((forward_value - vmin) / vspread)
    forward_arrowwidth = plot_style.flux_nan_arrowwidth if numpy.isnan(forward_value) else plot_style.flux_arrowwidth
    forward_edgewidth = plot_style.flux_nan_edgewidth if numpy.isnan(forward_value) else plot_style.flux_edgewidth
    forward_edgecolor = plot_style.flux_edgecolor if numpy.isnan(forward_value) else plot_style.flux_edgecolor
    forward_zorder = 0 if numpy.isnan(forward_value) else 1

    if plot_style.plot_double_arrows:
        arrow_separation = plot_style.flux_nan_separation if numpy.isnan(forward_value) else plot_style.flux_separation

        # Make space for two arrows, draw which ones are appropriate
        curve_forward, curve_backward = curve_forward.split(separation_distance=arrow_separation)
        ax.add_patch(FancyArrowPatch(path=curve_forward.to_path(), arrowstyle=ArrowStyle("-|>",
                     head_length=plot_style.flux_arrowsize, head_width=plot_style.flux_arrowsize / 2),
                     color=forward_color, linewidth=forward_arrowwidth, joinstyle=plot_style.flux_joinstyle,
                     capstyle=plot_style.flux_capstyle, zorder=forward_zorder,
                     path_effects=[patheffects.withStroke(linewidth=forward_arrowwidth + forward_edgewidth * 2, foreground=forward_edgecolor)]))

        if reaction.reaction_type == ReactionType.REVERSIBLE:
            # Also draw backwards arrow
            backward_value = kegg_map.get_reaction_score(reaction.reaction_id.backwards())
            backward_color = plot_style.flux_nan_color if numpy.isnan(backward_value) else plot_style.flux_cmap((backward_value - vmin) / vspread)
            backward_arrowwidth = plot_style.flux_nan_arrowwidth if numpy.isnan(backward_value) else plot_style.flux_arrowwidth
            backward_edgewidth = plot_style.flux_nan_edgewidth if numpy.isnan(backward_value) else plot_style.flux_edgewidth
            backward_edgecolor = plot_style.flux_nan_edgecolor if numpy.isnan(backward_value) else plot_style.flux_edgecolor
            backward_zorder = 0 if numpy.isnan(backward_value) else 1
            ax.add_patch(FancyArrowPatch(path=curve_backward.to_path(), arrowstyle=ArrowStyle("-|>",
                         head_length=plot_style.flux_arrowsize,
                         head_width=plot_style.flux_arrowsize / 2),
                         color=backward_color, linewidth=backward_arrowwidth,
                         joinstyle=plot_style.flux_joinstyle,
                         capstyle=plot_style.flux_capstyle, zorder=backward_zorder,
                         path_effects=[patheffects.withStroke(linewidth=backward_arrowwidth + backward_edgewidth * 2, foreground=backward_edgecolor)]))

    else:
        # Just draw a single arrow
        stylename = "<|-|>" if reaction.reaction_type == ReactionType.REVERSIBLE else "-|>"
        ax.add_patch(FancyArrowPatch(path=curve_forward.to_path(), arrowstyle=ArrowStyle(stylename,
                     head_length=plot_style.flux_arrowsize, head_width=plot_style.flux_arrowsize / 2),
                     color=forward_color, linewidth=forward_arrowwidth, joinstyle=plot_style.flux_joinstyle,
                     capstyle=plot_style.flux_capstyle, zorder=forward_zorder,
                     path_effects=[patheffects.withStroke(linewidth=forward_arrowwidth + forward_edgewidth * 2, foreground=forward_edgecolor)]))


def _find_reaction_curve(reaction: KeggReactionArrow, plot_style: PlotStyle) -> Curve2:
    # Build an initial box
    from_entry = _to_vector(reaction.substrate)
    to_entry = _to_vector(reaction.product)
    enzyme_entry = _to_vector(reaction.reaction_enzyme_in_map)
    box = Box2.enclosing(from_entry, to_entry, enzyme_entry)

    # If any of the points are very close to the box, just move them there
    # (KEGG isn't very precise with aligning things)
    from_entry = _snap_to_box(from_entry, box)
    to_entry = _snap_to_box(to_entry, box)
    enzyme_entry = _snap_to_box(enzyme_entry, box)
    box = Box2.enclosing(from_entry, to_entry, enzyme_entry)

    # Check if the enzyme is on the side of a box (necessary for determining the corners in our path)
    enzyme_on_min_or_max_x = enzyme_entry.x == box.min.x or enzyme_entry.x == box.max.x
    enzyme_on_min_or_max_y = enzyme_entry.y == box.min.y or enzyme_entry.y == box.max.y

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
        else:  # Some weird diagonal situation, start with straight line, then rounded line
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
        else:  # Some weird diagonal situation, start with rounded line, then straight line
            curve.append_curve_then_line_to(to_entry)

    # Make space for the compound itself, some margin, and the arrow head
    curve.shorten_both_sides(plot_style.compound_radius * 1.4 + plot_style.flux_arrowsize * 1.2)

    # Arrowhead is part of the line, so append a straight line segment to make space
    curve.extend_both_sides(plot_style.flux_arrowsize * 1.2)

    return curve


def _draw_maplink(ax: Axes, relation: KeggOtherPathwayRelationInMap, plot_style: PlotStyle) -> None:
    """Draws a link between a reference to another pathway, and a compound."""
    compound_entry = relation.compound
    map_entry = relation.pathway

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
    patch = PathPatch(curve.to_path(), lw=plot_style.pathway_link_linewidth, edgecolor=plot_style.pathway_link_edgecolor,
                      linestyle=plot_style.pathway_link_linestyle, zorder=-5, fill=False)
    ax.add_patch(patch)
