import os
import trimesh
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
import numpy as np

from shapely.geometry import LineString,MultiLineString
from shapely.affinity import rotate
from slice_topology import boundaries_and_holes

from shapely.ops import substring
import heapq
from shapely import Point
import matplotlib.pyplot as plt

import shutil

def list_idx_segment(slice_scanline_indices, infill_segments):
    segments=[]    #store idx and Line strings. if multi split into line string and give same idx for them
    
    for scanline_index, geometry in zip(slice_scanline_indices, infill_segments):

        if geometry.geom_type == "LineString":
            segments.append({
                "scanline": scanline_index,
                "geometry": geometry
            })

        elif geometry.geom_type == "MultiLineString":
            for line in geometry.geoms:
                segments.append({
                    "scanline": scanline_index,
                    "geometry": line
                })
    return segments


def dict_idx_segments(segments):
    remaining = {}
            
    for item in segments:
            if item["scanline"] not in remaining:
                remaining[item["scanline"]]=[item["geometry"]]
            else:
                remaining[item["scanline"]].append(item["geometry"])
    return remaining



def boundary_connector(start, end, final_region, slice_polygon):

    start = Point(start)
    end = Point(end)
   
    # Get all polygons
    if final_region.geom_type == "Polygon":
        polygons = [final_region]
    else:
        polygons = list(final_region.geoms)

    for poly in polygons:

        # Get exterior and interior rings
        rings = [poly.exterior] + list(poly.interiors)

        for ring in rings:

            ring_line = LineString(ring.coords)

            # Both points must lie on the same boundary ring
            if ring_line.distance(start) > 1e-6:
                continue

            if ring_line.distance(end) > 1e-6:
                continue

            # Get distances along the ring
            d1 = ring_line.project(start)
            d2 = ring_line.project(end)

            length = ring_line.length
            lo, hi = sorted((d1, d2))

            # Calculate both possible paths
            path1 = list(substring(ring_line, lo, hi).coords)

            path2 = (
                list(substring(ring_line, hi, length).coords)
                + list(substring(ring_line, 0, lo).coords)[1:]
            )
            # Select the shorter path
            if (hi - lo) <= (length - (hi - lo)):
                coords = path1
                forward = d1 <= d2
            else:
                coords = path2
                forward = d1 >= d2

            if not forward:
                coords.reverse()

            if len(coords) < 2:
                continue

            coords[0] = (start.x, start.y)
            coords[-1] = (end.x, end.y)

            path = LineString(coords)

            if slice_polygon.covers(path):
                return path
    return None




def plot(z, full_region, top_solid_toolpath,bottom_solid_toolpath, normal_toolpath,bounds,all_outer_walls,all_inner_walls):

    fig, ax = plt.subplots(figsize=(10, 10))

    for poly in all_outer_walls[z]:
            x, y = poly.exterior.xy
            ax.plot(x, y, color="black", linewidth=1.2)
    for poly in all_inner_walls[z]:
            x, y = poly.exterior.xy
            ax.plot(x, y, color="black", linewidth=1.2)

    # Plot solid toolpath in red
    for item in top_solid_toolpath:

        geometry = item["geometry"]

        if geometry.is_empty:
            continue

        x, y = geometry.xy

        if item["type"] == "INFILL" or item["type"] == "EXTRUSION":

            ax.plot(
                x, y,
                color="red",
                linewidth=1,
                label="Top solid infill"
            )


        elif item["type"] == "TRAVEL":

            ax.plot(
                x, y,
                color="orange",
                linewidth=1,
                linestyle="--",
                label="Travel"
            )
    # Plot solid toolpath in red
    for item in bottom_solid_toolpath:

        geometry = item["geometry"]

        if geometry.is_empty:
            continue

        x, y = geometry.xy

        if item["type"] == "INFILL" or item["type"] == "EXTRUSION":

            ax.plot(
                x, y,
                color="green",
                linewidth=1,
                label="Bottom solid infill"
            )

        elif item["type"] == "TRAVEL":

            ax.plot(
                x, y,
                color="orange",
                linewidth=1,
                linestyle="--",
                label="Travel"
            )

    # Plot normal toolpath in blue
    for item in normal_toolpath:

        geometry = item["geometry"]

        if geometry.is_empty:
            continue

        x, y = geometry.xy

        if item["type"] == "INFILL" or item["type"] == "EXTRUSION":

            ax.plot(
                x, y,
                color="blue",
                linewidth=1,
                label="Normal infill"
            )

        elif item["type"] == "TRAVEL":

            ax.plot(
                x, y,
                color="orange",
                linewidth=1,
                linestyle="--",
                label="Travel"
            )
  

    # Formatting
    ax.set_aspect("equal", adjustable="box")

    ax.set_title(f"Toolpath Visualization - Z = {z}")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_xlim(bounds[0][0]-1,bounds[1][0]+1)
    ax.set_ylim(bounds[0][1]-1,bounds[1][1]+1)
    


    ax.grid(True, linestyle=":", alpha=0.5)
    

    # Remove duplicate legend entries
    handles, labels = ax.get_legend_handles_labels()

    by_label = dict(zip(labels, handles))

    ax.legend(
        by_label.values(),
        by_label.keys(),
        loc="best"
    )

    plt.tight_layout()

    # Save image
    os.makedirs("toolpath", exist_ok=True)

    plt.savefig(
        f"toolpath/toolpath_Z_{z}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Toolpath image saved: toolpath/toolpath_Z_{z}.png")


def check_connectors(line,all_infill_lines,printed_connectors):
    clear=True
    for infill_line in all_infill_lines:
        if line.intersects(infill_line):
            clear = False
            break

    if clear:
        for previous_connector in printed_connectors:
            if line.intersects(previous_connector):
                clear = False
                break
    return clear



def check_infill_endpoint_collision(connector,all_infill_lines,selected_point,tolerance=1e-5):
    for infill_line in all_infill_lines:

        endpoints = [
            Point(infill_line.coords[0]),
            Point(infill_line.coords[-1])
        ]

        for endpoint in endpoints:

            # Allow the intended destination endpoint
            if endpoint.distance(Point(selected_point)) <= tolerance:
                continue

            # Reject contact with other infill endpoints
            if connector.distance(endpoint) <= tolerance:
                return False

    return True




def connect_infill_segments(
    segments,
    final_region,
    slice_polygon,
    max_connector_length=0,
    line_offset=1e-4,
    obstacle_lines=None
):


    if not segments:
        return []

    remaining = dict_idx_segments(segments)

    all_infill_lines = [
        item["geometry"] for item in segments
    ]

    # Include other infill as obstacles
    if obstacle_lines is not None:
        all_infill_lines += obstacle_lines

    printed_connectors = []

    first_scanline = min(remaining.keys())

    first_segment = remaining[first_scanline].pop(0)

    toolpath = []

    position = first_segment.coords[-1]

    toolpath.append({
        "type": "INFILL",
        "geometry": first_segment
    })

    current_scanline = first_scanline

    while sum(len(v) for v in remaining.values()) > 0:

        selected_segment_index = None
        selected_point = None
        selected_scanline = None
        selected_connector = None
        connector_type = None

        active_scanlines = [
            s for s, v in remaining.items() if v
        ]

        levels = sorted(
            set(abs(s-current_scanline) for s in active_scanlines)
        )

        # Search scanlines by increasing distance
        for diff in levels:

            if diff == 0:
                scanlines_to_check = [current_scanline]
            else:
                scanlines_to_check = [
                    current_scanline-diff,
                    current_scanline+diff
                ]

            nearest_points = []

            for scanline in scanlines_to_check:

                if scanline not in remaining:
                    continue

                for i, segment in enumerate(remaining[scanline]):

                    nearest_points.append(
                        (scanline, i, segment.coords[0])
                    )

                    nearest_points.append(
                        (scanline, i, segment.coords[-1])
                    )

            nearest_points.sort(
                key=lambda x: Point(position).distance(Point(x[2]))
            )

            for scanline, segment_index, point in nearest_points:

                connector = LineString([position, point])

                clear = False
                selected_connector = None
                connector_type = None
                core = None

                # STEP 1: Direct connector
                if (
                    connector.length <= max_connector_length
                    and final_region.buffer(1e-6).covers(connector)
                ):

                    if connector.length > 2*line_offset:

                        core = substring(
                            connector,
                            line_offset,
                            connector.length-line_offset
                        )

                        clear = check_connectors(
                            core,
                            all_infill_lines,
                            printed_connectors
                        )

                    if clear:
                        selected_connector = connector
                        connector_type = "EXTRUSION"

                # STEP 2: Boundary connector
                if not clear:

                    boundary_path = boundary_connector(
                        position,
                        point,
                        final_region,
                        slice_polygon
                    )

                    if boundary_path is not None:

                        if boundary_path.length > 2*line_offset:

                            core = substring(
                                boundary_path,
                                line_offset,
                                boundary_path.length-line_offset
                            )

                            clear = check_connectors(
                                core,
                                all_infill_lines,
                                printed_connectors
                            )

                        if clear:
                            selected_connector = boundary_path
                            connector_type = "EXTRUSION"

                # STEP 3: Check infill endpoint collisions
                if clear:

                    clear = check_infill_endpoint_collision(
                        core,
                        all_infill_lines,
                        point
                    )

                if not clear:
                    continue

                selected_segment_index = segment_index
                selected_point = point
                selected_scanline = scanline

                break

            if selected_segment_index is not None:
                break

        # STEP 4: If no extrusion connector is possible,
        # use a travel connector.
        if selected_segment_index is None:

            nearest_points = []

            for scanline in remaining:

                for i, segment in enumerate(remaining[scanline]):

                    nearest_points.append(
                        (scanline, i, segment.coords[0])
                    )

                    nearest_points.append(
                        (scanline, i, segment.coords[-1])
                    )

            selected_scanline, selected_segment_index, selected_point = min(
                nearest_points,
                key=lambda x: Point(position).distance(Point(x[2]))
            )

            connector = LineString([position, selected_point])

            connector_type = "TRAVEL"

        else:
            connector = selected_connector

        # STEP 5: Select and orient infill segment
        selected_segment = remaining[selected_scanline].pop(
            selected_segment_index
        )

        if selected_point == selected_segment.coords[-1]:

            selected_segment = LineString(
                list(selected_segment.coords)[::-1]
            )

        # Store connector
        if connector_type == "EXTRUSION":
            printed_connectors.append(connector)

        toolpath.append({
            "type": connector_type,
            "geometry": connector
        })

        toolpath.append({
            "type": "INFILL",
            "geometry": selected_segment
        })

        # Update current position
        position = selected_segment.coords[-1]
        current_scanline = selected_scanline

    return toolpath


def prepare_segments(lines):

    segments = []

    for i, line in enumerate(lines):
        segments.append({
            "scanline": i,
            "geometry": line
        })

    return segments


def generate_toolpath(
    infill_lines,
    top_solid_regions,
    bottom_solid_regions,
    normal_regions,
    printable_regions,
    full_regions,bounds,
    all_outer_walls,
    all_inner_walls
):

    shutil.rmtree("toolpath", ignore_errors=True)

    all_toolpaths = {}
  
    for z in sorted(infill_lines.keys()):

        # Get infill lines
        top_solid_lines = infill_lines[z]["top_solid"]
        bottom_solid_lines = infill_lines[z]["bottom_solid"]
        normal_lines = infill_lines[z]["normal"]

        # Prepare segments
        top_solid_segments = prepare_segments(top_solid_lines)
        bottom_solid_segments = prepare_segments(bottom_solid_lines)

        normal_segments = prepare_segments(normal_lines)

        # Collect obstacle lines
        solid_obstacles = [
            item["geometry"] for item in top_solid_segments+bottom_solid_segments
        ]

        normal_obstacles = [
            item["geometry"] for item in normal_segments
        ]

        # Connect solid infill within solid region only
        top_solid_toolpath = connect_infill_segments(
            segments=top_solid_segments,
            final_region=top_solid_regions[z],
            slice_polygon=full_regions[z],
            obstacle_lines=normal_obstacles
        )
        bottom_solid_toolpath = connect_infill_segments(
                    segments=bottom_solid_segments,
                    final_region=bottom_solid_regions[z],
                    slice_polygon=full_regions[z],
                    obstacle_lines=normal_obstacles
                )

        # Connect normal infill within normal region only
        normal_toolpath = connect_infill_segments(
            segments=normal_segments,
            final_region=normal_regions[z],
            slice_polygon=full_regions[z],
            obstacle_lines=solid_obstacles
        )

        # Store toolpaths separately
        all_toolpaths[z] = {
            "top solid": top_solid_toolpath,
            "bottom solid": bottom_solid_toolpath,
            "normal": normal_toolpath
        }

        

        plot(
            z,
            [full_regions[z]],
            top_solid_toolpath,
            bottom_solid_toolpath,
            normal_toolpath,bounds,
            all_outer_walls,
            all_inner_walls
        )
      
    return all_toolpaths

    
