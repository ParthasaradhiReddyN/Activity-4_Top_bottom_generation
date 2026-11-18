from shapely.geometry import Polygon,LineString
from shapely.affinity import rotate
import numpy as np

def generate_infill_lines(region,spacing,angle=0):
    if region.is_empty:
        return []
    min_x,min_y,max_x,max_y=region.bounds
    centre_x=(min_x+max_x)/2
    centre_y=(min_y+max_y)/2
    diagonal=np.sqrt((max_x-min_x)**2+(max_y-min_y)**2)
    lines=[]
    y_values=np.arange(centre_y-diagonal,centre_y+diagonal,spacing)
    for y in y_values:
        line=LineString([(centre_x-diagonal,y),(centre_x+diagonal,y)])
        line=rotate(line,angle,origin=(centre_x,centre_y))
        clipped=line.intersection(region)
        if clipped.is_empty:
            continue
        if clipped.geom_type=="LineString":
            lines.append(clipped)
        elif clipped.geom_type=="MultiLineString":
            for segment in clipped.geoms:
                if segment.length>0:
                    lines.append(segment)
        elif clipped.geom_type=="GeometryCollection":
            for geom in clipped.geoms:
                if geom.geom_type=="LineString" and geom.length>0:
                    lines.append(geom)
    return lines

def full_printable_regions(slices,W):
    buffer_offset=W
    full_regions={}
    printable_regions={}
    for z,data in slices.items():
        full_region=Polygon()
        printable_region=Polygon()
        all_polygons=data["outer"]+data["inner"]
        all_polygons=sorted(all_polygons,key=lambda poly:poly.area,reverse=True)
        for poly in all_polygons:
            if poly in data["outer"]:
                buffered_poly = poly.buffer(-buffer_offset)
                full_region=full_region.union(poly)
                printable_region=printable_region.union(buffered_poly)
            else:
                buffered_poly = poly.buffer(buffer_offset)
                full_region=full_region.difference(poly)
                printable_region=printable_region.difference(buffered_poly)
        full_regions[z]=full_region
        printable_regions[z]=printable_region
    return full_regions,printable_regions


def generate_infill_segments(top_solid_regions,bottom_solid_regions,normal_regions,line_angle,solid_line_spacing,normal_line_spacing):
    infill_lines={}
    z_values=sorted(top_solid_regions.keys())
    for i,z in enumerate(z_values):
        if i%2!=0:
            angle=line_angle
        else:
            angle=90+line_angle

        top_solid_lines=generate_infill_lines(top_solid_regions[z],spacing=solid_line_spacing,angle=angle)
        bottom_solid_lines=generate_infill_lines(bottom_solid_regions[z],spacing=solid_line_spacing,angle=angle)
        normal_lines=generate_infill_lines(normal_regions[z],spacing=normal_line_spacing,angle=angle)
        infill_lines[z]={"top_solid":top_solid_lines,"bottom_solid":bottom_solid_lines,"normal":normal_lines}
    return infill_lines