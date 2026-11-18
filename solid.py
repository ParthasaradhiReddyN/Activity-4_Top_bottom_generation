import os
import trimesh
from shapely.geometry import Polygon,LineString
from shapely.affinity import rotate
import numpy as np
from shapely.ops import unary_union

import slice_topology
from solid_toolpath import generate_toolpath
from walls import create_walls
from infill import generate_infill_segments,full_printable_regions

def find_surface_triangles(mesh,slices):
    top_triangles_by_layer={}
    bottom_triangles_by_layer={}
    face_normals=mesh.face_normals
    for i,face in enumerate(mesh.faces):
        normal=face_normals[i]
        if normal[2]>0.999:
            triangle_vertices=mesh.vertices[face]
            triangle_z=np.mean(triangle_vertices[:,2])
            closest_z=min(slices.keys(),key=lambda z: abs(z-triangle_z))
            if closest_z not in top_triangles_by_layer:
                top_triangles_by_layer[closest_z]=[]
            top_triangles_by_layer[closest_z].append(i)
  
        if normal[2]<-0.999:
            triangle_vertices=mesh.vertices[face]
            triangle_z=np.mean(triangle_vertices[:,2])
            closest_z=min(slices.keys(),key=lambda z: abs(z-triangle_z))
            if closest_z not in bottom_triangles_by_layer:
                bottom_triangles_by_layer[closest_z]=[]
            bottom_triangles_by_layer[closest_z].append(i)
    return top_triangles_by_layer,bottom_triangles_by_layer

def find_surface_regions(mesh,slices,top_triangles_by_layer,bottom_triangles_by_layer):
    top_regions={}
    bottom_regions={}
    for z in slices:
        top_polygons=[]
        for index in top_triangles_by_layer.get(z,[]):
            face=mesh.faces[index]
            vertices=mesh.vertices[face]
            triangle_polygon=Polygon(vertices[:,:2])
            if triangle_polygon.is_valid and triangle_polygon.area>0:
                top_polygons.append(triangle_polygon)
        if top_polygons:
            top_regions[z]=unary_union(top_polygons)

        else:
            top_regions[z]=Polygon()

        bottom_polygons=[]
        for index in bottom_triangles_by_layer.get(z,[]):
            face=mesh.faces[index]
            vertices=mesh.vertices[face]
            triangle_polygon=Polygon(vertices[:,:2])
            if triangle_polygon.is_valid and triangle_polygon.area>0:
                bottom_polygons.append(triangle_polygon)
        if bottom_polygons:
            bottom_regions[z]=unary_union(bottom_polygons)

        else:
            bottom_regions[z]=Polygon()
    return top_regions,bottom_regions

def top_bottom_regions(top_regions,bottom_regions,printable_regions,slices,n_layers=3):
    
    top_solid_regions={z:Polygon() for z in slices}
    bottom_solid_regions={z:Polygon() for z in slices}

    z_values=sorted(slices.keys())
    for i,z in enumerate(z_values):
        buffered_top=top_regions[z]
        valid_layers=[]
        for k in range(n_layers):
            target_index=i-k
            if target_index<0:
                continue
            target_z=z_values[target_index]
            projected_region=buffered_top.intersection(printable_regions[target_z])
            top_solid_regions[target_z]=(top_solid_regions[target_z].union(projected_region))
            if not projected_region.is_empty:
                valid_layers.append(target_index)
        if len(valid_layers) < n_layers:
                    extra_index = i-n_layers
                    if extra_index >= 0:
                        extra_z = z_values[extra_index]
                        extra_region = top_regions[z].intersection(
                            printable_regions[extra_z]
                        )
                        if not extra_region.is_empty:
                            top_solid_regions[extra_z] = (
                                top_solid_regions[extra_z].union(extra_region)
                            )
    
    for i,z in enumerate(z_values):
        buffered_bottom=bottom_regions[z]
        valid_layers = []
        for k in range(n_layers):
            target_index=i+k
            if target_index>=len(z_values):
                continue
            target_z=z_values[target_index]
            projected_region=buffered_bottom.intersection(printable_regions[target_z])
            bottom_solid_regions[target_z]=(bottom_solid_regions[target_z].union(projected_region))
            if not projected_region.is_empty:
                            valid_layers.append(target_index)
        if len(valid_layers) < n_layers:
                    extra_index = i+n_layers
                    if extra_index < len(z_values):
                        extra_z = z_values[extra_index]
                        extra_region = bottom_regions[z].intersection(
                            printable_regions[extra_z]
                        )
                        if not extra_region.is_empty:
                            bottom_solid_regions[extra_z] = (
                                bottom_solid_regions[extra_z].union(extra_region)
                            )
        
    return top_solid_regions,bottom_solid_regions

def seperate_normal_region(top_solid_regions,bottom_solid_regions,printable_regions,W):
    solid_gap=W
    normal_regions={}
    z_values=sorted(top_solid_regions.keys())
    for z in z_values:
        original_solid=(top_solid_regions[z].union(bottom_solid_regions[z])).intersection(printable_regions[z])
        normal_regions[z]=printable_regions[z].difference(original_solid.buffer(solid_gap))
    return normal_regions



def main(stl_path="models/cube.stl", W=0.2, D=0.25,line_angle=30,n_layers=3,n_walls=1):

    mesh=trimesh.load_mesh(stl_path)
    solid_line_spacing=W
    line_spacing=W/D

    extremes=mesh.bounds

    slices=slice_topology.boundaries_and_holes(stl_path ,z_increment=0.4,threshold_area=0)

    slices,all_outer_walls,all_inner_walls=create_walls(slices,n_walls,W)

    top_triangles_by_layer,bottom_triangles_by_layer=find_surface_triangles(mesh,slices)

    top_regions,bottom_regions=find_surface_regions(mesh,slices,top_triangles_by_layer,bottom_triangles_by_layer)

    full_regions,printable_regions=full_printable_regions(slices,W)

    top_solid_regions,bottom_solid_regions=top_bottom_regions(top_regions,bottom_regions,printable_regions,slices,n_layers)

    normal_regions=seperate_normal_region(top_solid_regions,bottom_solid_regions,printable_regions,W)

    infill_lines=generate_infill_segments(top_solid_regions,bottom_solid_regions,normal_regions,line_angle,solid_line_spacing,line_spacing)

    all_toolpaths = generate_toolpath(
        infill_lines,
        top_solid_regions,
        bottom_solid_regions,
        normal_regions,
        printable_regions,
        full_regions,
        extremes,
        all_outer_walls,
        all_inner_walls
    )
    return all_toolpaths

print(main())













