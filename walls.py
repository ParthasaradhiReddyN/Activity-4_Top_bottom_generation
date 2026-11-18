def create_walls(slices,n_walls,W):
    all_outer_walls={}
    all_inner_walls={}
    for z,data in slices.items():

        outer_walls=[]
        inner_walls=[]
        for outer in data["outer"]:
            for i in range(n_walls):
                if outer.is_empty:
                    break
                outer_walls.append(outer)
                outer=outer.buffer(-W)
        for inner in data["inner"]:
            for i in range(n_walls):
                if inner.is_empty:
                    break
                inner_walls.append(inner)
                inner=inner.buffer(W)

        outer_boundaries=[]
        inner_holes=[]
        for outer in data["outer"]:
            outer=outer.buffer(-(n_walls-1)*W)
            outer_boundaries.append(outer)
        for inner in data["inner"]:
            inner=inner.buffer((n_walls-1)*W)
            inner_holes.append(inner)
        
        slices[z]["outer"]=outer_boundaries
        slices[z]["inner"]=inner_holes
        all_outer_walls[z]=outer_walls
        all_inner_walls[z]=inner_walls

    return slices,all_outer_walls,all_inner_walls

