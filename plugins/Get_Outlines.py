# import matplotlib.pyplot as plt
# from matplotlib.patches import Rectangle, Ellipse, Polygon, PathPatch
# from matplotlib.transforms import Affine2D
# from matplotlib.path import Path
# from matplotlib.collections import PatchCollection
# import numpy as np
import shapely
import shapely.ops

def Get_Outlines(data, choosen_net):
    # figure, axes = plt.subplots()
    # axes.set_aspect(1)
    # axes.invert_yaxis()

    shape = {0: "Kreis", 1: "Oval", 2: "Rechteck"}

    Color = {0: "red", 1: "green", 2: "orange", 3: "cyan", 4: "pink", 31: "blue"}
    for i in range(0, 32):
        if i not in Color:
            Color[i] = "silver"
    
    shapes = [[] for _ in range(32)]
    drills = []

    for uuid, d in list(data.items()):
        if d["Netname"] != choosen_net:
            continue
        if d["type"] == "VIA":
            # add circle to shapes
            for layer in d["Layer"]:
                shapes[layer].append(shapely.geometry.Point(d["Position"]).buffer(d["Width"] / 2))
            # circ = plt.Circle(d["Position"], d["Width"] / 2, color="grey", alpha=0.5)
            # axes.add_artist(circ)
            if "Drill" in d and len(d["Layer"]) > 1: 
                # axes.add_artist(plt.Circle(d["Position"], d["Drill"] / 2, color="w"))
                drills.append({ "Center": d["Position"], "Width": d["Drill"], "Layer": d["Layer"]})
            # plt.plot(
            #     [d["Position"][0] - 0.5, d["Position"][0] + 0.5],
            #     [d["Position"][1], d["Position"][1]],
            #     color="grey",
            #     alpha=0.5,
            # )
            # plt.plot(
            #     [d["Position"][0], d["Position"][0]],
            #     [d["Position"][1] - 0.5, d["Position"][1] + 0.5],
            #     color="grey",
            #     alpha=0.5,
            # )
        elif d["type"] == "WIRE":
            for layer in d["Layer"]:
                shapes[layer].append(shapely.geometry.LineString([d["Start"], d["End"]]).buffer(d["Width"] / 2))
            # plt.arrow(
            #     d["Start"][0],
            #     d["Start"][1],
            #     d["End"][0] - d["Start"][0],
            #     d["End"][1] - d["Start"][1],
            #     width=d["Width"],
            #     head_length=0,
            #     head_width=d["Width"],
            #     color=Color[layer],
            #     alpha=0.5,
            # )
            # axes.add_artist(plt.Circle(d["Start"], d["Width"] / 2, color=Color[layer], alpha=0.25))
            # axes.add_artist(plt.Circle(d["End"], d["Width"] / 2, color=Color[layer], alpha=0.25))
            # data[uuid]["R"] = calcResWIRE(d["Start"], d["End"], d["Width"])
        elif d["type"] == "PAD":
            for l in d["Layer"]:
                if d["Shape"] in {0, 2}:  # oval
                    # ellip = Ellipse(
                    #     d["Position"],
                    #     *d["Size"],
                    #     color=Color[l],
                    #     alpha=0.5,
                    #     angle=d["Orientation"],
                    # )
                    # axes.add_patch(ellip)
                    shapes[l].append(shapely.geometry.Point(d["Position"]).buffer(1))
                    shapes[l][-1] = shapely.affinity.scale(shapes[l][-1], d["Size"][0] / 2, d["Size"][1] / 2)
                    shapes[l][-1] = shapely.affinity.rotate(shapes[l][-1], d["Orientation"], d["Position"])
                else:
                    # rec = plt.Rectangle(
                    #     np.array(d["Position"]) - np.array(d["Size"]) / 2,
                    #     width=d["Size"][0],
                    #     height=d["Size"][1],
                    #     color=Color[l],
                    #     alpha=0.5,
                    #     transform=Affine2D().rotate_deg_around(
                    #         *d["Position"], d["Orientation"]
                    #     )
                    #     + axes.transData,
                    # )
                    # axes.add_patch(rec)
                    # draw rectangle with shapely
                    line = shapely.geometry.LineString([
                        (d["Position"][0] - d["Size"][0] / 2, d["Position"][1] - d["Size"][1] / 2),
                        (d["Position"][0] + d["Size"][0] / 2, d["Position"][1] - d["Size"][1] / 2),
                        (d["Position"][0] + d["Size"][0] / 2, d["Position"][1] + d["Size"][1] / 2),
                        (d["Position"][0] - d["Size"][0] / 2, d["Position"][1] + d["Size"][1] / 2),
                        (d["Position"][0] - d["Size"][0] / 2, d["Position"][1] - d["Size"][1] / 2),
                    ])
                    line_r = shapely.affinity.rotate(line, d["Orientation"], d["Position"])
                    # convert to polygon
                    shapes[l].append(shapely.geometry.Polygon(line_r))
            if "Drill" in d and len(d["Layer"]) > 1: 
                # axes.add_artist(plt.Circle(d["Position"], d["Drill"] / 2, color="w"))
                drills.append({ "Center": d["Position"], "Width": d["Drill"], "Layer": d["Layer"]})
        elif d["type"] == "ZONE":
            polys = d["Polygons"]
            for layer in polys:
                for poly in polys[layer]:
                    shapes[layer].append(shapely.geometry.Polygon(poly))

    # merge the polygons of each layer
    merged_shapes = []
    for layer in shapes:
        if len(layer) == 0:
            merged_shapes.append(None)
        else:
            merged_shapes.append(shapely.ops.unary_union(layer))
    
    # def plot_polygon(ax, poly, **kwargs):
    #     # plot representative point
    #     point = poly.representative_point()
    #     ax.plot(point.x, point.y, "o", color=kwargs["color"], alpha=0.5)
    #     path = Path.make_compound_path(
    #         Path(np.asarray(poly.exterior.coords)[:, :2]),
    #         *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])
    #     patch = PathPatch(path, **kwargs)
    #     collection = PatchCollection([patch], **kwargs)
    #     ax.add_collection(collection, autolim=True)
    #     ax.autoscale_view()
    #     return collection
    # # plot the merged polygons
    # for layer, shape in enumerate(merged_shapes):
    #     if shape is None:
    #         continue
    #     if isinstance(shape, shapely.geometry.Polygon):
    #         plot_polygon(axes, shape, color=Color[layer], alpha=0.5)
    #     elif isinstance(shape, shapely.geometry.MultiPolygon):
    #         for poly in shape.geoms:
    #             plot_polygon(axes, poly, color=Color[layer], alpha=0.5)
    #     else:
    #         print("unknown shape", type(shape))
    # # plot drills like black crosses
    # for drill in drills:
    #     for layer in drill["Layer"]:
    #         axes.plot(
    #             [drill["Center"][0] - 0.5, drill["Center"][0] + 0.5],
    #             [drill["Center"][1], drill["Center"][1]],
    #             color=Color[layer],
    #             alpha=0.5,
    #         )
    #         axes.plot(
    #             [drill["Center"][0], drill["Center"][0]],
    #             [drill["Center"][1] - 0.5, drill["Center"][1] + 0.5],
    #             color=Color[layer],
    #             alpha=0.5,
    #         )
    # plt.show()
            
    return merged_shapes, drills