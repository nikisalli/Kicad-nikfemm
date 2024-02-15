import pcbnew
from importlib import reload
import sys
import os
import wx
import traceback
import numpy as np

# import pip
# def install(package):
#     if hasattr(pip, "main"):
#         pip.main(["install", package])
#     else:
#         pip._internal.main(["install", package])
# install("PySpice")
# import PySpice

debug = 0

def plot_polygon(poly, layer, simulation):
    # find a point that is inside the polygon
    point = poly.representative_point()
    # convert to a list
    point = [point.x, point.y]

    simulation.draw_region(point, 5.95e7, layer)

    # get exterior and interior coordinates
    exterior = np.asarray(poly.exterior.coords)[:, :2]
    interiors = [np.asarray(ring.coords)[:, :2] for ring in poly.interiors]
    # divide all by 1e3 to convert from mm to m
    exterior /= 1e3
    interiors = [interior / 1e3 for interior in interiors]

    # draw exterior
    # convert interior and exterior to a list of lists like [[x1, y1], [x2, y2], ...]
    exterior_coords = [[x, y] for x, y in exterior]
    interiors_coords = [[[x, y] for x, y in interior] for interior in interiors]

    # remove last point if it is the same as the first
    if exterior_coords[0] == exterior_coords[-1]:
        exterior_coords = exterior_coords[:-1]

    for interior in interiors_coords:
        if interior[0] == interior[-1]:
            interior = interior[:-1]
    
    simulation.draw_polygon(exterior_coords, layer)
    for interior in interiors_coords:
        simulation.draw_polygon(interior, layer)

class ActionKiCadPlugin(pcbnew.ActionPlugin):
    def defaults(self):
        self.name = "nikfemm PCB FEM analysis"
        self.category = "PCB FEM analysis"
        self.description = "A plugin for power density and voltage drop analysis of PCBs"
        self.show_toolbar_button = True
        self.plugin_path = os.path.dirname(__file__)
        self.icon_file_name = os.path.join(self.plugin_path, "icon.png")
        self.dark_icon_file_name = os.path.join(self.plugin_path, "icon.png")

        current_dir = os.path.dirname(os.path.abspath(__file__))
        sys.path.append(current_dir)

    def Run(self):
        import shapely
        from nikfemm import MultiLayerCurrentDensitySimulation
        import numpy as np

        try:
            print("###############################################################")

            from Get_PCB_Elements import Get_PCB_Elements
            from Connect_Nets import Connect_Nets
            from Get_Outlines import Get_Outlines

            board = pcbnew.GetBoard()
            connect = board.GetConnectivity()
            ItemList = Get_PCB_Elements(board, connect)
            data = Connect_Nets(ItemList)


            Selected = [d for uuid, d in list(data.items()) if d["IsSelected"]]

            # get the net of the first selected item
            net = Selected[0]["Netname"]

            merged_shapes, drills = Get_Outlines(data, choosen_net=net)

            simulation = MultiLayerCurrentDensitySimulation(2, [35e-6, 35e-6])

            # since some layers are empty, create a map that maps a list of non-empty layers to the layer number
            layer_map = {}
            new_layer = 0
            for layer, shape in enumerate(merged_shapes):
                if shape is not None:
                    if layer not in layer_map:
                        layer_map[layer] = new_layer
                        new_layer += 1
            
            # merged_shapes are MULTIPOLYGONS
            # we get a list of points for each polygon in the multipolygon and a layer number
            for layer, shape in enumerate(merged_shapes):
                if shape is None:
                    continue
                if isinstance(shape, shapely.geometry.Polygon):
                    plot_polygon(shape, layer_map[layer])
                elif isinstance(shape, shapely.geometry.MultiPolygon):
                    for poly in shape.geoms:
                        plot_polygon(poly, layer_map[layer])
                else:
                    print("unknown shape", type(shape))
            
            # drills are like { "Center": d["Position"], "Width": d["Drill"], "Layer": d["Layer"]}
            for drill in drills:
                # divide by 1e3 to convert from mm to m
                simulation.add_interconnection([drill["Center"][0] / 1e3, drill["Center"][1] / 1e3],
                                            [drill["Center"][0] / 1e3, drill["Center"][1] / 1e3],
                                            layer_map[drill['Layer'][0]], layer_map[drill['Layer'][1]], 0.01)

            system, triangles, vertices = simulation.generate_system(False, 1, 10)

            dlg = wx.MessageDialog(
                None,
                f"Net: {net}",
                "Info",
                wx.OK | wx.ICON_INFORMATION,
            )
            dlg.ShowModal()
            dlg.Destroy()

        except Exception as e:
            dlg = wx.MessageDialog(
                None,
                traceback.format_exc(),
                "Fatal Error",
                wx.OK | wx.ICON_ERROR,
            )
            dlg.ShowModal()
            dlg.Destroy()

        pcbnew.Refresh()


if not __name__ == "__main__":
    ActionKiCadPlugin().register()
