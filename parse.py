import numpy as np
import pcbnew
# from pprint import pprint
from pathlib import Path
# from importlib import reload
import shapely
from Get_PCB_Elements import Get_PCB_Elements
from Connect_Nets import Connect_Nets
# from Get_PCB_Stackup import Get_PCB_Stackup
from Get_Outlines import Get_Outlines

from nikfemm import MultiLayerCurrentDensitySimulation

import time

time1 = time.time()

path = "nik-servo-driver-v1.kicad_pcb"
debug = False

board = pcbnew.LoadBoard(path)

time2 = time.time()

connect = board.GetConnectivity()

time3 = time.time()

ItemList = Get_PCB_Elements(board, connect)

time4 = time.time()

data = Connect_Nets(ItemList)

time5 = time.time()

# PhysicalLayerStack, CuStack = Get_PCB_Stackup(ProjectPath=board_FileName)

merged_shapes, drills = Get_Outlines(data, choosen_net="GND")

time6 = time.time()

simulation = MultiLayerCurrentDensitySimulation(2, [35e-6, 35e-6])

time7 = time.time()

def plot_polygon(poly, layer):
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

# since some layers are empty, create a map that maps a list of non-empty layers to the layer number
layer_map = {}
new_layer = 0
for layer, shape in enumerate(merged_shapes):
    if shape is not None:
        if layer not in layer_map:
            layer_map[layer] = new_layer
            new_layer += 1
print(layer_map)

time8 = time.time()

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

time9 = time.time()

# drills are like { "Center": d["Position"], "Width": d["Drill"], "Layer": d["Layer"]}
for drill in drills:
    # divide by 1e3 to convert from mm to m
    simulation.add_interconnection([drill["Center"][0] / 1e3, drill["Center"][1] / 1e3],
                                   [drill["Center"][0] / 1e3, drill["Center"][1] / 1e3],
                                   layer_map[drill['Layer'][0]], layer_map[drill['Layer'][1]], 0.01)

time10 = time.time()

system, triangles, vertices = simulation.generate_system(False, 1, 10)

time11 = time.time()

# plot meshes
# import matplotlib.pyplot as plt
# for num, trilist in enumerate(triangles):
#     # trilist is a list of tuples of indices [(vertex1, vertex2, vertex3),...]
#     # vertices is a list of coordinates [[x1, y1], [x2, y2],...]
# 
#     # plot the current mesh using matplotlib
# 
#     fig, ax = plt.subplots()
#     ax.set_aspect('equal')
#     for triangle in trilist:
#         # get the coordinates of the triangle
#         x = [vertices[triangle[0]][0], vertices[triangle[1]][0], vertices[triangle[2]][0], vertices[triangle[0]][0]]
#         y = [vertices[triangle[0]][1], vertices[triangle[1]][1], vertices[triangle[2]][1], vertices[triangle[0]][1]]
#         ax.plot(x, y, 'k-')
# 
#     plt.show()

# set voltages
simulation.set_voltage(system, [149.77e-3, 73e-3], 1, 0)
simulation.set_voltage(system, [147e-3, 55e-3], -1, 1)

time12 = time.time()

# simulation.set_voltage(system, [150e-3, 32e-3], 1, 0)
# simulation.set_voltage(system, [160e-3, 30e-3], -1, 1)

# solve
verts, triangles = simulation.solve(system)

time13 = time.time()

# print(verts)
# print(triangles)

# 3d plot using plotly
# verts are the vertices of the mesh [(layer, x, y, voltage), ...]
# triangles are the indices of the vertices that make up the triangles of the mesh [(layer, vertex1, vertex2, vertex3), ...]

import plotly.graph_objects as go
import plotly.express as px
import plotly

fig = go.Figure()

# Unpack the vertices and triangles
layers, xs, ys, voltages = zip(*verts)
_, v1s, v2s, v3s = zip(*triangles)

# divide layers by 1e3 to convert from m to mm
layers = [layer / 1e3 for layer in layers]

# Create a 3D surface plot for the triangles
mesh = go.Mesh3d(
    x=xs,
    y=ys,
    z=layers,
    i=v1s,
    j=v2s,
    k=v3s,
    intensity=voltages,
    intensitymode='vertex',
    colorscale='Jet',
    opacity=0.9,
)

fig.update_layout(
    scene=dict(
        aspectmode='data',
        aspectratio=dict(x=1, y=1, z=1)
    )
)

# Add the mesh to the figure
fig.add_trace(mesh)

# Show the plot
plotly.offline.plot(fig, filename='nik-servo-driver-v1.html', auto_open=False)

time14 = time.time()

# print the data above but in milliseconds and in a more readable alignment and rounded to 3 decimal places
print("LoadBoard:                          ", round((time2 - time1) * 1e3, 3), "ms")
print("GetConnectivity:                    ", round((time3 - time2) * 1e3, 3), "ms")
print("GetPCBElements:                     ", round((time4 - time3) * 1e3, 3), "ms")
print("ConnectNets:                        ", round((time5 - time4) * 1e3, 3), "ms")
print("GetOutlines:                        ", round((time6 - time5) * 1e3, 3), "ms")
print("MultiLayerCurrentDensitySimulation: ", round((time7 - time6) * 1e3, 3), "ms")
print("layer_map:                          ", round((time8 - time7) * 1e3, 3), "ms")
print("draw_polygon:                       ", round((time9 - time8) * 1e3, 3), "ms")
print("add_interconnection:                ", round((time10 - time9) * 1e3, 3), "ms")
print("generate_system:                    ", round((time11 - time10) * 1e3, 3), "ms")
print("set_voltage:                        ", round((time12 - time11) * 1e3, 3), "ms")
print("solve:                              ", round((time13 - time12) * 1e3, 3), "ms")
print("@ Simulation total:                 ", round((time13 - time1) * 1e3, 3), "ms")
print("plotly:                             ", round((time14 - time13) * 1e3, 3), "ms")
print("@ plot and simulation Total:        ", round((time14 - time1) * 1e3, 3), "ms")
