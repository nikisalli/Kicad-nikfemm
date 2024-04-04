import pcbnew
from importlib import reload
import sys
import os
import wx
import wx.grid as grid
import wx.lib.agw.floatspin as floatspin
import traceback
import numpy as np
import time
import networkx as nx
import shapely
import shapely.ops
import plotly.graph_objects as go
from nikfemm import MultiLayerCurrentDensitySimulation

ToMM = pcbnew.ToMM

# Überprüfe, ob der Punkt sich innerhalb des Polygons befindet
def IsPointInPolygon(point_, polygon_):
    point = np.array(point_)
    polygon = np.array(polygon_)

    n = len(polygon)
    inside = False

    p1x, p1y = polygon[0]
    for i in range(n + 1):
        p2x, p2y = polygon[i % n]
        if point[1] > min(p1y, p2y):
            if point[1] <= max(p1y, p2y):
                if point[0] <= max(p1x, p2x):
                    if p1y != p2y:
                        x_intersect = (point[1] - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or point[0] <= x_intersect:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

def getHash(obj):
    return obj.m_Uuid.Hash()

def getHashList(objlist):
    return [getHash(obj) for obj in objlist]

def getPolygon(obj):
    poly_obj = obj.GetEffectivePolygon()
    Polygon = [ToMM(poly_obj.CVertex(p)) for p in range(poly_obj.FullPointCount())]
    return Polygon

def getLayer(obj, PossibleLayer=set([0, 31])):
    return sorted(set(obj.GetLayerSet().CuStack()) & PossibleLayer)

def getConnections(track, connect):
    def getVectorLen(vector):
        return np.sqrt(vector.dot(vector))

    def getDistance(point1, point2):
        return getVectorLen(np.array(point2) - np.array(point1))

    def MoveToObjCenter(wirePos, width, objPos):
        objPos = np.array(objPos)
        wirePos = np.array(wirePos)

        diffVector = objPos - wirePos
        # if getVectorLen(diffVector) > width / 2:
        #     return wirePos + width / 2 * diffVector / getVectorLen(diffVector)
        # else:
        #     return wirePos

        x = np.sign(diffVector[0]) * min([abs(diffVector[0]), width / 2])
        y = np.sign(diffVector[1]) * min([abs(diffVector[1]), width / 2])
        return wirePos + np.array([x, y])

    ConnStart = []
    ConnEnd = []

    Start = ToMM(track.GetStart())
    End = ToMM(track.GetEnd())

    for con in connect.GetConnectedTracks(track):
        if type(con) is pcbnew.PCB_VIA:
            print(ToMM(con.GetWidth()))
            print(ToMM(con.GetPosition()))
        elif type(con) is pcbnew.PCB_TRACK:
            conStart = ToMM(con.GetStart())
            conEnd = ToMM(con.GetEnd())
            if Start == conStart:
                ConnStart.append(getHash(con))
            if Start == conEnd:
                ConnStart.append(getHash(con))
            if End == conStart:
                ConnEnd.append(getHash(con))
            if End == conEnd:
                ConnEnd.append(getHash(con))

            if getHash(con) not in ConnStart + ConnEnd:
                distance = [
                    getDistance(Start, conStart),
                    getDistance(Start, conEnd),
                    getDistance(End, conStart),
                    getDistance(End, conEnd),
                ]
                minDis = min(distance)

                if distance[0] == minDis or distance[1] == minDis:
                    ConnStart.append(getHash(con))
                else:
                    ConnEnd.append(getHash(con))

    for con in connect.GetConnectedPads(track):
        Polygon = getPolygon(con)
        Start_ = MoveToObjCenter(Start, ToMM(track.GetWidth()), ToMM(con.GetPosition()))
        End_ = MoveToObjCenter(End, ToMM(track.GetWidth()), ToMM(con.GetPosition()))

        if IsPointInPolygon(Start_, Polygon):
            ConnStart.append(getHash(con))
        if IsPointInPolygon(End_, Polygon):
            ConnEnd.append(getHash(con))

    return ConnStart, ConnEnd

def Get_PCB_Elements(board: pcbnew.BOARD, connect: pcbnew.CONNECTIVITY_DATA):
    # RunSimulation()

    start = time.time()

    DesignSettings = board.GetDesignSettings()
    BoardThickness = ToMM(DesignSettings.GetBoardThickness())
    print("BoardThickness", BoardThickness)
    PossibleLayer = set(DesignSettings.GetEnabledLayers().CuStack())

    print("GetTracks", len(board.GetTracks()))
    print("GetAreaCount", board.GetAreaCount())
    print("GetPads", len(board.GetPads()))
    print("AllConnectedItems", len(board.AllConnectedItems()))
    print("GetFootprints", len(board.GetFootprints()))
    print("GetDrawings", len(board.GetDrawings()))
    print("GetAllNetClasses", len(board.GetAllNetClasses()))

    ItemList = {}

    for track in board.GetTracks():
        temp = {"Layer": getLayer(track, PossibleLayer)}
        if type(track) is pcbnew.PCB_VIA:
            temp["type"] = "VIA"
            temp["Position"] = ToMM(track.GetStart())
            temp["Drill"] = ToMM(track.GetDrill())
            temp["Width"] = ToMM(track.GetWidth())
            temp["connStart"] = sorted(
                getHashList(connect.GetConnectedPads(track))
                + getHashList(connect.GetConnectedTracks(track))
            )
            temp["Area"] = 0
        elif type(track) is pcbnew.PCB_TRACK:
            temp["type"] = "WIRE"
            temp["Start"] = ToMM(track.GetStart())
            temp["End"] = ToMM(track.GetEnd())
            temp["Width"] = ToMM(track.GetWidth())
            temp["Length"] = ToMM(track.GetLength())
            temp["Area"] = temp["Width"] * temp["Length"]
            if track.GetLength() == 0:
                continue
            temp["Layer"] = [track.GetLayer()]
            temp["connStart"], temp["connEnd"] = getConnections(track, connect)
        elif type(track) is pcbnew.PCB_ARC:
            temp["type"] = "WIRE"
            temp["Start"] = ToMM(track.GetStart())
            temp["End"] = ToMM(track.GetEnd())
            temp["Radius"] = ToMM(track.GetRadius())
            temp["Width"] = ToMM(track.GetWidth())
            temp["Length"] = ToMM(track.GetLength())
            temp["Area"] = temp["Width"] * temp["Length"]
            if track.GetLength() == 0:
                continue
            temp["Layer"] = [track.GetLayer()]
            temp["connStart"], temp["connEnd"] = getConnections(track, connect)
        else:
            print("type", type(track), "is not considered!")
            continue

        temp["Netname"] = track.GetNetname()
        temp["NetCode"] = track.GetNetCode()
        temp["id"] = getHash(track)
        temp["IsSelected"] = track.IsSelected()
        ItemList[temp["id"]] = temp
    
    # build a list of segments and their layer and net
    # this is used to remove duplicate segments
    segments = set()
    duplicate_segments = set()
    for Pad in board.AllConnectedItems():
        temp = {"Layer": getLayer(Pad, PossibleLayer)}
        if type(Pad) is pcbnew.ZONE:
            for layer in temp["Layer"]:
                poly = Pad.GetFill(layer)
                if poly.OutlineCount() > 0:
                    for i in range(poly.OutlineCount()):
                        for j in range(poly.Outline(i).SegmentCount()):
                            segment = poly.Outline(i).CSegment(j)
                            elem = ((segment.A[0], segment.A[1]), (segment.B[0], segment.B[1]), layer, Pad.GetNetname())
                            inv_elem = ((segment.B[0], segment.B[1]), (segment.A[0], segment.A[1]), layer, Pad.GetNetname())
                            if elem in segments or inv_elem in segments:
                                duplicate_segments.add(elem)
                            else:
                                segments.add(elem)

    for Pad in board.AllConnectedItems():
        temp = {"Layer": getLayer(Pad, PossibleLayer)}
        if type(Pad) is pcbnew.PAD:
            temp["type"] = "PAD"
            temp["Shape"] = Pad.GetShape()
            # temp["PadAttr"] = Pad.ShowPadAttr()
            # temp["IsFlipped"] = Pad.IsFlipped()
            temp["Position"] = ToMM(Pad.GetPosition())
            temp["Size"] = ToMM(Pad.GetSize())
            temp["Orientation"] = Pad.GetOrientation().AsDegrees()
            temp["DrillSize"] = ToMM(Pad.GetDrillSize())
            temp["Drill"] = temp["DrillSize"][0]
            temp["Area"] = ToMM(ToMM(Pad.GetEffectivePolygon().Area()))
            temp["PadName"] = Pad.GetPadName()
            
            if Pad.GetParentFootprint():
                temp["FootprintReference"] = Pad.GetParentFootprint().GetReference()

        elif type(Pad) is pcbnew.ZONE:
            # pcbnew.ZONE().GetZoneName
            if "teardrop" in Pad.GetZoneName():
                continue
            temp["type"] = "ZONE"
            temp["Position"] = ToMM(Pad.GetPosition())
            temp["Area"] = ToMM(ToMM(Pad.GetFilledArea()))
            temp["NumCorners"] = Pad.GetNumCorners()
            temp["ZoneName"] = Pad.GetZoneName()
            temp["Segments"] = {}
            for layer in temp["Layer"]:
                poly = Pad.GetFill(layer)
                if poly.OutlineCount() > 0:
                    temp["Segments"][layer] = []
                for i in range(poly.OutlineCount()):
                    for j in range(poly.Outline(i).SegmentCount()):
                        segment = poly.Outline(i).CSegment(j)
                        elem = ((segment.A[0], segment.A[1]), (segment.B[0], segment.B[1]), layer, Pad.GetNetname())
                        inv_elem = ((segment.B[0], segment.B[1]), (segment.A[0], segment.A[1]), layer, Pad.GetNetname())
                        if elem in duplicate_segments or inv_elem in duplicate_segments:
                            continue
                        else:
                            temp["Segments"][layer].append(
                                (ToMM(segment.A), ToMM(segment.B))
                            )

        else:
            if not type(track) == pcbnew.PCB_TRACK:
                print("type", type(track), "is not considered!")
            continue

        temp["Netname"] = Pad.GetNetname()
        temp["NetCode"] = Pad.GetNetCode()
        temp["id"] = getHash(Pad)
        temp["IsSelected"] = Pad.IsSelected()
        temp["connStart"] = sorted(
            getHashList(connect.GetConnectedPads(Pad))
            + getHashList(connect.GetConnectedTracks(Pad))
        )
        ItemList[temp["id"]] = temp

    for uuid, d in list(ItemList.items()):  # TODO: WIRES still need to be considered
        if d["type"] == "ZONE":
            for item in d["connStart"]:
                if not "connEND" in ItemList[item]:
                    ItemList[item]["connStart"].append(uuid)

    # now for every ZONE we need to find the closed polygons and the holes in it
    for uuid, d in list(ItemList.items()):
        ItemList[uuid]["Polygons"] = {}
        if d["type"] == "ZONE":
            for layer in d["Segments"].keys():
                # get all segments for this layer
                segments = d["Segments"][layer].copy()
                # get all segments that are connected to other segments
                # we basically need to find all cycles in an undirected graph
                G = nx.Graph()
                for segment in segments:
                    G.add_edge(segment[0], segment[1])
                cycles = list(nx.cycle_basis(G))
                ItemList[uuid]["Polygons"][layer] = cycles

    print("time", time.time() - start)

    return ItemList


OK = 0
NotYetConnected = 0
ErrorConnection = -1

def Connect_Nets(data):
    for uuid, d in list(data.items()):
        layer = {x: NotYetConnected for x in d["Layer"]}  # connection with start

        if "connStart" in d:
            data[uuid]["netStart"] = dict(layer)
        if "connEnd" in d:
            data[uuid]["netEnd"] = dict(layer)

    def getNet(uuid, conn_uuid, layer, pos: (0, 0)):
        if layer not in data[uuid]["Layer"]:
            return ErrorConnection

        temp = NotYetConnected

        if data[uuid]["type"] == "WIRE" and data[conn_uuid]["type"] == "WIRE":
            if pos == data[uuid]["Start"]:
                temp = data[uuid]["netStart"][layer]
            if temp > NotYetConnected:
                return temp
            if pos == data[uuid]["End"]:
                temp = data[uuid]["netEnd"][layer]
        else:
            if "netStart" in data[uuid] and conn_uuid in data[uuid]["connStart"]:
                temp = data[uuid]["netStart"][layer]
            if temp > NotYetConnected:
                return temp
            if "netEnd" in data[uuid] and conn_uuid in data[uuid]["connEnd"]:
                return data[uuid]["netEnd"][layer]
        return temp

    def setNet(uuid, conn_uuid, layer, newNet: NotYetConnected, pos: (0, 0)):
        if layer not in data[uuid]["Layer"]:
            return ErrorConnection

        if data[uuid]["type"] == "WIRE" and data[conn_uuid]["type"] == "WIRE":
            if pos == data[uuid]["Start"]:
                data[uuid]["netStart"][layer] = newNet
            if pos == data[uuid]["End"]:
                data[uuid]["netEnd"][layer] = newNet
        else:
            if "netStart" in data[uuid] and conn_uuid in data[uuid]["connStart"]:
                data[uuid]["netStart"][layer] = newNet

            if "netEnd" in data[uuid] and conn_uuid in data[uuid]["connEnd"]:
                data[uuid]["netEnd"][layer] = newNet
        return OK

    nodeCounter = 0
    for uuid, d in list(data.items()):
        for layer in data[uuid]["Layer"]:
            if "netStart" not in data[uuid]:
                continue
            if d["netStart"][layer] > NotYetConnected:
                continue

            if "Start" in d:
                pos = d["Start"]
            else:
                pos = (0, 0)

            tempNet = NotYetConnected
            for conn in d["connStart"]:
                tmp = getNet(conn, uuid, layer, pos)
                if tmp > NotYetConnected:
                    tempNet = tmp
                    continue
            if tempNet == NotYetConnected:
                nodeCounter += 1
                tempNet = nodeCounter

            for conn in d["connStart"]:
                setNet(conn, uuid, layer, newNet=tempNet, pos=pos)
            data[uuid]["netStart"][layer] = tempNet

    for uuid, d in list(data.items()):
        for layer in data[uuid]["Layer"]:
            if "netEnd" not in data[uuid]:
                continue
            if d["netEnd"][layer] > NotYetConnected:
                continue

            if "End" in d:
                pos = d["End"]
            else:
                pos = (0, 0)

            tempNet = NotYetConnected
            for conn in d["connEnd"]:
                tmp = getNet(conn, uuid, layer, d["End"])
                if tmp > 0:
                    tempNet = tmp
            if tempNet == NotYetConnected:
                nodeCounter += 1
                tempNet = nodeCounter

            for conn in d["connEnd"]:
                setNet(conn, uuid, layer, newNet=tempNet, pos=d["End"])
            data[uuid]["netEnd"][layer] = tempNet

    return data


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

class NikFEMMFrame(wx.Frame):
    def __init__(self, data, *args, **kw):
        super(NikFEMMFrame, self).__init__(*args, **kw)

        self.InitUI(data)

    def InitUI(self, data):
        self.data = data
        # get the nets
        nets = set()
        for uuid, d in list(data.items()):
            if d["Netname"] != "" and not d["Netname"].startswith("unconnected"):
                nets.add(d["Netname"])
        
        nets = list(nets)
        nets.sort()

        # get the pads
        pads = {}
        for net in nets:
            pads[net] = []
            for uuid, d in list(data.items()):
                if d["Netname"] == net and d["type"] == "PAD":
                    if "FootprintReference" in d:
                        if "PadName" in d:
                            name = d["FootprintReference"] + " Pin " + d["PadName"]
                            # name, position, layer
                            pads[net].append((name, d["Position"], d["Layer"]))
        self.pads = pads

        self.SetSize((1200, 1000))
        
        # the gui should have a widget to select the net
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        topgridbox = wx.GridBagSizer(2, 2)

        # add a label to select the net
        st1 = wx.StaticText(panel, label="Select net:")
        topgridbox.Add(st1, (0, 0), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)

        # add a spinner to select the net
        self.spinner = wx.Choice(panel, choices=[])
        topgridbox.Add(self.spinner, (1, 0), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)
        # set callback when the selection changes
        self.spinner.Bind(wx.EVT_CHOICE, self.OnSelectNet)
        # add the nets to the spinner
        self.spinner.Clear()
        self.spinner.AppendItems(nets)
        # select the first net
        self.spinner.SetSelection(0)

        # add a label for the start simulation button
        st2 = wx.StaticText(panel, label="click to start simulation:")
        topgridbox.Add(st2, (0, 1), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)

        # add a button to start the simulation
        btn = wx.Button(panel, label="Start simulation")
        topgridbox.Add(btn, (1, 1), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)
        btn.Bind(wx.EVT_BUTTON, self.OnStartSimulation)

        # add a widget to select the resistance for the vias and pads
        # label
        st3 = wx.StaticText(panel, label="Select Via and Pad resistance (in Ohm):")
        topgridbox.Add(st3, (0, 2), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)
        # spinner
        self.resistance_spinner = floatspin.FloatSpin(panel, value=0.01, increment=0.01, digits=6, min_val=1e-6, max_val=1e6)
        topgridbox.Add(self.resistance_spinner, (1, 2), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)
        self.via_resistance = 0.01
        # add a callback when the resistance changes
        self.resistance_spinner.Bind(floatspin.EVT_FLOATSPIN, self.OnResistanceChange)

        # add a widget to select the maximum triangle area for the mesher
        st4 = wx.StaticText(panel, label="Select max triangle area (in mm²):")
        topgridbox.Add(st4, (0, 3), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)
        self.max_triangle_area_spinner = floatspin.FloatSpin(panel, value=1, increment=1, digits=4, min_val=1e-4, max_val=1e6)
        topgridbox.Add(self.max_triangle_area_spinner, (1, 3), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)
        self.max_triangle_area = 1
        self.max_triangle_area_spinner.Bind(floatspin.EVT_FLOATSPIN, self.OnMaxTriangleAreaChange)

        # add a widget to select the minimum triangle angle for the mesher
        st5 = wx.StaticText(panel, label="Select min triangle angle (in degrees):")
        topgridbox.Add(st5, (0, 4), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)
        self.min_triangle_angle_spinner = floatspin.FloatSpin(panel, value=30, increment=1, digits=0, min_val=1, max_val=33)
        topgridbox.Add(self.min_triangle_angle_spinner, (1, 4), (1, 1), flag=wx.EXPAND | wx.ALL, border=5)
        self.min_triangle_angle = 15
        self.min_triangle_angle_spinner.Bind(floatspin.EVT_FLOATSPIN, self.OnMinTriangleAngleChange)

        vbox.Add(topgridbox, 0, wx.ALL | wx.EXPAND, 5)

        # add a table to show the pads and select the voltage
        self.table = grid.Grid(panel)
        self.table.CreateGrid(0, 2)
        self.table.SetColLabelValue(0, "Pad")
        self.table.SetColLabelValue(1, "Voltage")
        vbox.Add(self.table, 1, wx.ALL | wx.EXPAND, 5)

        # make table scrollable and stretchable
        self.table.EnableScrolling(True, True)
        self.table.EnableDragColSize(True)
        
        panel.SetSizer(vbox)
    
    def OnSelectNet(self, e):
        # get the selected net
        choosen_net = self.spinner.GetStringSelection()
        
        # get the pads for the selected net
        pads = self.pads[choosen_net]

        # update the table
        self.table.ClearGrid()
        self.table.DeleteRows(0, self.table.GetNumberRows())
        self.table.AppendRows(len(pads))
        for i, (name, pos, layer) in enumerate(pads):
            self.table.SetCellValue(i, 0, name)
            self.table.SetCellValue(i, 1, "-")
    
    def OnResistanceChange(self, e):
        self.via_resistance = self.resistance_spinner.GetValue()
    
    def OnMaxTriangleAreaChange(self, e):
        self.max_triangle_area = self.max_triangle_area_spinner.GetValue() / 1e6
        print("max triangle area mm2", self.max_triangle_area * 1e6)
        print("max triangle area m2", self.max_triangle_area)

    def OnMinTriangleAngleChange(self, e):
        self.min_triangle_angle = self.min_triangle_angle_spinner.GetValue()
    
    def OnStartSimulation(self, e):
        # get the selected net
        choosen_net = self.spinner.GetStringSelection()

        # get the pads for the selected net that have a voltage different from "-"
        choosen_pads = []
        for i in range(self.table.GetNumberRows()):
            name = self.table.GetCellValue(i, 0)
            voltage = self.table.GetCellValue(i, 1)
            if voltage != "-":
                choosen_pads.append([name, voltage])
                # find the pad in the pads list and append the data
                for pad in self.pads[choosen_net]:
                    if pad[0] == name:
                        choosen_pads[-1].append(pad[1]) # position
                        choosen_pads[-1].append(pad[2]) # layer
                        break
        
        print("choosen_net")
        print(choosen_net)
        print("pads")
        print(choosen_pads)

        merged_shapes, drills = Get_Outlines(self.data, choosen_net=choosen_net)

        layer_map = {}
        new_layer = 0
        for layer, shape in enumerate(merged_shapes):
            if shape is not None:
                if layer not in layer_map:
                    layer_map[layer] = new_layer
                    new_layer += 1
        print(layer_map)

        # simulate
        simulation = MultiLayerCurrentDensitySimulation(len(layer_map), [35e-6 for x in layer_map])

        for layer, shape in enumerate(merged_shapes):
            if shape is None:
                continue
            if isinstance(shape, shapely.geometry.Polygon):
                plot_polygon(shape, layer_map[layer], simulation)
            elif isinstance(shape, shapely.geometry.MultiPolygon):
                for poly in shape.geoms:
                    plot_polygon(poly, layer_map[layer], simulation)
            else:
                print("unknown shape", type(shape))
        
        for drill in drills:
            # divide by 1e3 to convert from mm to m
            simulation.add_interconnection([drill["Center"][0] / 1e3, drill["Center"][1] / 1e3],
                                        [drill["Center"][0] / 1e3, drill["Center"][1] / 1e3],
                                        layer_map[drill['Layer'][0]], layer_map[drill['Layer'][1]], self.via_resistance)
    
        system, triangles, vertices = simulation.generate_system(False, self.max_triangle_area, self.min_triangle_angle) # max triangle area, min triangle angle

        # set the voltage for the pads
        for pad in choosen_pads:
            voltage = float(pad[1])
            position = [pos / 1e3 for pos in pad[2]]
            layer = layer_map[pad[3][0]]
            simulation.set_voltage(system, position, voltage, layer)
        
        # solve the system
        voltages, triangles, power_densities = simulation.solve(system)

        # normalize power densities
        power_densities = [p[1] for p in power_densities]
        # clip power densities to 95th percentile
        power_densities = np.clip(power_densities, 0, np.quantile(power_densities, 0.95))

        # plot the results
        fig = go.Figure()
        # Unpack the vertices and triangles
        layers, xs, ys, voltages = zip(*voltages)
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
        fig.show()

        # plot the power densities
        fig = go.Figure()
        # Create a 3D surface plot for the triangles
        mesh = go.Mesh3d(
            x=xs,
            y=ys,
            z=layers,
            i=v1s,
            j=v2s,
            k=v3s,
            intensity=power_densities,
            intensitymode='cell',
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
        fig.show()
    
    def OnClose(self, e):
        self.Destroy()

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
        try:
            reload(pcbnew)

            print("###############################################################")

            board = pcbnew.GetBoard()
            connect = board.GetConnectivity()
            ItemList = Get_PCB_Elements(board, connect)
            data = Connect_Nets(ItemList)

            # make a small gui to select the net and the voltage for the pads
            app = wx.App()
            frame = NikFEMMFrame(data, None, title="nikfemm PCB FEM analysis")

            frame.Show()
            app.MainLoop()

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
