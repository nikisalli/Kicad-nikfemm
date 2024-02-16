import numpy as np
import pcbnew
import time
import networkx as nx

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


# def getMember(obj, small=True):
#     import inspect
#     from pprint import pprint

#     members = inspect.getmembers(obj, predicate=inspect.ismethod)
#     if small:
#         for m in members:
#             print(m[0])
#     else:
#         pprint(members)


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
            # temp["FootprintUUID"] = getHash(Pad.GetParent())
            # if Pad.GetParent():
            #     temp["FootprintReference"] = Pad.GetParent().GetReference()

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
