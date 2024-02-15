import pcbnew
from importlib import reload
import sys
import os
import wx
import traceback
from pprint import pprint
from pathlib import Path
import math


# import pip
# def install(package):
#     if hasattr(pip, "main"):
#         pip.main(["install", package])
#     else:
#         pip._internal.main(["install", package])
# install("PySpice")
# import PySpice

debug = 0


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
            print("###############################################################")

            from Get_PCB_Elements import Get_PCB_Elements
            from Connect_Nets import Connect_Nets
            from Connect_Nets import Connect_Nets

            board = pcbnew.GetBoard()
            connect = board.GetConnectivity()
            ItemList = Get_PCB_Elements(board, connect)
            data = Connect_Nets(ItemList)


            Selected = [d for uuid, d in list(data.items()) if d["IsSelected"]]

            dlg = wx.MessageDialog(
                None,
                "Selected: {}".format(len(Selected)),
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
