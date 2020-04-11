import OCC
# import OCC.Adaptor2d
import OCC as occ
import inspect
import importlib

# i = importlib.import_module("matplotlib.text")
""" 
g = globals()
# def get_imports():
# inspect.getmembers(inspect.getmembers(OCC)[0][1])[0]

for mod_name, mod in inspect.getmembers(occ):

    # print('OCC.'+ mod_name)
    if mod_name in ['VERSION'] or mod_name[0] == '_':
        continue

    md = importlib.import_module('OCC.'+ mod_name)

    for name, sub in inspect.getmembers(mod):
        if str(name)[0] == '_' or 'swig' in name.lower() or \
                name in ['new_instancemethod', 'register_handle']:
            continue

        if callable(sub) is True:

            print('OCC.' + mod_name + '.' + name)
            g[mod_name + '.' + name] = sub

        # importlib.import_module(sub)
"""

from OCC.BRepBuilderAPI import (
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_Transform,
    BRepBuilderAPI_Sewing,
    BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeSolid,
    BRepBuilderAPI_MakeShell,
    BRepBuilderAPI_MakeEdge2d,
    BRepBuilderAPI_NurbsConvert,
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeVertex,
    BRepBuilderAPI_MakeShape,

    BRepBuilderAPI_FindPlane,
    BRepBuilderAPI_Collect,
    BRepBuilderAPI_ModifyShape,
    BRepBuilderAPI_RightCorner,
    BRepBuilderAPI_Trimmed,
    BRepBuilderAPI_VertexInspector,
    BRepBuilderAPI_BndBoxTreeSelector,
    BRepBuilderAPI_GTransform,
    BRepBuilderAPI_NurbsConvert,
    BRepBuilderAPI_Command,

)

from OCC.TopoDS import TopoDS_Face, TopoDS_Vertex, TopoDS_Edge

from OCC.BRepExtrema import (
    BRepExtrema_DistShapeShape,
    BRepExtrema_ExtCC,
    BRepExtrema_DistanceSS,
    BRepExtrema_ExtCF,
    BRepExtrema_ExtFF,
    BRepExtrema_ExtPC,
    BRepExtrema_ExtPF,
    BRepExtrema_Poly,
    BRepExtrema_Poly_Distance,
    BRepExtrema_ExtFF,
    BRepExtrema_IsInFace,
    BRepExtrema_IsOnEdge,
    BRepExtrema_IsVertex,


)

from OCC.BRepFill import (
    BRepFill_NSections,
    BRepFill_Filling,
    brepfill_Face,
    BRepFill_CurveConstraint
)

from OCC.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakePrism,
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeHalfSpace,
    BRepPrimAPI_MakeOneAxis,
    BRepPrimAPI_MakeSweep,

)
from OCC.BRepOffsetAPI import (
    BRepOffsetAPI_MakeEvolved,
    BRepOffsetAPI_MakePipe,
    BRepOffsetAPI_ThruSections,
    BRepOffsetAPI_MakeOffset,
    BRepOffsetAPI_NormalProjection
)

from OCC.BRep import (
    BRep_Tool,
    BRep_Builder,

    BRep_Curve3D,
    BRep_CurveOn2Surfaces,
    BRep_CurveOnClosedSurface,
    BRep_CurveOnSurface,
    BRep_GCurve,
    BRep_CurveRepresentation,

    BRep_PointOnCurve,
    BRep_PointOnCurveOnSurface,
    BRep_PointsOnSurface,
    BRep_PointOnSurface,
    BRep_PointRepresentation,

    BRep_ListNodeOfListOfPointRepresentation,
    BRep_ListIteratorOfListOfCurveRepresentation,
    BRep_ListIteratorOfListOfPointRepresentation,
    BRep_ListNodeOfListOfCurveRepresentation,
    BRep_ListOfCurveRepresentation,
    BRep_ListOfPointRepresentation,

    BRep_TFace,
    BRep_TVertex,
    BRep_TEdge,

    BRep_Polygon3D,
    BRep_PolygonOnSurface,
    BRep_PolygonOnClosedSurface,
    BRep_PolygonOnTriangulation,
    BRep_PolygonOnClosedTriangulation

)


from OCC.BRepFeat import brepfeat
from OCC.BRepTools import (
    breptools,
    BRepTools_Modifier,
    BRepTools_Quilt,
    BRepTools_DataMapNodeOfMapOfVertexPnt2d,
    BRepTools_DataMapIteratorOfMapOfVertexPnt2d,
    BRepTools_MapOfVertexPnt2d,
    BRepTools_ReShape,
    BRepTools_ShapeSet,
    BRepTools_Substitution,
    BRepTools_WireExplorer,
    BRepTools_Modification,
    BRepTools_TrsfModification,
    BRepTools_GTrsfModification,
    BRepTools_NurbsConvertModification
)

# Boolean
from typing import Union


Brep_Prim = Union[
BRep_TFace,
    BRep_TVertex,
    BRep_TEdge,
]


Brep_OnSurface = Union[
    BRep_CurveOnSurface,
    BRep_CurveOnClosedSurface,
    BRep_PointOnCurveOnSurface,
    BRep_PointsOnSurface,
    BRep_PointOnSurface,
]


Brep_Point = [
    BRep_PointOnCurve,
    BRep_PointOnSurface
]

Brep_Curve = Union[
    BRep_Curve3D,
    BRep_CurveOn2Surfaces,
    BRep_GCurve,
    BRep_CurveRepresentation,
]

Brep_Representation = Union[
    BRep_PointRepresentation,
]


Brep_OnTriangulation = Union[
    BRep_PolygonOnTriangulation,
    BRep_PolygonOnClosedTriangulation
]


Brep_Polygon = Union[
    BRep_PolygonOnClosedSurface,
    BRep_Polygon3D,
    BRep_PolygonOnSurface,
    BRep_PolygonOnTriangulation,
    BRep_PolygonOnClosedTriangulation
]

Brep_List = Union[
    BRep_ListNodeOfListOfPointRepresentation,
    BRep_ListIteratorOfListOfCurveRepresentation,
    BRep_ListIteratorOfListOfPointRepresentation,
    BRep_ListNodeOfListOfCurveRepresentation,
    BRep_ListOfCurveRepresentation,
    BRep_ListOfPointRepresentation
]

