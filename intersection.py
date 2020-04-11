from OCC.BRepAlgoAPI import (BRepAlgoAPI_Common,
                             BRepAlgoAPI_Section,
                             BRepAlgoAPI_Fuse,
                             BRepAlgoAPI_Cut)
from OCC.TopTools import TopTools_ListOfShape, TopTools_ListIteratorOfListOfShape, TopTools_ListNodeOfListOfShape
from .oci import topods as tds
from OCC.TopExp import topexp
from OCC.TopoDS import topods
from .oci.topods import TopoDS_Shape
from .oci.topabs import *
from OCC.Geom import Handle_Geom_Surface, Handle_Geom_Curve
from OCC.GeomAPI import GeomAPI_IntSS, GeomAPI_IntCS
from OCC.IntTools import (IntTools_FaceFace,
                          IntTools_EdgeFace,
                          IntTools_EdgeEdge,
                          IntTools_CommonPrt)
# from OCC import Topo
from OCC.BRepFeat import BRepFeat_SplitShape
from .Common import assert_isdone
from .oci import brep
from collections import defaultdict
from .oci import gp
from OCC.IntCurvesFace import IntCurvesFace_ShapeIntersector
# from OCC.IntWalk import
from . import Construct, Topo, ShapeDict
from typing import Union, Optional

"""
# GeomAPI_IntCS
# src/IntCurvesFace/IntCurvesFace_Intersector.cxx

# BRepIntCurveSurface

BRepAlgoAPI_Common
 	GeomInt
 	
https://www.opencascade.com/content/use-brepalgosection-instead-brepalgoapisection
BRepAlgo_Section

"""


# --------------------------------------------------------------------
# topoDS
def _bool_op(op) -> TopoDS_Shape:
    if op.BuilderCanWork() is True:
        op.Build()
        op.RefineEdges()
        op.FuseEdges()
        shp = op.Shape()
        op.Destroy()
        return shp


def union(shp1: TopoDS_Shape, shp2: TopoDS_Shape, merge=False) -> TopoDS_Shape:
    """ shp1 + shape2 """
    intrs = BRepAlgoAPI_Fuse(shp1, shp2)
    return _bool_op(intrs)


def difference(shp1: TopoDS_Shape, shp2: TopoDS_Shape) -> TopoDS_Shape:
    """ shp1 cuts shape2 """
    intrs = BRepAlgoAPI_Cut(shp1, shp2)
    return _bool_op(intrs)


def intersection_solid(shp1: tds.TopoDS_Shape, shp2: TopoDS_Shape) -> TopoDS_Shape:
    """ Solid Solid intersection
    :param: shp1 - TopoDS_Shape1
    :param: shp2 - TopoDS_Shape1
    """
    return _bool_op(BRepAlgoAPI_Common(shp1, shp2))


def intersection(shp1: TopoDS_Shape,
                 shp2: Union[TopoDS_Shape, gp.GP3d]) -> TopoDS_Shape:
    """
    Most Robust TopoDS intersection

    BRepAlgoAPI_Common will only return if the intersection is solid.
    BRepAlgoAPI_Section will work for face-on-face

    similar issue with GeomAPI is documented here:
        https://www.opencascade.com/content/use-brepalgosection-instead-brepalgoapisection

    :param: shp1 - TopoDS_Shape1
    :param: shp2 - TopoDS_Shape2

    BRepAlgoAPI_Section(TopoDS_Shape const &,TopoDS_Shape const &,BOPAlgo_PaveFiller const &,Standard_Boolean const)
    BRepAlgoAPI_Section(TopoDS_Shape const &,TopoDS_Shape const &,Standard_Boolean const)
    BRepAlgoAPI_Section(TopoDS_Shape const &,gp_Pln const &,Standard_Boolean const)
    BRepAlgoAPI_Section(TopoDS_Shape const &,Handle_Geom_Surface const &,Standard_Boolean const)
    BRepAlgoAPI_Section(Handle_Geom_Surface const &,TopoDS_Shape const &,Standard_Boolean const)
    BRepAlgoAPI_Section(Handle_Geom_Surface const &,Handle_Geom_Surface const &,Standard_Boolean const)

    returns wires representing the intersection
    """
    intrs = BRepAlgoAPI_Section(shp1, shp2)
    if intrs.BuilderCanWork() is True:
        intrs.Build()
        # todo add isDone check (maybe ??)
        # intrs.FuseEdges()
        intrs.RefineEdges()
        shp = intrs.Shape()
        intrs.Destroy()
        return shp



def split2(base, cutters):
    # https://www.opencascade.com/doc/occt-7.0.0/overview/html/
    # occt_user_guides__boolean_operations.html#occt_algorithms_10a
    builder = BRepAlgoAPI_Section()
    tools = TopTools_ListOfShape()
    for cutter in cutters:
        tools.Append(cutter)
    builder.ComputePCurveOn1()
    # builder.SetArguments(base)
    # builder.SetTools(tools)



def split(base, plane):
    splt = BRepFeat_SplitShape()
    splt.Init(base)

    sect = BRepAlgoAPI_Section(base, plane, False)
    sect.ComputePCurveOn1(True)
    sect.Approximation(True)
    sect.Build()
    edge = sect.Shape()

    rdict = set()
    # print(Topo(edge).number_of_edges())

    Ex = TopExp_Explorer(edge, TopAbs_EDGE)

    while Ex.More():
        # print('edge', Ex.Current())
        base_iter = TopExp_Explorer(base, TopAbs_FACE)
        curr = Ex.Current()

        while base_iter.More():

            # print('face', base_iter.Current())
            bface = base_iter.Current()

            if sect.HasAncestorFaceOn1(curr, bface):
                # print('has1', curr, bface)
                k, v = hash(bface), hash(curr)
                if (k, v) not in rdict:

                    rdict.add((k, v))
                    e = topods.Edge(curr)
                    f = topods.Face(bface)
                    splt.Add(e, f)

            if sect.HasAncestorFaceOn2(curr, bface):
                # print('has2', curr, bface)
                pass
            base_iter.Next()
        Ex.Next()
    splt.Build()
    return splt.Shape()


def split_solid(base, plane):
    splt = BRepFeat_SplitShape()
    splt.Init(base)

    sect = BRepAlgoAPI_Section(base, plane, False)
    sect.ComputePCurveOn1(True)
    sect.Approximation(True)
    sect.Build()
    edge = sect.Shape()

    rdict = set()
    # print(Topo(edge).number_of_edges())

    Ex = TopExp_Explorer(edge, TopAbs_EDGE)

    while Ex.More():
        # print('edge', Ex.Current())
        base_iter = TopExp_Explorer(base, TopAbs_FACE)
        curr = Ex.Current()

        while base_iter.More():

            # print('face', base_iter.Current())
            bface = base_iter.Current()

            if sect.HasAncestorFaceOn1(curr, bface):
                # print('has1', curr, bface)
                k, v = hash(bface), hash(curr)
                if (k, v) not in rdict:

                    rdict.add((k, v))
                    e = topods.Edge(curr)
                    f = topods.Face(bface)
                    splt.Add(e, f)

            if sect.HasAncestorFaceOn2(curr, bface):
                # print('has2', curr, bface)
                pass
            base_iter.Next()
        Ex.Next()
    splt.Build()
    return splt.Shape()


class ModificationStep:
    def __init__(self):
        self.replaced = []
        self.added = []

    def edge_dict_updates(self, edge_dict, new_shape, base_shape):
        pass


class Intersector:
    """
    given shape1 and shape2, form the intersection
    need to be able to return which 'topo' of
    """
    def __init__(self, shp1, shp2):
        # self.shp1 = shp1
        # self.shp2 = shp2
        self.algo = BRepAlgoAPI_Section(shp1, shp2, True)
        self.algo.ComputePCurveOn1(True)
        self.algo.ComputePCurveOn2(True)
        # self.algo.Approximation(True)

    @property
    def shp1(self):
        return self.algo.Shape1()

    @property
    def shp2(self):
        return self.algo.Shape2()

    def Shape(self):
        return self.algo.Shape()

    def faces_on(self):
        self.algo.Build()
        from OCC.TopTools import TopTools_ListIteratorOfListOfShape
        seen = set()

        edge_list = self.algo.Modified(self.shp1)
        itre = TopTools_ListIteratorOfListOfShape(edge_list)
        while itre.More():
            print(itre.Value())
            itre.Next()
        print('-')
        edge_list = self.algo.Modified(self.shp2)
        itre = TopTools_ListIteratorOfListOfShape(edge_list)
        while itre.More():
            print(itre.Value())
            itre.Next()
        # edge_list = self.algo.SectionEdges()
        # itr = TopTools_ListIteratorOfListOfShape(edge_list)
        # edge_list = self.algo.SectionEdges()


        res = self.Shape()
        itr = TopExp_Explorer(res, TopAbs_EDGE)



        faces1, faces2 = [], []
        while itr.More():
            curr_edge = itr.Current()

            s1_iter = TopExp_Explorer(self.shp1, TopAbs_FACE)
            s2_iter = TopExp_Explorer(self.shp2, TopAbs_FACE)

            while s1_iter.More():
                curr_face1 = s1_iter.Current()

                if self.algo.HasAncestorFaceOn1(curr_edge, curr_face1):
                    k, v = hash(curr_face1), hash(curr_edge)
                    if (k, v) not in seen:
                        seen.add((k, v))
                        faces1.append(curr_face1)
                s1_iter.Next()

            while s2_iter.More():
                curr_face2 = s2_iter.Current()
                if self.algo.HasAncestorFaceOn2(curr_edge, curr_face2):
                    k, v = hash(curr_face2), hash(curr_edge)
                    if (k, v) not in seen:
                        seen.add((k, v))
                        faces2.append(curr_face2)

                s2_iter.Next()

            # s2_iter.ReInit()
            # s1_iter.ReInit()
            itr.Next()
        return faces1, faces2

    def commonface(self):
        w = Construct.make_wirex(*Topo(self.Shape()).edges())
        f = Construct.make_face(w)
        return f


class Splitter(ModificationStep):
    def __init__(self, base_shape):
        ModificationStep.__init__(self)
        self._base = base_shape
        self._data = {}
        self.cutting_edge = None
        self.ancestors = set()

    def edge_dict_updates(self, edge_dict, new_shape, base_shape):
        # includes 3 new edges,
        # 2 edges - 'replacement' for cut edge
        # 1 edge -  'cutter' the edge
        # and 1 edge that is 'removed'
        new_hash = ShapeDict.of_edges(new_shape)

        for new in self.added:
            # remove the 'cutter' edge from new_hash
            hnew = hash(new)
            new_hash.remove(hnew)
            # add 'cutter' to state dict
            edge_dict[hnew] = -1

        # find the 'removed' edge
        prev_hash = ShapeDict.of_edges(base_shape)
        removed_edges = prev_hash.difference(new_hash)

        assert len(removed_edges) == 1
        removed_id = list(removed_edges.keys())[0]
        for eid in new_hash:
            if eid not in edge_dict:
                # updated the 'replacement' edges with
                edge_dict[eid] = removed_edges[removed_id]

        # remove the 'removed' edge from state dict
        edge_dict.remove(removed_id)

    def __call__(self, plane: gp.gp_Pln, vertex=None, edge_dict=None, **kwargs):
        splt = BRepFeat_SplitShape()
        if isinstance(self._base, TopoDS_Shape):
            base_shape = self._base
        else:
            base_shape = self._base.Shape()

        splt.Init(base_shape)

        sect = BRepAlgoAPI_Section(base_shape, plane, False)
        sect.ComputePCurveOn1(True)
        sect.Approximation(True)
        sect.Build()

        self.cutting_edge = sect.Shape()

        ancestors = set()
        new_faces = []
        edge_iter = TopExp_Explorer(self.cutting_edge, TopAbs_EDGE)

        while edge_iter.More():
            base_iter = TopExp_Explorer(base_shape, TopAbs_FACE)
            curr_edge = edge_iter.Current()

            while base_iter.More():

                curr_face = base_iter.Current()

                if sect.HasAncestorFaceOn1(curr_edge, curr_face):
                    k, v = hash(curr_face), hash(curr_edge)
                    if (k, v) not in self.ancestors:
                        ancestors.add((k, v))
                        e = topods.Edge(curr_edge)
                        f = topods.Face(curr_face)
                        splt.Add(e, f)
                        # todo - only add the closest one !!!!
                        new_faces.append(f)
                        self.added.append(e)
                        break
                # if sect.HasAncestorFaceOn2(curr_edge, curr_face):
                    # print('has2', curr_edge, curr_face)
                #     pass
                base_iter.Next()
            edge_iter.Next()

        # -------------------------------------
        splt.Build()
        new_shape = splt.Shape()
        sect.Destroy()

        return new_shape


class SplitQuilt(Splitter):
    def __init__(self, base_shape):
        Splitter.__init__(self, base_shape)

    def edge_dict_updates(self, edge_dict, new_shape, base_shape):
        # includes 3 new edges,
        # 2 edges - 'replacement' for cut edge
        # 1 edge -  'cutter' the edge
        # and 1 edge that is 'removed'
        new_hash = ShapeDict.of_edges(new_shape)

        for new in self.added:
            # remove the 'cutter' edge from new_hash
            hnew = hash(new)
            new_hash.remove(hnew)
            # add 'cutter' to state dict
            edge_dict[hnew] = -1

        # find the 'removed' edge
        prev_hash = ShapeDict.of_edges(base_shape)
        removed_edges = prev_hash.difference(new_hash)

        assert len(removed_edges) == 1
        removed_id = list(removed_edges.keys())[0]
        for eid in new_hash:
            if eid not in edge_dict:
                # updated the 'replacement' edges with
                edge_dict[eid] = removed_edges[removed_id]

        # remove the 'removed' edge from state dict
        edge_dict.remove(removed_id)

    def __call__(self, plane: gp.gp_Pln, vertex=None, edge_dict=None, **kwargs):
        if isinstance(self._base, TopoDS_Shape):
            base_shape = self._base
        else:
            base_shape = self._base.Shape()

        new_shape = Splitter.__call__(self, plane)
        # -------------------------------------
        # dict updates
        if edge_dict:
            # includes 3 new edges,
            # 2 edges - 'replacement' for cut edge
            # 1 edge -  'cutter' the edge
            # and 1 edge that is 'removed'
            new_hash = ShapeDict.of_edges(new_shape)

            for new in self.added:
                # remove the 'cutter' edge from new_hash
                hnew = hash(new)
                new_hash.remove(hnew)
                # add 'cutter' to state dict
                edge_dict[hnew] = -1

            # find the 'removed' edge
            prev_hash = ShapeDict.of_edges(base_shape)
            removed_edges = prev_hash.difference(new_hash)

            assert len(removed_edges) == 1
            removed_id = list(removed_edges.keys())[0]
            for eid in new_hash:
                if eid not in edge_dict:
                    # updated the 'replacement' edges with
                    edge_dict[eid] = removed_edges[removed_id]

            # remove the 'removed' edge from state dict
            edge_dict.remove(removed_id)

        # -------------------------------------
        # find the edge that is out of order on one of the faces
        # technically, this should be self.new_edge
        FACES = list(Topo(new_shape).faces())
        to_swap = []
        for i, face in enumerate(FACES):
            wire = brep.breptools.OuterWire(face)
            edge_iter = brep.BRepTools_WireExplorer(wire, face)
            orient1 = defaultdict(int)
            cnt = 0
            while edge_iter.More():
                edge_in_new = edge_iter.Current()
                cnt += 1
                orient1[edge_in_new.Orientation()] += 1
                edge_iter.Next()

            bestk, bestv = None, float('inf')
            for k, v in orient1.items():
                if v < bestv:
                    bestk, bestv = k, v

            # only one edge should be disoriented....
            if bestv != cnt:
                edge_iter.Clear()
                edge_iter = brep.BRepTools_WireExplorer(wire, face)
                while edge_iter.More():
                    edge_in_new = edge_iter.Current()
                    if edge_in_new.Orientation() == bestk:
                        to_swap.append([i, edge_in_new])
                    edge_iter.Next()

        print(len(to_swap))
        # once again this should only ever be length one ...
        for i, edge in to_swap:
            va = topexp.FirstVertex(edge)
            vb = topexp.LastVertex(edge)
            new_edge = Construct.make_edge(va, vb)
            if edge_dict:
                edge_dict[hash(new_edge)] = -1

            # for edge in
            # not quite working ....
            reshaper = brep.BRepTools_ReShape()
            reshaper.Replace(edge, new_edge)
            fi = FACES[i]
            FACES[i] = reshaper.Apply(fi)
            print('swap', fi, FACES[i])

        return FACES


def splitwire(base, in_edge):
    sect = BRepAlgoAPI_Section(base, in_edge)
    sect.Build()
    sect.RefineEdges()
    edge = sect.Shape()

    splt = BRepFeat_SplitShape(base)
    Ex = TopExp_Explorer(edge, TopAbs_VERTEX)
    while Ex.More():
        # print(Ex.Current())
        Sx = TopExp_Explorer(base, TopAbs_EDGE)
        while Sx.More():
            if sect.HasAncestorFaceOn1(Ex.Current(), Sx.Current()):
                print('add', Ex.Current(), Sx.Current())
                splt.Add(Ex.Current(), Sx.Current())
            Sx.Next()
        Ex.Next()
    splt.Build()
    return splt.Shape()


def intersect_shape_with_ptdir(occtopology, pypt, pydir):
    """
    This function projects a point in a direction and calculates the at
    which point does the point intersects the OCCtopology.

    Parameters
    ----------
    occtopology : OCCtopology
        The OCCtopology to be projected on.
        OCCtopology includes: OCCshape, OCCcompound, OCCcompsolid,
         OCCsolid, OCCshell, OCCface, OCCwire, OCCedge, OCCvertex

    pypt : tuple of floats
        The point to be projected. A pypt is a tuple that documents the xyz coordinates of a pt e.g. (x,y,z)

    pydir : tuple of floats
        The direction of the point to be projected. A pydir is a tuple that documents the xyz vector of a dir e.g. (x,y,z)

    Returns
    -------
    intersection point : pypt
        The point in which the projected point intersect the OCCtopology. If None means there is no intersection.

    intersection face : OCCface
        The OCCface in which the projected point hits. If None means there is no intersection.
    """
    occ_line = gp_Lin(gp_Ax1(gp_Pnt(pypt[0], pypt[1], pypt[2]), gp_Dir(pydir[0], pydir[1], pydir[2])))
    shape_inter = IntCurvesFace_ShapeIntersector()
    shape_inter.Load(occtopology, 1e-6)
    shape_inter.PerformNearest(occ_line, 0.0, float("+inf"))
    if shape_inter.IsDone():
        npts = shape_inter.NbPnt()
        if npts != 0:
            return None, None# shape_inter.Pnt(1)), shape_inter.Face(1)
        else:
            return None, None
    else:
        return None, None


def touches(shp1: TopoDS_Shape, shp2: TopoDS_Shape) -> bool:
    topo = Topo(intersection(shp1, shp2))
    if topo.number_of_edges() > 1:
        return True
    return False


def _Xtouches(shp1: TopoDS_Shape, shp2: TopoDS_Shape, tol=0.0001):
    """
    So far the fastest one ive found for doing this ...
    :param shp1:
    :param shp2:
    :param tol:
    :return:
    """
    bdss = brep.BRepExtrema_DistShapeShape(shp1, shp2)
    bdss.Perform()
    with assert_isdone(bdss, 'failed computing minimum distances'):
        min_dist = bdss.Value()
        if min_dist <= tol:
            return True

    return False


# --------------------------------------------------------------------
# GEOM
def intersect_shape_by_line(
        shape: TopoDS_Shape,
        line: gp.gp_Lin,
        low_parameter: Optional[float]=0.0,
        hi_parameter: Optional[float]=float("+inf"),
        tol:Optional[float]=0.001
        ):
    """
    finds the intersection of a shape and a line

    :param shape: any TopoDS_*
    :param line: gp_Lin
    :param low_parameter:
    :param hi_parameter:

    IntCurvesFace_ShapeIntersector

    :return: a list with a number of tuples that corresponds to the number
    of intersections found
    the tuple contains ( gp_Pnt, TopoDS_Face, u,v,w ), respectively the
    intersection point, the intersecting face
    and the u,v,w parameters of the intersection point
    :raise:
    """
    shape_inter = IntCurvesFace_ShapeIntersector()
    shape_inter.Load(shape, tol)
    shape_inter.PerformNearest(line, low_parameter, hi_parameter)

    with assert_isdone(shape_inter, "failed to computer shape / line intersection"):
        try:
            return (shape_inter.Pnt(1),
                    shape_inter.Face(1),
                    shape_inter.UParameter(1),
                    shape_inter.VParameter(1),
                    shape_inter.WParameter(1))
        except:
            return None, None, None, None, None


def plane_plane(plane1, plane2) -> gp.gp_Lin:
    # todo plane and plane -> Line
    pass


def line_wire() -> gp.gp_Pnt:
    # from OCC.GeomAPI import Ge
    # GeomAPI_InterCurveCurve
    pass


def line_wire2() -> gp.gp_Pnt:
    # Geom2dAPI_InterCurveCurve
    pass


class BoolOp:
    argmap = {
        (Handle_Geom_Surface, Handle_Geom_Surface): GeomAPI_IntSS,
        (Handle_Geom_Curve, Handle_Geom_Surface): GeomAPI_IntCS,
        (TopoDS_Shape, gp.gp_Pln): IntCurvesFace_ShapeIntersector,

    }
    """
    todod: 
        GeomAbs_Intersection
        IntAna, 
        IntTools
        Intf
        IntSurf
        IntCurvesFace
        BRepIntCurveSurface_Inter
    """
    def __init__(self, obj1, obj2):
        self.obj1 = obj1
        self.obj2 = obj2
        t1 = type(obj1)
        t2 = type(obj2)
        self.cls = None
        if (t1, t2) in self.argmap:
            self.cls = self.argmap[(t1, t2)]

        elif (t2, t1) in self.argmap:
            self.cls = self.argmap[(t2, t1)]

        if self.cls is None:
            raise Exception(
                'could not find a valid BOP Algo for types {} {}'.format(t1, t2)
            )
        self.algo = self.cls()
        # if isinstance(obj1, ()):
        #     Handle_Geom_Surface

    def intersection(self):
        return

    def union(self):
        return

    def difference(self):
        return



