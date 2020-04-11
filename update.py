from OCC.ShapeFix import ShapeFix_Wire
from OCC.ShapeAnalysis import ShapeAnalysis_Wire
from OCC.ShapeExtend import ShapeExtend_WireData
from OCC.BRepTools import BRepTools_WireExplorer

from lib.OCCUtils import Topo, TopoMap, Wire, Edge, Construct, Vertex, Common
from lib.OCCUtils import intersection as bop
from .oci.topabs import *
from .oci import shape
from lib.OCCUtils.oci.topods import *

from lib.OCCUtils.oci import gp, brep
from collections import defaultdict
from typing import Tuple, Union


def copy_reversed(wire, edge_dict=None):
    if edge_dict is not None:
        return copy_reversed_params(wire, edge_dict)
    itr = BRepTools_WireExplorer(wire)
    verts = []
    while itr.More():
        current_item = itr.Current()
        out_fst = topexp.FirstVertex(current_item)
        # out_lst = topexp.LastVertex(current_item)
        verts.append(out_fst)
        itr.Next()

    verts.reverse()
    return Wire.create(*verts)


def copy_reversed_params(wire, edge_dict):
    itr = BRepTools_WireExplorer(wire)
    verts = []
    while itr.More():
        current_item = itr.Current()
        out_fst = topexp.FirstVertex(current_item)
        out_lst = topexp.LastVertex(current_item)

        new_edge = Construct.make_edge(
            Construct.copy_vertex(out_lst),
            Construct.copy_vertex(out_fst)
        )

        verts.append(new_edge)
        edge_dict[hash(new_edge)] = edge_dict[hash(current_item)]
        itr.Next()

    # verts.reverse()
    return Construct.make_wirex(*verts)


def follow(wire, face=None):
    if face:
        itr = BRepTools_WireExplorer(wire,  face)
    else:
        itr = BRepTools_WireExplorer(wire)
    verts = []
    print('------------')
    while itr.More():
        current_item = itr.Current()

        out_fst = topexp.FirstVertex(current_item)
        out_lst = topexp.LastVertex(current_item)
        v1 = Common.vertex2pnt(out_fst)
        v2 = Common.vertex2pnt(out_lst)

        print(gp.repr3(v1), gp.repr3(v2), current_item.Orientation())

        # verts.append(out_fst)
        itr.Next()
    print('------------')
    return


class RemoveHole(bop.ModificationStep):
    def compute_from_wires(self, inner_wire, outer_wire, edge_dict=None):
        """
        assumes the outer wire and inner wire are oriented away from the face
        they will be creating
        outer - cw
        inner - ccw
        """
        rev_outer = outer_wire  # copy_reversed(outer_wire, edge_dict=edge_dict)
        rev_inner = copy_reversed(inner_wire, edge_dict=edge_dict)

        bdss = brep.BRepExtrema_DistShapeShape(rev_inner, rev_outer)

        support_inner = bdss.SupportOnShape1(1)

        support_outer = bdss.SupportOnShape2(1)
        point_outer = bdss.PointOnShape2(1)

        print(support_inner, support_outer)
        inner_vertex = topods.Vertex(support_inner)
        support_outer = topods.Edge(support_outer)
        both = self(rev_outer, support_outer, point_outer, rev_inner, inner_vertex)
        if edge_dict:
            self.edge_dict_updates(edge_dict)
        return both

    def edge_dict_updates(self, edge_dict, *args):

        for added in self.added:
            edge_dict[hash(added)] = -1

        for old, new in self.replaced:
            old_hash = hash(old)
            for new_edge in new:
                edge_dict[hash(new_edge)] = edge_dict.get(old_hash)

        # for old, new in self.replaced:
        #    edge_dict.remove(hash(old))

    def __call__(self,
                 outer_wire: TopoDS_Shape,
                 support_outer: TopoDS_Edge,
                 new_point: gp.gp_Pnt,
                 inner_wire: TopoDS_Wire,
                 support_inner: TopoDS_Vertex,
                 ) -> TopoDS_Shape:
        """
        Im sure theres a builtin command for this, but i cant find it...

        'outer_wire' is a wire which contains a hole defined by the reverse
        of 'inner_wire'

        Remove the hole by replacing the edge 'support_outer' with a loop that
        includes inner_wire

        Produces a single Wire oriented to its interior face.

        :param outer_wire:
        :param support_outer: Edge on the outer wire
        :param new_point:
        :param inner_wire:
        :param support_inner:
        :return:
        """
        new_vert_out2in = Construct.make_vertex(new_point)
        new_vert_in2out = Construct.make_vertex(new_point)

        v1 = topexp.FirstVertex(support_outer)
        v2 = topexp.LastVertex(support_outer)

        new_edge_a = Construct.make_edge(v1, new_vert_in2out)
        new_edge_d = Construct.make_edge(new_vert_out2in, v2)

        self.replaced.append([support_outer, [new_edge_a, new_edge_d]])

        pt = Common.vertex2pnt(support_inner)
        vert_inner = Construct.make_vertex(pt)

        new_edge_b = Construct.make_edge(new_vert_in2out, support_inner)
        new_edge_c = Construct.make_edge(vert_inner, new_vert_out2in)
        self.added.extend([new_edge_b, new_edge_c])

        # get the edges which are on inner support
        edges_before = []
        edges_after = []

        found = False
        itr = BRepTools_WireExplorer(inner_wire)
        while itr.More():

            if itr.CurrentVertex() == support_inner:
                # the edge leaving the selected support
                found = True

            if found is True:
                edges_after.append(itr.Current())
            else:
                edges_before.append(itr.Current())
            itr.Next()

        new_order = edges_after + edges_before
        last = new_order.pop()
        # -------------------------------------------------
        builder = brep.BRepBuilderAPI_MakeWire()

        builder.Add(new_edge_a)
        builder.Add(new_edge_b)
        for edge in new_order:
            builder.Add(edge)

        last_vert = topexp.FirstVertex(last)
        new_last = Construct.make_edge(last_vert, vert_inner)

        self.replaced.append([last, [new_last]])
        builder.Add(new_last)
        builder.Add(new_edge_c)
        builder.Add(new_edge_d)

        builder.Build()
        replacement = builder.Wire()
        reshaper = brep.BRepTools_ReShape()
        reshaper.Replace(support_outer, replacement)
        return reshaper.Apply(outer_wire)


class SplitFace(bop.ModificationStep):
    def __call__(self,
                 outer_wire: TopoDS_Shape,
                 source_vertex: TopoDS_Vertex,
                 support_outer: TopoDS_Edge,
                 new_point: gp.gp_Pnt,
                 ) -> Tuple[TopoDS_Wire, TopoDS_Wire]:
        """

        'outer_wire' is a wire which contains a hole defined by the reverse
        of 'inner_wire'

        Remove the hole by replacing the edge 'support_outer' with a loop that
        includes inner_wire

        :param outer_wire:
        :param support_outer:
        :param new_pnt:
        :param inner_wire:
        :param support_inner:
        :return:
        """
        vert_new_left = Construct.make_vertex(new_point)
        vert_new_rght = Construct.make_vertex(new_point)

        right_src_vertex = Construct.make_vertex(Common.vertex2pnt(source_vertex))

        vsrc1 = topexp.FirstVertex(support_outer)
        vsrc2 = topexp.LastVertex(support_outer)

        rght_edge_1 = Construct.make_edge(vsrc1, vert_new_rght)
        rght_edge_2 = Construct.make_edge(vert_new_rght, right_src_vertex)

        left_edge_1 = Construct.make_edge(source_vertex, vert_new_left)
        left_edge_2 = Construct.make_edge(vert_new_left, vsrc2)

        self.replaced.append([support_outer, [rght_edge_1, left_edge_2]])
        self.added = [left_edge_1, rght_edge_2]

        #
        build_left = brep.BRepBuilderAPI_MakeWire()
        build_rght= brep.BRepBuilderAPI_MakeWire()

        found = False
        itr = ShapeExtend_WireData(outer_wire)
        z = itr.Index(support_outer)
        print(z)
        itr.SetLast(z)
        print(itr.Index(support_outer))
        N = itr.NbEdges()

        for i in range(1, N):
            edge = itr.Edge(i)
            if found is True:
                build_rght.Add(edge)
            else:
                build_left.Add(edge)
                if topexp.FirstVertex(edge) == source_vertex:
                    print('found ')
                    found = True

        build_left.Add(left_edge_1)
        build_left.Add(left_edge_2)
        build_left.Build()

        build_rght.Add(rght_edge_1)
        build_rght.Add(rght_edge_2)
        build_rght.Build()

        left = build_left.Wire()
        rght = build_rght.Wire()

        return left, rght


class InsertVertex:
    def __init__(self):
        self.edge = None
        self.new_vert = None
        self.e1 = None
        self.e2 = None

    def __call__(self,
                 wire: TopoDS_Shape,
                 sup_edge: Union[TopoDS_Shape, TopoDS_Edge],
                 new_pnt: gp.gp_Pnt):
        if isinstance(sup_edge, TopoDS_Shape):
            sup_edge = topods.Edge(sup_edge)

        self.edge = sup_edge
        self.new_vert = Construct.make_vertex(new_pnt)
        v1 = topexp.FirstVertex(sup_edge)
        v2 = topexp.LastVertex(sup_edge)
        self.e1 = Construct.make_edge(v1, self.new_vert)
        self.e2 = Construct.make_edge(self.new_vert, v2)

        replacement = Construct.make_wire(self.e1, self.e2)
        reshaper = brep.BRepTools_ReShape()
        reshaper.Replace(sup_edge, replacement)
        return reshaper.Apply(wire)

    @property
    def vertex(self):
        return self.new_vert


def add_vertex_to_wire_edge(wire: TopoDS_Shape,
                            sup_edge: Union[TopoDS_Shape, TopoDS_Edge],
                            new_pnt: gp.gp_Pnt):
    """

    :param wire: wire which contains the edge
    :param sup_edge: The edge to be modified
    :param new_pnt:
    :return:
    """
    # if support.ShapeType() == TopAbs_EDGE:
    return wire


def split_hole(inner_wire: TopoDS_Wire, outer_wire: TopoDS_Wire) -> Tuple[TopoDS_Face, TopoDS_Face]:
    """
    Two wires

    **Assumptions**
        - wires are coplanar
        - inner_wire.isInside(outer_wire)

    :param inner_wire:  TopoDS_Wire outer loop
    :param outer_wire: TopoDS_Wire
    :return: new_wire: TopoDS_Wire Th
    """
    _map = TopTools_IndexedDataMapOfShapeListOfShape()
    topexp.MapShapesAndAncestors(inner_wire, TopAbs_VERTEX, TopAbs_EDGE, _map)

    bdss = brep.BRepExtrema_DistShapeShape(inner_wire, outer_wire)
    bdss.Perform()
    # print('numsols', bdss.NbSolution())

    sols_outer = defaultdict(set)
    sols_inner = defaultdict(set)
    sl = defaultdict(set)
    si, sj = 0, 0
    # edge_i, edge_o = None, None

    if bdss.NbSolution() == 1:
        raise Exception('FUCK YOU EDGE CASE!!!!')

    # need an edge-edge to match twice
    for i in range(1, bdss.NbSolution() + 1):
        sup_inner = bdss.SupportOnShape1(i)
        sup_outer = bdss.SupportOnShape2(i)
        hsh_outer = hash(sup_outer)

        if sup_inner.ShapeType() == TopAbs_VERTEX:
            edges_for_v = _map.FindFromKey(sup_inner)
            if edges_for_v.IsEmpty():
                continue
            _itr = TopTools_ListIteratorOfListOfShape(edges_for_v)
            while _itr.More():
                inner_edge = _itr.Value()
                hsh_inner = hash(inner_edge)

                if hsh_outer in sl[hsh_inner]:
                    sols_j = sols_inner[hsh_inner].intersection(sols_outer[hsh_outer])
                    assert len(sols_j) == 1
                    sol_j = sols_j.pop()
                    # edge_i = inner_edge
                    # edge_o = sup_outer
                    si = i
                    sj = sol_j
                    break
                else:
                    sols_inner[hsh_inner].add(i)
                    sl[hsh_inner].add(hsh_outer)

                _itr.Next()

        sols_outer[hsh_outer].add(i)

    if si == 0:
        raise Exception('FUCK!!!!!')

    p1 = bdss.PointOnShape1(si)
    p2 = bdss.PointOnShape2(si)
    p3 = bdss.PointOnShape2(sj)
    p4 = bdss.PointOnShape1(sj)
    cc1 = Wire.create(p1, p2, p3, p4)

    new_face = Construct.make_face(cc1)

    face_out = Construct.make_face(outer_wire)
    face_inn = Construct.make_face(inner_wire)
    face_real = bop.difference(face_out, face_inn)

    face_real = bop.difference(face_real, new_face)

    return next(Topo(face_real).faces()), new_face

