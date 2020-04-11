#!/usr/bin/env python

##Copyright 2008-2013 Jelle Feringa (jelleferinga@gmail.com)
##
##This file is part of pythonOCC.
##
##pythonOCC is free software: you can redistribute it and/or modify
##it under the terms of the GNU Lesser General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
##
##pythonOCC is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU Lesser General Public License for more details.
##
##You should have received a copy of the GNU Lesser General Public License
##along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>

from OCC.TopoDS import TopoDS_Wire
# from OCC.Topo
from .base import BaseObject
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge
from lib.OCCUtils.Construct import make_edge, assert_isdone
from .Topology import Topo
from lib.OCCUtils.edge import Edge, TopoDS_Edge

from OCC.TopExp import topexp
from OCC.TopoDS import TopoDS_Vertex
from OCC.ShapeFix import ShapeFix_Wire
from OCC.ShapeAnalysis import ShapeAnalysis_Wire
from OCC.ShapeExtend import ShapeExtend_WireData
from OCC.BRepTools import BRepTools_WireExplorer

from OCC.TopExp import TopExp_Explorer, topexp_MapShapesAndAncestors
from .oci.topabs import *
from .vertex import Vertex
from lib.OCCUtils.oci import gp


class Wire(TopoDS_Wire, BaseObject):
    def __init__(self, wire, data=None):
        """

        :param wire:
        """
        assert isinstance(wire, TopoDS_Wire), 'need a TopoDS_Wire, got a %s' % wire.__class__
        assert not wire.IsNull()
        super(Wire, self).__init__()
        BaseObject.__init__(self, 'wire', data=data)
        # we need to copy the base shape using the following three
        # lines
        assert self.IsNull()
        self.TShape(wire.TShape())
        self.Location(wire.Location())
        self.Orientation(wire.Orientation())
        assert not self.IsNull()

    def __reversed__(self):
        self.Reverse()

    @classmethod
    def validate(cls, wire, face, tol=1e-6):
        fix = ShapeFix_Wire(wire, face, tol)
        saw = ShapeAnalysis_Wire(wire, face, tol)
        if saw.CheckSmall(tol):
            print('found small')

        if saw.CheckOrder():
            print('found order')
            fix.FixReorder()
            fix.Perform()

        if saw.CheckSelfIntersection():
            print('found inters')

        # if saw.CheckS

        w = fix.Wire()

    @classmethod
    def fix_collinear(cls, wire, tol):
        """

        :return:
        """
        if wire.Reversed() is 1:
            wire.Reverse()

        explorer = BRepTools_WireExplorer(wire)
        prev = None

        while explorer.More():

            curr = explorer.Current()

            if prev is None:
                prev = curr
                explorer.Next()
                continue

            v1 = gp.gp_Vec(topexp.FirstVertex(curr), topexp.LastVertex(curr))
            v2 = gp.gp_Vec(topexp.FirstVertex(prev), topexp.LastVertex(prev))

            if v1.Angle(v2) < tol:
                # remove common vertex
                comm = topexp.CommonVertex(curr, prev)


            # if angle between these is close to 0 or pi, remove the vertex ...

        return

    @classmethod
    def create(cls, *args, close=True):
        builder = BRepBuilderAPI_MakeWire()
        if isinstance(args, (list, tuple)):
            if isinstance(args[0], (gp.gp_Pnt, TopoDS_Vertex)):
                s = 0 if close is True else 1
                for i in range(s, len(args)):
                    edge = make_edge(args[i-1], args[i])
                    builder.Add(edge)

            elif isinstance(args[0], TopoDS_Edge):
                for e in args:
                    builder.Add(e)
        # while not builder.IsDone():

        with assert_isdone(builder, 'f'):
            builder.Build()
            shp = builder.Wire()
            builder.Delete()
            return shp

    def vertices(self):
        return Topo(self).vertices()

    def __repr__(self):
        return 'Wire'





class WireData(ShapeExtend_WireData):
    pass
