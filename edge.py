##Copyright 2008-2015 Jelle Feringa (jelleferinga@gmail.com)
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
from OCC.GCPnts import GCPnts_UniformAbscissa, GCPnts_AbscissaPoint
from OCC.gp import gp_Vec, gp_Dir, gp_Pnt

from OCC.TopExp import topexp
from OCC.TopoDS import TopoDS_Edge, TopoDS_Vertex, TopoDS_Face

from OCC.GeomLProp import GeomLProp_CurveTool
from OCC.Geom import Geom_OffsetCurve, Geom_TrimmedCurve
from OCC.GeomLib import geomlib
from OCC.GeomAPI import GeomAPI_ProjectPointOnCurve

from OCC.ShapeAnalysis import ShapeAnalysis_Edge

from OCC.BRep import BRep_Tool, BRep_Tool_Continuity
from OCC.BRepIntCurveSurface import BRepIntCurveSurface_Inter
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.BRepLProp import BRepLProp_CLProps
from OCC.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_HCurve

# high-level
from .Common import vertex2pnt, minimum_distance, assert_isdone, fix_continuity
from .Construct import make_edge
from .types_lut import geom_lut
from .base import BaseObject
from typing import List, Union, Optional


class IntersectCurve(object):
    def __init__(self, instance):
        self.instance = instance

    def intersect(self, other, tolerance=1e-2):
        """Intersect self with a point, curve, edge, face, solid
        method wraps dealing with the various topologies
        """
        if isinstance(other, TopoDS_Face):
            face_curve_intersect = BRepIntCurveSurface_Inter()
            face_curve_intersect.Init(other, self.instance.adaptor.Curve(), tolerance)
            pnts = []
            while face_curve_intersect.More():
                next(face_curve_intersect)
                pnts.append(face_curve_intersect.Pnt())
            return pnts


class DiffGeomCurve(object):
    def __init__(self, instance):
        self.instance = instance
        self._local_props = BRepLProp_CLProps(self.instance.adaptor, 2, self.instance.tolerance)

    @property
    def _curvature(self):
        return self._local_props

    def radius(self, u):
        """returns the radius at u
        """
        # NOT SO SURE IF THIS IS THE SAME THING!!!
        self._curvature.SetParameter(u)
        pnt = gp_Pnt()
        self._curvature.CentreOfCurvature(pnt)
        return pnt

    def curvature(self, u):
        # ugly
        self._curvature.SetParameter(u)
        return self._curvature.Curvature()

    def tangent(self, u):
        """sets or gets ( iff vector ) the tangency at the u parameter
        tangency can be constrained so when setting the tangency,
        you're constrainting it in fact
        """
        self._curvature.SetParameter(u)
        if self._curvature.IsTangentDefined():
            ddd = gp_Dir()
            self._curvature.Tangent(ddd)
            return ddd
        else:
            raise ValueError('no tangent defined')

    def normal(self, u):
        """returns the normal at u

        computes the main normal if no normal is found
        see:
        www.opencascade.org/org/forum/thread_645+&cd=10&hl=nl&ct=clnk&gl=nl
        """
        try:
            self._curvature.SetParameter(u)
            a_dir = gp_Dir()
            self._curvature.Normal(a_dir)
            return a_dir
        except:
            raise ValueError('no normal was found')

    def derivative(self, u, n):
        """
        returns n derivatives at parameter b
        """
        self._curvature.SetParameter(u)
        deriv = {1: self._curvature.D1,
                 2: self._curvature.D2,
                 3: self._curvature.D3,
                 }
        try:
            return deriv[n]
        except KeyError:
            raise AssertionError('n of derivative is one of [1,2,3]')

    def points_from_tangential_deflection(self):
        pass


# ===========================================================================
#    Curve.Construct
# ===========================================================================
class ConstructFromCurve():
    def __init__(self, instance):
        self.instance = instance

    def make_offset(self, offset, vec):
        """
        returns an offsetted curve
        @param offset: the distance between self.crv and the curve to offset
        @param vec:    offset direction
        """
        return Geom_OffsetCurve(self.instance.h_crv, offset, vec)


class Edge(TopoDS_Edge, BaseObject):
    def __init__(self, edge, data=None):
        assert isinstance(edge, TopoDS_Edge), 'need a TopoDS_Edge, got a %s' % edge.__class__
        assert not edge.IsNull()
        super(Edge, self).__init__()
        BaseObject.__init__(self, 'edge', data=data)
        # we need to copy the base shape using the following three
        # lines
        assert self.IsNull()
        self.TShape(edge.TShape())
        self.Location(edge.Location())
        self.Orientation(edge.Orientation())
        assert not self.IsNull()

        # tracking state
        self._local_properties_init = False
        self._curvature_init = False
        self._geometry_lookup_init = False
        self._curve_handle = None
        self._curve = None
        self._adaptor = None
        self._adaptor_handle = None

        # instantiating cooperative classes
        # cooperative classes are distinct through CamelCaps from
        # normal method -> pep8
        self.DiffGeom = DiffGeomCurve(self)
        self.Intersect = IntersectCurve(self)
        self.Construct = ConstructFromCurve(self)

        # GeomLProp object
        self._curvature = None

    def Shape(self):
        return list(Topo(self).edges()).pop()

    def is_closed(self):
        return self.adaptor.IsClosed()

    def is_periodic(self):
        return self.adaptor.IsPeriodic()

    def is_rational(self):
        return self.adaptor.IsRational()

    def continuity(self):
        return self.adaptor.Continuity

    def degree(self):
        if 'line' in self.type:
            return 1
        elif 'curve' in self.type:
            return self.adaptor.Degree()
        else:
            # hyperbola, parabola, circle
            return 2

    def nb_knots(self):
        return self.adaptor.NbKnots()

    def nb_poles(self):
        return self.adaptor.NbPoles()

    @property
    def curve(self):
        if self._curve is not None and not self.is_dirty:
            pass
        else:
            self._curve_handle = BRep_Tool().Curve(self)[0]
            self._curve = self._curve_handle.GetObject()
        return self._curve

    @property
    def curve_handle(self):
        if self._curve_handle is not None and not self.is_dirty:
            return self._curve_handle
        else:
            return None

    @property
    def adaptor(self):
        if self._adaptor is not None and not self.is_dirty:
            pass
        else:
            self._adaptor = BRepAdaptor_Curve(self)
            self._adaptor_handle = BRepAdaptor_HCurve(self._adaptor)
        return self._adaptor

    @property
    def adaptor_handle(self):
        if self._adaptor_handle is not None and not self.is_dirty:
            pass
        else:
            self.adaptor
        return self._adaptor_handle

    @property
    def geom_curve_handle(self):
        """
        :return: Handle_Geom_Curve adapted from `self`
        """
        if self._adaptor_handle is not None and not self.is_dirty:
            return self._adaptor.Curve().Curve()
        else:
            return None

    @property
    def type(self):
        return geom_lut[self.adaptor.Curve().GetType()]

    def pcurve(self, face):
        """
        computes the 2d parametric spline that lies on the surface of the face
        :return: Geom2d_Curve, u, v
        """
        crv, u, v = BRep_Tool().CurveOnSurface(self, face)
        return crv.GetObject(), u, v

    def _local_properties(self):
        self._lprops_curve_tool = GeomLProp_CurveTool()
        self._local_properties_init = True

    def domain(self):
        """returns the u,v domain of the curve"""
        return self.adaptor.FirstParameter(), self.adaptor.LastParameter()

    # ===========================================================================
    #    Curve.GlobalProperties
    def length(self, lbound=None, ubound=None, tolerance=1e-5):
        """returns the curve length
        if either lbound | ubound | both are given, than the length
        of the curve will be measured over that interval
        """
        _min, _max = self.domain()
        if _min < self.adaptor.FirstParameter():
            raise ValueError('the lbound argument is lower than the first parameter of the curve: %s ' % (self.adaptor.FirstParameter()))
        if _max > self.adaptor.LastParameter():
            raise ValueError('the ubound argument is greater than the last parameter of the curve: %s ' % (self.adaptor.LastParameter()))

        lbound = _min if lbound is None else lbound
        ubound = _max if ubound is None else ubound
        return GCPnts_AbscissaPoint().Length(self.adaptor, lbound, ubound, tolerance)

    # ===========================================================================
    #    Curve.modify
    def trim(self, lbound, ubound):
        """
        trim the curve
        @param lbound:
        @param ubound:
        """
        a, b = sorted([lbound, ubound])
        tr = Geom_TrimmedCurve(self.adaptor.Curve().Curve(), a, b).GetHandle()
        return Edge(make_edge(tr))

    def extend_by_point(self, pnt, degree=3, beginning=True):
        """extends the curve to point

        does not extend if the degree of self.curve > 3
        @param pnt:
        @param degree:
        @param beginning:
        """
        if self.degree > 3:
            raise ValueError('to extend you self.curve should be <= 3, is %s' % (self.degree))
        return geomlib.ExtendCurveToPoint(self.curve, pnt, degree, beginning)

    # ===========================================================================
    #    Curve.
    def closest(self, other):
        return minimum_distance(self, other)

    def project_vertex(self, pnt_or_vertex):
        """ returns the closest orthogonal project on `pnt` on edge
        """
        if isinstance(pnt_or_vertex, TopoDS_Vertex):
            pnt_or_vertex = vertex2pnt(pnt_or_vertex)

        poc = GeomAPI_ProjectPointOnCurve(pnt_or_vertex, self.curve_handle)
        return poc.LowerDistanceParameter(), poc.NearestPoint()

    def distance_on_curve(self, distance, close_parameter, estimate_parameter):
        """returns the parameter if there is a parameter
        on the curve with a distance length from u
        raises OutOfBoundary if no such parameter exists
        """
        gcpa = GCPnts_AbscissaPoint(self.adaptor, distance, close_parameter, estimate_parameter, 1e-5)
        with assert_isdone(gcpa, 'couldnt compute distance on curve'):
            return gcpa.Parameter()

    def mid_point(self):
        """
        :return: the parameter at the mid point of the curve, and
        its corresponding gp_Pnt
        """
        _min, _max = self.domain()
        _mid = (_min+_max) / 2.
        return _mid, self.adaptor.Value(_mid)

    def divide_by_number_of_points(self, n_pts, lbound=None, ubound=None):
        """returns a nested list of parameters and points on the edge
        at the requested interval [(param, gp_Pnt),...]
        """
        _lbound, _ubound = self.domain()
        if lbound:
            _lbound = lbound
        elif ubound:
            _ubound = ubound

        # minimally two points or a Standard_ConstructionError is raised
        if n_pts <= 1:
            n_pts = 2

        try:
            npts = GCPnts_UniformAbscissa(self.adaptor, n_pts, _lbound, _ubound)
        except:
            print("Warning : GCPnts_UniformAbscissa failed")
        if npts.IsDone():
            tmp = []
            for i in range(1, npts.NbPoints()+1):
                param = npts.Parameter(i)
                pnt = self.adaptor.Value(param)
                tmp.append((param, pnt))
            return tmp
        else:
            return None

    def __eq__(self, other):
        if hasattr(other, 'topo'):
            return self.IsEqual(other)
        else:
            return self.IsEqual(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def first_vertex(self) -> TopoDS_Vertex:
        return topexp.FirstVertex(self)

    def last_vertex(self) -> TopoDS_Vertex:
        return topexp.LastVertex(self)

    def get_common_vertex(self, edge: TopoDS_Edge) -> Union[TopoDS_Vertex, None]:
        return self.common_vertex(self, edge)

    @classmethod
    def common_vertex(cls, edge1, edge2):
        vert = TopoDS_Vertex()
        if topexp.CommonVertex(edge1, edge2, vert):
            return vert
        else:
            return None

    def as_vec(self) -> gp_Vec:
        if self.is_line():
            first, last = map(vertex2pnt, [self.first_vertex(), self.last_vertex()])
            return gp_Vec(first, last)
        else:
            raise ValueError("edge is not a line, hence no meaningful vector can be returned")

    @classmethod
    def create(cls, *args):
        """
        https://www.opencascade.com/doc/occt-7.3.0/refman/html/class_b_rep_builder_a_p_i___make_edge.html
        :type V1: TopoDS_Vertex &
        :type V2: TopoDS_Vertex &

        :type P1: gp_Pnt
        :type P2: gp_Pnt

        :type L: gp_Lin
        :type p1: (Optional) float or gp_Pnt or TopoDS_Vertex
        :type p2: (Optional) float or gp_Pnt or TopoDS_Vertex

        :type L: gp_Circ
        :type p1: (Optional) float or gp_Pnt or TopoDS_Vertex
        :type p2: (Optional) float or gp_Pnt or TopoDS_Vertex

        :type L: gp_Elips
        :type p1: (Optional) float or gp_Pnt or TopoDS_Vertex
        :type p2: (Optional) float or gp_Pnt or TopoDS_Vertex

        :type L: gp_Hypr
        :type p1: (Optional) float or gp_Pnt or TopoDS_Vertex
        :type p2: (Optional) float or gp_Pnt or TopoDS_Vertex

        :type L: gp_Parab
        :type p1: (Optional) float or gp_Pnt or TopoDS_Vertex
        :type p2: (Optional) float or gp_Pnt or TopoDS_Vertex

        :type L: Handle_Geom_Curve &
        :type p1: (Optional) float or gp_Pnt or TopoDS_Vertex
        :type p2: (Optional) float or gp_Pnt or TopoDS_Vertex

        :type L: Handle_Geom_Curve &
        :type p1: (Optional) gp_Pnt or TopoDS_Vertex
        :type p2: (Optional) gp_Pnt or TopoDS_Vertex
        :type p1: (Optional) float
        :type p2: (Optional) float

        :type L: Handle_Geom2d_Curve &
        :type S: Handle_Geom_Surface &
        :type p1: (Optional) float or gp_Pnt or TopoDS_Vertex
        :type p2: (Optional) float or gp_Pnt or TopoDS_Vertex

        :type L: Handle_Geom2d_Curve
        :type S: Handle_Geom_Surface &
        :type P1: (Optional) gp_Pnt or TopoDS_Vertex
        :type P2: (Optional) gp_Pnt or TopoDS_Vertex
        :type p1: (Optional) float
        :type p2: (Optional) float

        * The general method to directly create an edge is to give -
        a 3D curve C as the support (geometric domain) of the edge, -
        two vertices V1 and V2 to limit the curve
        (definition of the restriction of the edge),
        and - two real values p1 and p2 which are the
        parameters for the vertices V1 and V2 on the curve.
        The curve may be defined as a 2d curve in the parametric
        space of a surface: a pcurve.
        The surface on which the edge is built is then kept at the level
        of the edge.
        The default tolerance will be associated with this edge.
        """
        edge = BRepBuilderAPI_MakeEdge(*args)
        with assert_isdone(edge, 'failed to produce edge'):
            result = edge.Edge()
            edge.Delete()
            return result

    # ===========================================================================
    #  Curve.
    def parameter_to_point(self, u):
        """returns the coordinate at parameter u
        """
        return self.adaptor.Value(u)

    def fix_continuity(self, continuity):
        """
        splits an edge to achieve a level of continuity
        :param continuity: GeomAbs_C*
        """
        return fix_continuity(self, continuity)

    def continuity_from_faces(self, f1, f2):
        return BRep_Tool_Continuity(self, f1, f2)

    def is_line(self):
        """checks if the curve is planar
        """
        if self.nb_knots() == 2 and self.nb_poles() == 2:
            return True
        else:
            return False

    def is_seam(self, face):
        """
        :return: True if the edge has two pcurves on one surface
        ( in the case of a sphere for example... )
        """
        sae = ShapeAnalysis_Edge()
        return sae.IsSeam(self, face)

    def is_edge_on_face(self, face):
        """checks whether 'self' lies on a surface or a face
        """
        shp = ShapeAnalysis_Edge()
        return shp.HasPCurve(self, face)

    # ===========================================================================
    #    Curve.graphic
    def show(self):
        """
        poles, knots, should render all slightly different.
        here's how...

        http://www.opencascade.org/org/forum/thread_1125/
        """
        super(Edge, self).show()


if __name__ == '__main__':
    from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
    from .Topology import Topo
    b = BRepPrimAPI_MakeBox(10, 20, 30).Shape()
    t = Topo(b)
    ed = next(t.edges())
    my_e = Edge(ed)
    print(my_e.tolerance)
