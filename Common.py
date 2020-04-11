#! /usr/bin/python

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
##along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

import random

from OCC.Bnd import Bnd_Box
from OCC.TColgp import (TColgp_HArray1OfPnt,
                        TColgp_Array1OfPnt,
                        TColgp_Array1OfPnt2d,
                        TColgp_Array1OfVec)
from OCC.TColStd import TColStd_HArray1OfBoolean

# --- BREP
from OCC.BRep import BRep_Tool

from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.BRepMesh import BRepMesh_IncrementalMesh
from OCC.BRepBndLib import brepbndlib_Add # , brepbndlib_AddOBB
from OCC.BRepAdaptor import (BRepAdaptor_Curve,
                             BRepAdaptor_HCurve,
                             BRepAdaptor_CompCurve,
                             BRepAdaptor_HCompCurve)
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.BRepExtrema import BRepExtrema_DistShapeShape #, \
                # BRepExtrema_ShapeProximity
from OCC.BRepClass3d import BRepClass3d_SolidClassifier
from OCC.BRepGProp import (brepgprop_LinearProperties,
                           brepgprop_SurfaceProperties,
                           brepgprop_VolumeProperties)
# geom
from OCC.GeomAPI import (GeomAPI_Interpolate,
                         GeomAPI_PointsToBSpline,
                         GeomAPI_IntSS,
                         GeomAPI_ProjectPointOnCurve)

from OCC.GeomAbs import GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.Geom import Geom_Curve, Geom_Plane

from OCC.gp import gp_Pnt, gp_Dir, gp_Ax2, gp_Vec, gp_Trsf, gp_XYZ, gp_Lin

from OCC.TopExp import topexp_CommonVertex, topexp
from OCC.TopAbs import TopAbs_ON, TopAbs_OUT, TopAbs_IN
from OCC.TopoDS import TopoDS_Edge, TopoDS_Shape, TopoDS_Wire, TopoDS_Vertex
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB

from OCC.GProp import GProp_GProps
from OCC.Bnd import Bnd_Box, Bnd_B3d
from OCC import Graphic3d
from OCC.Approx import Approx_Curve3d
from OCC.ProjLib import projlib_Project
from OCC.ShapeFix import ShapeFix_ShapeTolerance
from OCC.ShapeUpgrade import ShapeUpgrade_ShapeDivideContinuity
from OCC.GCPnts import GCPnts_UniformDeflection
from OCC.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCC.IntAna import IntAna_Int3Pln
from OCC.TCollection import TCollection_ExtendedString

# from OCC.BOPInt import
import numpy as np
from typing import Tuple, Optional, List
from .oci.gp import repr3
from . import types_lut
# ===========================================================================
# No PythonOCC dependencies...
# ===========================================================================


class assert_isdone(object):
    """
    raises an assertion error when IsDone() returns false, with the error
    specified in error_statement
    """
    def __init__(self, to_check, error_statement):
        self.to_check = to_check
        self.error_statement = error_statement

    def __enter__(self, ):
        if self.to_check.IsDone():
            pass
        else:
            raise AssertionError(self.error_statement)

    def __exit__(self, assertion_type, value, traceback):
        pass


def reprv(vertex):
    return '<V>:' + str(hash(vertex)) + ' '  + repr3(vertex2pnt(vertex))


def roundlist(li, n_decimals=3):
    return [round(i, n_decimals) for i in li]


# ===========================================================================
# CONSTANTS
# ===========================================================================
TOLERANCE = 1e-6


def get_boundingbox(shape: TopoDS_Shape, tol=TOLERANCE):
    """
    :param shape: TopoDS_Shape such as TopoDS_Face
    :param tol: tolerance
    :return: [xmin, ymin, zmin, xmax, ymax, zmax]
    """
    bbox = Bnd_Box()
    bbox.SetGap(tol)
    brepbndlib_Add(shape, bbox)
    # print(bbox.IsVoid())
    if bbox.IsVoid() is True:
        return None
    # xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    cmin = bbox.CornerMin()
    cmax = bbox.CornerMax()

    xmin, ymin, zmin = cmin.X(), cmin.Y(), cmin.Z()
    xmax, ymax, zmax = cmax.X(), cmax.Y(), cmax.Z()
    return xmin, ymin, zmin, xmax, ymax, zmax


def bbx(shape: TopoDS_Shape, tol=TOLERANCE):
    bbox = Bnd_Box()
    bbox.SetGap(tol)
    brepbndlib_Add(shape, bbox)
    # print(bbox.IsVoid())
    if bbox.IsVoid() is True:
        return None
    return bbox


def aabb(shape, optimal_OBB=True):
    """ return the oriented bounding box of the TopoDS_Shape `shape`
    Parameters
    ----------
    shape : TopoDS_Shape or a subclass such as TopoDS_Face
        the shape to compute the bounding box from
    optimal_OBB : bool, True by default. If set to True, compute the
        optimal (i.e. the smallest oriented bounding box). Optimal OBB is
        a bit longer.
    Returns
    -------
        a list with center, x, y and z sizes
        a shape
    """
    obb = Bnd_OBB()
    if optimal_OBB:
        is_triangulationUsed = True
        is_optimal = True
        is_shapeToleranceUsed = False
        brepbndlib_AddOBB(shape, obb, is_triangulationUsed, is_optimal, is_shapeToleranceUsed)
    else:
        brepbndlib_AddOBB(shape, obb)

    # converts the bounding box to a shape
    aBaryCenter = obb.Center()
    aXDir = obb.XDirection()
    aYDir = obb.YDirection()
    aZDir = obb.ZDirection()
    aHalfX = obb.XHSize()
    aHalfY = obb.YHSize()
    aHalfZ = obb.ZHSize()

    ax = gp_XYZ(aXDir.X(), aXDir.Y(), aXDir.Z())
    ay = gp_XYZ(aYDir.X(), aYDir.Y(), aYDir.Z())
    az = gp_XYZ(aZDir.X(), aZDir.Y(), aZDir.Z())
    p = gp_Pnt(aBaryCenter.X(), aBaryCenter.Y(), aBaryCenter.Z())
    anAxes = gp_Ax2(p, gp_Dir(aZDir), gp_Dir(aXDir))
    anAxes.SetLocation(gp_Pnt(p.XYZ() - ax*aHalfX - ay*aHalfY - az*aHalfZ))
    aBox = BRepPrimAPI_MakeBox(anAxes, 2.0*aHalfX, 2.0*aHalfY, 2.0*aHalfZ).Shape()
    return aBaryCenter, [aHalfX, aHalfY, aHalfZ], aBox


def get_boundingbox_slow(shape, tol=1e-6, use_mesh=True):
    """ return the bounding box of the TopoDS_Shape `shape`
    Parameters
    ----------
    shape : TopoDS_Shape or a subclass such as TopoDS_Face
        the shape to compute the bounding box from
    tol: float
        tolerance of the computed boundingbox
    use_mesh : bool
        a flag that tells whether or not the shape has first to be meshed before the bbox
        computation. This produces more accurate results
    """
    bbox = Bnd_Box()
    bbox.SetGap(tol)
    if use_mesh:
        mesh = BRepMesh_IncrementalMesh()
        mesh.SetParallelDefault(True)
        mesh.SetShape(shape)
        mesh.Perform()
        if not mesh.IsDone():
            raise AssertionError("Mesh not done.")
    brepbndlib_Add(shape, bbox, use_mesh)

    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    return xmin, ymin, zmin, xmax, ymax, zmax


def smooth_pnts(pnts):
    smooth = [pnts[0]]
    for i in range(1, len(pnts)-1):
        prev = pnts[i-1]
        this = pnts[i]
        next_pnt = pnts[i+1]
        pt = (prev+this+next_pnt) / 3.0
        smooth.append(pt)
    smooth.append(pnts[-1])
    return smooth


# ===========================================================================
# Data type utilities
# ===========================================================================
def color(r, g, b):
    return Quantity_Color(r, g, b, Quantity_TOC_RGB)


def to_string(_string):
    return TCollection_ExtendedString(_string)


def to_tcol_(_list, collection_type):
    array = collection_type(1, len(_list)+1)
    for n, i in enumerate(_list):
        array.SetValue(n+1, i)
    return array.GetHandle()


def _Tcol_dim_1(li, _type):
    """function factory for 1-dimensional TCol* types"""
    pts = _type(0, len(li)-1)
    for n, i in enumerate(li):
        pts.SetValue(n, i)
    pts.thisown = False
    return pts


def point_list_to_TColgp_Array1OfPnt(li):
    pts = TColgp_Array1OfPnt(0, len(li)-1)
    for n, i in enumerate(li):
        pts.SetValue(n, i)
    return pts


def point2d_list_to_TColgp_Array1OfPnt2d(li):
    return _Tcol_dim_1(li, TColgp_Array1OfPnt2d)


# ===========================================================================
# --- INTERPOLATION ---
# ===========================================================================
def filter_points_by_distance(list_of_point, distance=0.1):
    """
    get rid of those point that lie within tolerance of a
    consequtive series of points
    """
    tmp = [list_of_point[0]]
    for a in list_of_point[1:]:
        if any([a.IsEqual(i, distance) for i in tmp]):
            continue
        else:
            tmp.append(a)
    return tmp


def points_to_bspline(pnts):
    """
    Points to bspline
    """
    pnts = point_list_to_TColgp_Array1OfPnt(pnts)
    crv = GeomAPI_PointsToBSpline(pnts)
    return crv.Curve()


def interpolate_points_to_spline(list_of_points, start_tangent, end_tangent, filter_pts=True, tolerance=TOLERANCE):
    """
    GeomAPI_Interpolate is buggy: need to use `fix` in order
    to get the right points in...
    """
    def fix(li, _type):
        """function factory for 1-dimensional TCol* types"""
        pts = _type(1, len(li))
        for n, i in enumerate(li):
            pts.SetValue(n+1, i)
        pts.thisown = False
        return pts

    if filter_pts:
        list_of_points = filter_points_by_distance(list_of_points, 0.1)

    fixed_points = fix(list_of_points, TColgp_HArray1OfPnt)
    try:
        interp = GeomAPI_Interpolate(fixed_points.GetHandle(), False, tolerance)
        interp.Load(start_tangent, end_tangent, False)
        interp.Perform()
        if interp.IsDone():
            return interp.Curve()
    except RuntimeError:
        print("Failed to interpolate the shown points")


def interpolate_points_vectors_to_spline(list_of_points, list_of_vectors, vector_mask=None, tolerance=TOLERANCE):
    """
    build a curve from a set of points and vectors
    the vectors describe the tangent vector at the corresponding point
    """
    # GeomAPI_Interpolate is buggy: need to use `fix` in order to
    # get the right points in...
    assert len(list_of_points) == len(list_of_vectors), 'vector and point list not of same length'

    def fix(li, _type):
        """function factory for 1-dimensional TCol* types"""
        pts = _type(1, len(li))
        for n, i in enumerate(li):
            pts.SetValue(n+1, i)
        pts.thisown = False
        return pts

    if vector_mask is not None:
        assert len(vector_mask) == len(list_of_points), 'length vector mask is not of length points list nor []'
    else:
        vector_mask = [True for i in range(len(list_of_points))]

    fixed_mask = fix(vector_mask, TColStd_HArray1OfBoolean)
    fixed_points = fix(list_of_points, TColgp_HArray1OfPnt)
    fixed_vectors = fix(list_of_vectors, TColgp_Array1OfVec)

    try:
        interp = GeomAPI_Interpolate(fixed_points.GetHandle(), False, tolerance)
        interp.Load(fixed_vectors, fixed_mask.GetHandle(), False)
        interp.Perform()
        if interp.IsDone():
            return interp.Curve()
    except RuntimeError:
        # the exception was unclear
        raise RuntimeError('FAILED TO INTERPOLATE THE POINTS')


def interpolate_points_to_spline_no_tangency(list_of_points, filter_pts=True, closed=False, tolerance=TOLERANCE):
    """
    GeomAPI_Interpolate is buggy: need to use `fix`
    in order to get the right points in...
    """
    def fix(li, _type):
        """function factory for 1-dimensional TCol* types"""
        pts = _type(1, len(li))
        for n, i in enumerate(li):
            pts.SetValue(n+1, i)
        pts.thisown = False
        return pts

    if filter_pts:
        list_of_points = filter_points_by_distance(list_of_points, 0.1)

    fixed_points = fix(list_of_points, TColgp_HArray1OfPnt)
    try:
        interp = GeomAPI_Interpolate(fixed_points.GetHandle(), closed, tolerance)
        interp.Perform()
        if interp.IsDone():
            return interp.Curve()

    except RuntimeError:
        # the exception was unclear
        raise RuntimeError('FAILED TO INTERPOLATE THE POINTS')


# ===========================================================================
# --- BUILD PATCHES ---
# ===========================================================================
def common_vertex(edg1, edg2):
    vert = TopoDS_Vertex()
    # topexp.Vertices()
    if topexp.CommonVertex(edg1, edg2, vert):
        return vert
    else:
        raise ValueError('no common vertex found')


def other_vertex(edge: TopoDS_Edge, vert: TopoDS_Vertex):
    if hash(vert) == hash(topexp.FirstVertex(edge)):
        return topexp.LastVertex(edge)
    elif hash(vert) == hash(topexp.LastVertex(edge)):
        return topexp.FirstVertex(edge)
    return None


def midpoint(pnt1: gp_Pnt, pnt2: gp_Pnt) -> gp_Pnt:
    """
    computes the point that lies in the middle between pntA and pntB
    @param pnt1:    gp_Pnt
    @param pnt2:    gp_Pnt
    """
    vec1 = gp_Vec(pnt1.XYZ())
    vec2 = gp_Vec(pnt2.XYZ())
    veccie = (vec1+vec2)/2.
    return gp_Pnt(veccie.XYZ())


def point_in_boundingbox(solid, pnt, tolerance=1e-5):
    """returns True if *pnt* lies in *boundingbox*, False if not
    this is a much speedier test than checking the TopoDS_Solid
    Args:
        solid   TopoDS_Solid
        pnt:    gp_Pnt

    Returns: bool
    """
    bbox = Bnd_Box()
    bbox.SetGap(tolerance)
    brepbndlib_Add(solid, bbox)
    return not bbox.IsOut(pnt)


def point_in_solid(solid, pnt, tolerance=1e-5):
    """returns True if *pnt* lies in *solid*, False if not
    Args:
        solid   TopoDS_Solid
        pnt:    gp_Pnt

    Returns: bool
    """
    _in_solid = BRepClass3d_SolidClassifier(solid, pnt, tolerance)
    # print("State", _in_solid.State())
    if _in_solid.State() == TopAbs_ON:
        return None, 'on'
    elif _in_solid.State() == TopAbs_OUT:
        return False, 'out'
    elif _in_solid.State() == TopAbs_IN:
        return True, 'in'
    else:
        return None, None


def intersection_from_three_planes(planeA, planeB, planeC):
    """
    intersection from 3 planes
    accepts both Geom_Plane and gp_Pln
    @param planeA:
    @param planeB:
    @param planeC:
    """
    planeA = planeA if not hasattr(planeA, 'Pln') else planeA.Pln()
    planeB = planeB if not hasattr(planeB, 'Pln') else planeB.Pln()
    planeC = planeC if not hasattr(planeC, 'Pln') else planeC.Pln()

    intersection_planes = IntAna_Int3Pln(planeA, planeB, planeC)
    pnt = intersection_planes.Value()
    return pnt


def intersect_shape_by_line(shape: TopoDS_Shape, line: gp_Lin, low_parameter=0.0, hi_parameter=float("+inf")):
    """
    finds the intersection of a shape and a line

    :param shape: any TopoDS_*
    :param line: gp_Lin
    :param low_parameter:
    :param hi_parameter:

    :return: a list with a number of tuples that corresponds to the number
    of intersections found
    the tuple contains ( gp_Pnt, TopoDS_Face, u,v,w ), respectively the
    intersection point, the intersecting face
    and the u,v,w parameters of the intersection point
    :raise:
    """
    shape_inter = IntCurvesFace_ShapeIntersector()
    shape_inter.Load(shape, TOLERANCE)
    shape_inter.PerformNearest(line, low_parameter, hi_parameter)

    with assert_isdone(shape_inter, "failed to computer shape / line intersection"):
        return (shape_inter.Pnt(1),
                shape_inter.Face(1),
                shape_inter.UParameter(1),
                shape_inter.VParameter(1),
                shape_inter.WParameter(1))


def normal_vector_from_plane(plane, vec_length=1.) -> gp_Vec:
    """
    returns a vector normal to the plane of length vec_length
    @param plane:
    """
    trns = gp_Vec(plane.Axis().Direction())
    return trns.Normalized() * vec_length


# ===========================================================================
# FIX
# ===========================================================================
def fix_tolerance(shape, tolerance=TOLERANCE):
    ShapeFix_ShapeTolerance().SetTolerance(shape, tolerance)


def fix_continuity(edge, continuity=1):
    su = ShapeUpgrade_ShapeDivideContinuity(edge)
    su.SetBoundaryCriterion(eval('GeomAbs_C'+str(continuity)))
    su.Perform()
    te = st(su.Result())
    return te


def resample_curve_with_uniform_deflection(curve, deflection=0.5, degreeMin=3, degreeMax=8, continuity=GeomAbs_C2, tolerance=1e-4):
    """
    fits a bspline through the samples on `curve`
    @param curve: TopoDS_Wire, TopoDS_Edge, curve
    @param n_samples:
    """
    crv = to_adaptor_3d(curve)
    defl = GCPnts_UniformDeflection(crv, deflection)
    with assert_isdone(defl, 'failed to compute UniformDeflection'):
        print("Number of points:", defl.NbPoints())
    sampled_pnts = [defl.Value(i) for i in range(1, defl.NbPoints())]
    resampled_curve = GeomAPI_PointsToBSpline(point_list_to_TColgp_Array1OfPnt(sampled_pnts), degreeMin, degreeMax, continuity, tolerance)
    return resampled_curve.Curve().GetObject()


# ===========================================================================
# global properties
# ===========================================================================
class GpropsFromShape(object):
    def __init__(self, shape, tolerance=1e-5):
        self.shape = shape
        self.tolerance = tolerance

    def volume(self):
        """returns the volume of a solid"""
        prop = GProp_GProps()
        brepgprop_VolumeProperties(self.shape, prop, self.tolerance)
        return prop

    def surface(self):
        """returns the area of a surface"""
        prop = GProp_GProps()
        brepgprop_SurfaceProperties(self.shape, prop, self.tolerance)
        return prop

    def linear(self):
        """returns the length of a wire or edge"""
        prop = GProp_GProps()
        brepgprop_LinearProperties(self.shape, prop)
        return prop


def area(shape, tolerance=1e-5):
    prop = GProp_GProps()
    brepgprop_SurfaceProperties(shape, prop, tolerance)
    return prop.Mass()


def curve_length(crv):
    """
    get the length from a TopoDS_Edge or TopoDS_Wire
    """
    assert isinstance(crv, (TopoDS_Wire, TopoDS_Edge)), 'either a wire or edge...'
    gprop = GpropsFromShape(crv)
    return gprop.linear().Mass()


# =======================================================================
# Distance
# =======================================================================
def minimum_distance(shp1, shp2):
    """
    compute minimum distance between 2 BREP's
    @param shp1:    any TopoDS_*
    @param shp2:    any TopoDS_*

    @return: minimum distance,
             minimum distance points on shp1
             minimum distance points on shp2
    """
    bdss = BRepExtrema_DistShapeShape(shp1, shp2)
    bdss.Perform()
    with assert_isdone(bdss, 'failed computing minimum distances'):
        min_dist = bdss.Value()
        min_dist_shp1, min_dist_shp2 = [], []
        for i in range(1, bdss.NbSolution()+1):
            min_dist_shp1.append(bdss.PointOnShape1(i))
            min_dist_shp2.append(bdss.PointOnShape2(i))
    return min_dist, min_dist_shp1, min_dist_shp2


def vertex2pnt(vertex) -> gp_Pnt:
    """returns a gp_Pnt from a TopoDS_Vertex"""
    return BRep_Tool.Pnt(vertex)


def vertex_to_tuple(vertex, rnd=3) -> Tuple:
    pnt = BRep_Tool.Pnt(vertex)
    return tuple(map(lambda x: round(x, rnd), [pnt.X(), pnt.Y(), pnt.Z()]))


def pnt_to_np(pnt):
    """returns a gp_Pnt from a TopoDS_Vertex"""
    return np.array([pnt.X(), pnt.Y(), pnt.Z()])


def pnt2list(pnt):
    """returns a gp_Pnt from a TopoDS_Vertex"""
    return [pnt.X(), pnt.Y(), pnt.Z()]

def vertex_to_np(vertex):
    """returns a gp_Pnt from a TopoDS_Vertex"""
    return pnt_to_np(BRep_Tool.Pnt(vertex))


def adapt_edge_to_curve(edg):
    """
    returns a curve adaptor from an edge
    @param edg: TopoDS_Edge
    """
    return BRepAdaptor_Curve(edg)


def adapt_edge_to_hcurve(edg):
    c = BRepAdaptor_HCurve()
    c.ChangeCurve().Initialize(edg)
    return c


def to_adaptor_3d(curveType):
    """
    abstract curve like type into an adaptor3d
    @param curveType:
    """
    if isinstance(curveType, TopoDS_Wire):
        return BRepAdaptor_CompCurve(curveType)
    elif isinstance(curveType, TopoDS_Edge):
        return BRepAdaptor_Curve(curveType)
    elif issubclass(curveType.__class__, Geom_Curve):
        return GeomAdaptor_Curve(curveType.GetHandle())
    elif hasattr(curveType, 'GetObject'):
        _crv = curveType.GetObject()
        if issubclass(_crv.__class__, Geom_Curve):
            return GeomAdaptor_Curve(curveType)
    else:
        raise TypeError('allowed types are Wire, Edge or a subclass of Geom_Curve\nGot a %s' % (curveType.__class__))

from OCC.TopoDS import topods

def nearest_vertex(edge1, edge2):

    # v11 = topexp.FirstVertex(edge1)
    # v12 = topexp.LastVertex(edge1)

    lut = types_lut.brep_extrema_lut
    bdss = BRepExtrema_DistShapeShape(edge1, edge2)
    bdss.Perform()
    with assert_isdone(bdss, 'failed computing minimum distances'):
        for i in range(1, bdss.NbSolution() + 1):
            if ['vertex', 'vertex'] == [lut[bdss.SupportTypeShape1(i)],
                                        lut[bdss.SupportTypeShape2(i)]]:

                return topods.Vertex(bdss.SupportOnShape1(i)), \
                       topods.Vertex(bdss.SupportOnShape2(i))
    # vertex2pnt()


def project_point_on_curve(crv:TopoDS_Edge, pnt:gp_Pnt) -> Tuple[float, gp_Pnt]:
    """

    :param crv: TopoDS_Edge
    :param pnt:
    :return:
    """
    if isinstance(crv, TopoDS_Shape):
        # get the curve handle...
        adaptor = BRepAdaptor_Curve(crv)
        crv = adaptor.Curve().Curve()
        rrr = GeomAPI_ProjectPointOnCurve(pnt, crv)
        return rrr.LowerDistanceParameter(), rrr.NearestPoint()
    else:
        raise NotImplementedError('expected a TopoDS_Edge...')


def project_point_on_plane(plane: Geom_Plane, point: gp_Pnt) -> gp_Pnt:
    """
    project point on plane
    @param plane: Geom_Plane
    @param point: gp_Pnt
    """
    pl = plane.Pln()
    aa, bb = projlib_Project(pl, point).Coord()
    point = plane.Value(aa, bb)
    return point


def wire_to_curve(wire, tolerance=TOLERANCE, order=GeomAbs_C2, max_segment=200, max_order=12):
    """
    a wire can consist of many edges.
    these edges are merged given a tolerance and a curve
    @param wire:
    """
    adap = BRepAdaptor_CompCurve(wire)
    hadap = BRepAdaptor_HCompCurve(adap)
    approx = Approx_Curve3d(hadap.GetHandle(), tolerance, order, max_segment, max_order)
    with assert_isdone(approx, 'not able to compute approximation from wire'):
        return approx.Curve().GetObject()


# ===========================================================================
# --- RANDOMNESS ---
# ===========================================================================
def random_vec():
    x, y, z = [random.uniform(-1, 1) for i in range(3)]
    return gp_Vec(x, y, z)


def random_colored_material_aspect():
    clrs = [i for i in dir(Graphic3d) if i.startswith('Graphic3d_NOM_')]
    color = random.sample(clrs, 1)[0]
    print("color", color)
    return Graphic3d.Graphic3d_MaterialAspect(getattr(Graphic3d, color))


def random_color():
    return color(random.random(), random.random(), random.random())

