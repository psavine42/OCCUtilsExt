# from .oci.brep import BRepClass3d_SolidClassifier
from OCC.BRepClass3d import BRepClass3d_SolidClassifier, BRepClass3d_SClassifier, \
    BRepClass3d_SolidExplorer
from OCC.BRepClass import BRepClass_FaceExplorer, BRepClass_FacePassiveClassifier
from OCC.BRepApprox import BRepApprox_Approx
from OCC.BRepAdaptor import BRepAdaptor_Surface
from OCC.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane, GeomAbs_Circle
from . import Topo, Face, Construct
from OCC.TopoDS import TopoDS_Edge
from OCC.GProp import GProp_GProps
from OCC.BRepGProp import brepgprop
from .oci.gp import gp_Pln, gp_Ax2, gp_Pnt, repr3
import math
from typing import Union
from src.bim import *

def recognize_face(a_face):
    """ Takes a TopoDS shape and tries to identify its nature
    whether it is a plane a cylinder a torus etc.
    if a plane, returns the normal
    if a cylinder, returns the radius
    """
    surf = BRepAdaptor_Surface(a_face, True)
    surf_type = surf.GetType()
    if surf_type == GeomAbs_Plane:
        print("--> plane")
        # look for the properties of the plane
        # first get the related gp_Pln
        gp_pln = surf.Plane()
        location = gp_pln.Location()  # a point of the plane
        normal = gp_pln.Axis().Direction()  # the plane normal
        # then export location and normal to the console output
        print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
        print("--> Normal (global coordinates)", normal.X(), normal.Y(), normal.Z())
        return surf.Plane()
    elif surf_type == GeomAbs_Cylinder:
        print("--> cylinder")
        # look for the properties of the cylinder
        # first get the related gp_Cyl
        gp_cyl = surf.Cylinder()
        location = gp_cyl.Location()  # a point of the axis
        axis = gp_cyl.Axis().Direction()  # the cylinder axis
        # then export location and normal to the console output
        print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
        print("--> Axis (global coordinates)", axis.X(), axis.Y(), axis.Z())

    elif surf_type == GeomAbs_Circle:
        print("--> circle")
        # look for the properties of the cylinder
        # first get the related gp_Cyl
        # gp_cyl = surf.()
        # location = gp_cyl.Location()  # a point of the axis
        # axis = gp_cyl.Axis().Direction()  # the cylinder axis
        # then export location and normal to the console output
        # print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
        # print("--> Axis (global coordinates)", axis.X(), axis.Y(), axis.Z())

    else:
        # TODO there are plenty other type that can be checked
        # see documentation for the BRepAdaptor class
        # https://www.opencascade.com/doc/occt-6.9.1/refman/html/class_b_rep_adaptor___surface.html
        print("not implemented")


def major_faces(geom):
    """
    returns the 2 largest by area faces
    """
    topo = Topo(geom)
    face_arr = []
    for face in topo.faces():
        system = GProp_GProps()
        brepgprop.SurfaceProperties(face, system)
        area = system.Mass()
        face_arr.append((area, face))

    face_arr.sort(reverse=True, key=lambda x: x[0])
    return [x[1] for x in face_arr[:2]]


def pipe_segment_to_edge(compound, uid=None):
    """ if it is two circles with linear face segments,
        the end sections will have N - 2 edges"""
    the_topo = Topo(compound)
    faces = [Face(f) for f in the_topo.faces()]
    faces.sort(reverse=True, key=lambda x: len(x.edges()))
    areas = []
    verts = []
    for f1 in faces[:2]:
        surf = BRepAdaptor_Surface(f1, True)
        surf_type = surf.GetType()
        areas.append(f1.area)

        uv, pnt = f1.mid_point()

        if surf_type == GeomAbs_Plane:
            gp_pln = surf.Plane()
            normal = f1.normal

            gppln2 = gp_Pln(pnt, normal)
            v = MEPConnector(gppln2)
            verts.append(v)

    radius = math.sqrt(sum(areas)) / (2*math.pi)
    edge = MEPCurve(*verts)
    return MEPSegment(edge, radius, uid=uid)



def pipe_segment_to_edge_convll(compound, uid=None):
    """ if it is two circles with linear face segments,
        the end sections will have N - 2 edges"""
    the_topo = Topo(compound)
    faces = [Face(f) for f in the_topo.faces()]
    faces.sort(reverse=True, key=lambda x: len(x.edges()))
    areas = []
    verts = []
    for f1 in faces[:2]:
        surf = BRepAdaptor_Surface(f1, True)
        surf_type = surf.GetType()
        areas.append(f1.area)

        if surf_type == GeomAbs_Plane:
            gp_pln = surf.Plane()
            v = MEPConnector(gp_pln)
            verts.append(v)

    radius = math.sqrt(sum(areas)) / (2*math.pi)
    edge = MEPCurve(*verts)
    return MEPSegment(edge, radius, uid=uid)


def pipe_segment_to_edge_2plnopt(compound, uid=None):
    """ if it is two circles with linear face segments,
        the end sections will have N - 2 edges"""
    the_topo = Topo(compound)
    faces = [Face(f) for f in the_topo.faces()]
    faces.sort(reverse=True, key=lambda x: len(x.edges()))
    areas = []
    verts = []
    for f1 in faces[:2]:
        surf = BRepAdaptor_Surface(f1, True)
        surf_type = surf.GetType()
        areas.append(f1.area)

        if surf_type == GeomAbs_Plane:
            gp_pln = surf.Plane()
            v = MEPConnector(gp_pln)
            verts.append(v)

    radius = math.sqrt(sum(areas)) / (2*math.pi)
    edge = MEPCurve(*verts)
    return MEPSegment(edge, radius, uid=uid)
