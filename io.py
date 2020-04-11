from OCC.TopoDS import TopoDS_Shape
from OCC.BRepMesh import BRepMesh_IncrementalMesh
from OCC.StlAPI import stlapi_Read, StlAPI_Writer
from OCC.BRep import BRep_Builder
from OCC.gp import gp_Pnt, gp_Dir, gp_Pnt2d
from OCC.Bnd import Bnd_Box2d
from OCC.TopoDS import TopoDS_Compound
from OCC.IGESControl import IGESControl_Reader, IGESControl_Writer
from OCC.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCC.Interface import Interface_Static_SetCVal
from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
from OCC.TDocStd import TDocStd_Document
from OCC.XCAFDoc import (XCAFDoc_DocumentTool_ShapeTool,
                              XCAFDoc_DocumentTool_ColorTool)
from OCC.STEPCAFControl import STEPCAFControl_Reader
from OCC.TDF import TDF_LabelSequence, TDF_Label
from OCC.TCollection import TCollection_ExtendedString
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.TopLoc import TopLoc_Location
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
import os
from lib.OCCUtils import Construct, Topo, Common
from OCC.gp import gp_Pnt
# from OCC.Extend.TopologyUtils import (discretize_edge, get_sorted_hlr_edges,
#                                       list_of_shapes_to_compound)

from lib.OCCUtils.Construct import \
    make_vertex, make_edge, make_wirex, make_face, \
    make_sewed, make_solid
from typing import List



def read_step_file(filename, as_compound=True, verbosity=True):
    """ read the STEP file and returns a compound
    filename: the file path
    verbosity: optional, False by default.
    as_compound: True by default. If there are more than one shape at root,
    gather all shapes into one compound. Otherwise returns a list of shapes.
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError("%s not found." % filename)

    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)

    if status == IFSelect_RetDone:  # check status
        if verbosity:
            failsonly = False
            step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
            step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)
        transfer_result = step_reader.TransferRoots()
        if not transfer_result:
            raise AssertionError("Transfer failed.")
        _nbs = step_reader.NbShapes()
        if _nbs == 0:
            raise AssertionError("No shape to transfer.")
        elif _nbs == 1:  # most cases
            return step_reader.Shape(1)
        elif _nbs > 1:
            print("Number of shapes:", _nbs)
            shps = []
            # loop over root shapes
            for k in range(1, _nbs + 1):
                new_shp = step_reader.Shape(k)
                if not new_shp.IsNull():
                    shps.append(new_shp)
            if as_compound:
                builder = BRep_Builder()
                compound = TopoDS_Compound()
                builder.MakeCompound(compound)
                for s in shps:
                    builder.Add(compound, s)
                # shps = compound
                # compound, result = list_of_shapes_to_compound(shps)
                # if not result:
                #    print("Warning: all shapes were not added to the compound")
                return compound
            else:
                print("Warning, returns a list of shapes.")
                return shps
    else:
        raise AssertionError("Error: can't read file.")
    return None


def write_step_file(a_shape, filename, application_protocol="AP203"):
    """ exports a shape to a STEP file
    a_shape: the topods_shape to export (a compound, a solid etc.)
    filename: the filename
    application protocol: "AP203" or "AP214IS" or "AP242DIS"
    """
    # a few checks
    if a_shape.IsNull():
        raise AssertionError("Shape %s is null." % a_shape)
    if application_protocol not in ["AP203", "AP214IS", "AP242DIS"]:
        raise AssertionError("application_protocol must be either AP203 or AP214IS. You passed %s." % application_protocol)
    if os.path.isfile(filename):
        print("Warning: %s file already exists and will be replaced" % filename)
    # creates and initialise the step exporter
    step_writer = STEPControl_Writer()
    Interface_Static_SetCVal("write.step.schema", application_protocol)

    # transfer shapes and write file
    step_writer.Transfer(a_shape, STEPControl_AsIs)
    status = step_writer.Write(filename)

    if not status == IFSelect_RetDone:
        raise IOError("Error while writing shape to STEP file.")
    if not os.path.isfile(filename):
        raise IOError("File %s was not saved to filesystem." % filename)


def make_face_holes(*wires):
    if len(wires) == 1:
        return make_face(*wires)
    # print('wires', len(wires))
    faces = [make_face(w) for w in wires]
    shell = make_sewed(*faces)
    return shell


class TopoIO:
    TYPES = ['verts', 'edges', 'wires', 'faces', 'shells', 'solids']
    fns = [make_vertex, make_edge, make_wirex, make_face_holes, make_sewed, make_solid]

    def __init__(self):
        self._seen = [set() for _ in range(6)]
        self._built = [{} for _ in range(6)]

    @classmethod
    def dict_to_topo(cls, tdict, uid=None, solid=False):
        _built = [{} for _ in range(6)]
        used = [set() for _ in range(6)]

        try:

            for k, v in tdict['verts'].items():
                _built[0][k] = Construct.make_vertex(gp_Pnt(*v))

            for i in range(1, 6):
                for k, v in tdict[cls.TYPES[i]].items():
                    try:
                        _built[i][k] = cls.fns[i](*[_built[i-1][x] for x in v if x in _built[i-1]])
                        used[i-1] = used[i-1].union(set(v))
                    except:
                        # print('f', cls.TYPES[i], v)
                        pass
            res = []
            for i in range(6):
                typ = cls.TYPES[i]
                in_keys = set(tdict[typ].keys())
                not_seen = in_keys.difference(used[i])
                res += [_built[i][k] for k in not_seen if k in _built[i]]

            if len(res) == 1:
                # Clicked view:  2 523273845 <Topo>(v:32 e:40 w:16, f:13, c:0, s:1)
                return res.pop()

            if solid is True:
                tops = Topo(Construct.compound(res))
                ls = list(tops.shells())
                ls.sort(key=lambda x: Topo(x).number_of_faces())
                last = ls.pop()
                return make_solid(last)
            else:
                return Construct.compound(res)
        except Exception as e:
            print('error', uid, e)
            return None

    @staticmethod
    def topo_to_dict(shape):
        if not isinstance(shape, Topo):
            topo = Topo(shape)
        else:
            topo = shape
        verts = {}
        edges = {}
        wires = {}
        faces = {}
        shells = {}
        solids = {}

        for v in topo.vertices():
            verts[hash(v)] = Common.vertex_to_tuple(v, rnd=10)

        for edge in topo.edges():
            edges[hash(edge)] = [hash(x) for x in Topo(edge).vertices()]

        for wire in topo.wires():
            wires[hash(wire)] = [hash(x) for x in Topo(wire).edges()]

        for face in topo.faces():
            faces[hash(face)] = [hash(x) for x in Topo(face).wires()]

        for shell in topo.shells():
            shells[hash(shell)] = [hash(x) for x in Topo(shell).faces()]

        for solid in topo.solids():
            solids[hash(solid)] = [hash(x) for x in Topo(solid).shells()]

        return {'verts': verts,
                'edges': edges,
                'wires': wires,
                'faces': faces,
                'shells': shells,
                'solids': solids
                }




