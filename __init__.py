from . import types_lut
from . import Common
from . import Construct
from . import edge
from . import face
from . import vertex
from . import Topology
from . import solid
from . import shell
from . import Iteration
from . import wire

from .types_lut import ShapeToTopology
from .Common import get_boundingbox
from .Topology import Topo, TopoMap, ShapeSet, ShapeDict
from .face import Face
from .edge import Edge
from .vertex import Vertex
from .wire import Wire
from .solid import Shell, Solid

# import sanely wrapped OCC stuff
from .oci import brep
from OCC.gp import *
from .oci import topods

from .extrema import Extrema, Extremas

from OCC.TopoDS import (
    TopoDS_Vertex,
    TopoDS_Shape,
    TopoDS_Compound,
    TopoDS_Face,
    TopoDS_Shell,
    TopoDS_Wire,
    TopoDS_Solid,
    TopoDS_Edge,
    TopoDS_CompSolid
)

from .io import TopoIO