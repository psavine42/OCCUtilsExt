import numpy as np
from lib.OCCUtils import (
     brep, types_lut, topods, Common
)
from OCC.TopoDS import TopoDS_Shape, TopoDS_Edge, TopoDS_Face, TopoDS_Vertex
from lib.OCCUtils.oci import gp
# from OCC.BRep import BRep_Tool
from typing import List, Dict, Union, Tuple, NoReturn, Iterable


ExtremaType = Union[TopoDS_Edge, TopoDS_Face, TopoDS_Vertex]


class Extrema:
    _MP = {'vertex': topods.topods.Vertex,
           'edge': topods.topods.Edge,
           'face': topods.topods.Face
           }

    def __init__(self, bdss: brep.BRepExtrema_DistShapeShape, i, value, s1, s2):
        self.i = i
        self.value = value
        self.support_types = [
            types_lut.brep_extrema_lut[bdss.SupportTypeShape1(i)],
            types_lut.brep_extrema_lut[bdss.SupportTypeShape2(i)]
        ]
        self.bds = bdss
        self.s1, self.s2 = s1, s2

    @classmethod
    def cast(cls, support_type: str, shape) -> ExtremaType:
        return cls._MP[support_type](shape)

    def _cast(self, i) -> ExtremaType:
        return self._MP[self.support_types[i]](self.support_shps[i])

    @property
    def support_shps(self) -> Tuple[TopoDS_Shape, TopoDS_Shape]:
        """ Topo_DS Shapes """
        return self.bds.SupportOnShape1(self.i), self.bds.SupportOnShape2(self.i)

    @property
    def points(self) -> Tuple[gp.gp_Pnt, gp.gp_Pnt]:
        """ Topo_DS Shapes """
        return self.bds.PointOnShape1(self.i), self.bds.PointOnShape2(self.i)

    @property
    def supports(self):
        """ Topo_DS Shapes cast to Vertex, Edge, or Face"""
        return [self._cast(0), self._cast(1)]

    @property
    def params_on_edge(self):
        return [
            self.bds.ParOnEdgeS1(self.i) if self.support_types[0] == 'edge' else None,
            self.bds.ParOnEdgeS2(self.i) if self.support_types[1] == 'edge' else None
        ]

    def _normd_param(self, param, shp):
        # dom = Edge(shp).domain()

        return

    @property
    def edges(self):
        return [self._cast(0) if self.support_types[0] == 'edge' else None,
                self._cast(1) if self.support_types[1] == 'edge' else None]
    @property
    def faces(self):
        return [self._cast(0) if self.support_types[0] == 'face' else None,
                self._cast(1) if self.support_types[1] == 'face' else None]

    def norm_params_edge(self):
        p1, p2 = self.params_on_edge
        if p1 is None and p2 is None:
            return [None, None]

        sup1, sup2 = self.supports
        return

    @property
    def params_on_face(self):
        return [
            self.bds.ParOnFaceS1(self.i) if self.support_types[0] == 'face' else None,
            self.bds.ParOnFaceS2(self.i) if self.support_types[1] == 'face' else None
        ]

    def __repr__(self):
        return 'Extrema[{}]: {} supports: {}'.format(
            self.i, round(self.value, 3), str(self.support_types)
        )

    @classmethod
    def create(cls, shp1, shp2, bdss=None):
        if bdss is None:
            bdss = brep.BRepExtrema_DistShapeShape(shp1, shp2)
        val = bdss.Value()
        for i in range(1, bdss.NbSolution() + 1):
            yield Extrema(bdss, i, val, shp1, shp2)


class Extremas:
    def __init__(self, shp1, shp2):
        self.bdss = brep.BRepExtrema_DistShapeShape(shp1, shp2)
        self.bdss.Perform()
        self.extremas = list(Extrema.create(shp1, shp2, self.bdss))

    @property
    def value(self):
        return self.bdss.Value()

    def points(self, idx=None):
        return [[x.points[idx]] if idx else x.points for x in self.extremas]

    def by_supports(self, stype1, stype2):
        return

    @property
    def faces(self):
        fcs1, fcs2 = [], []
        for ex in self.extremas:
            p1, p2 = ex.faces
            if p1 is not None:
                fcs1.append(p1)
            if p2 is not None:
                fcs2.append(p2)
        return fcs1, fcs2

    @property
    def edges(self):
        fcs1, fcs2 = [], []
        for ex in self.extremas:
            p1, p2 = ex.edges
            if p1 is not None:
                fcs1.append(p1)
            if p2 is not None:
                fcs2.append(p2)
        return fcs1, fcs2
    # ----------------------------------------
    # python
    def __repr__(self):
        return '{}: N={}'.format(self.__class__.__name__, len(self.extremas))

    def __getitem__(self, item):
        return self.extremas[item]

    def __len__(self):
        return self.bdss.NbSolution()

    def __iter__(self):
        for i in self.extremas:
            yield i
