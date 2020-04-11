from OCC.gp import (
    gp_Pln,  gp_Pnt2d, gp_Vec, gp_Lin, gp_Pnt,

    gp_Dir,
    gp_Ax2, gp_Ax1,  gp_Ax2d, gp_Ax3, gp_Circ, gp_Dir2d, gp_Lin2d,

    gp_XYZ, gp_Rotation, gp_Scale, gp_Trsf, gp_Trsf2d, gp_CompoundTrsf,

    gp_Origin, gp_Origin2d, gp_Hypr, gp_Elips, gp_Elips2d,
    gp_Parab, gp_Parab2d, gp_Cone, gp_Cylinder, gp_Torus, gp_Sphere,

    gp_XY, gp_DX, gp_DZ, gp_DY, gp_DY2d, gp_DX2d
    )

from typing import Union

GP3d = Union[gp_Pln, gp_Pnt, gp_Lin, gp_Circ]
GP2d = Union[gp_Pnt2d, gp_Lin2d]

GPEdge = Union[gp_Lin, gp_Circ, gp_Hypr, gp_Elips, gp_Parab]
GPFace = Union[gp_Cone, gp_Cylinder, gp_Torus, gp_Pln, gp_Sphere]


class GPBase:
    def __init__(self, base=None):
        # self._base = base
        pass

    # def __getattr__(self, item):
    #     if hasattr(self._base, item):
    #         return getattr(self._base, item)
    #     elif item in self.__dict__:
    #         return self.__dict__[item]

def repr3c(self):
    return '({} {} {})'.format(
        round(self.X(), 3), round(self.Y(), 3), round(self.Z(), 3)
    )

def repr3(self):
    return '<{}>({} {} {})'.format(
        self.__class__.__name__, round(self.X(), 3), round(self.Y(), 3), round(self.Z(), 3)
    )


def repr6(self):
    return '<{}>({} {} {})'.format(
        self.__class__.__name__, round(self.X(), 3), round(self.Y(), 3), round(self.Z(), 3)
    )


class GP3(GPBase):
    def __repr__(self):
        return repr3(self)


class Pnt(gp_Pnt, GP3):
    def __init__(self, pnt):
        gp_Pnt.__init__(self, pnt)
        GP3.__init__(self)
        # self.SetXYZ(pnt.XYZ())

    def __repr__(self):
        return repr3(self)


class Dir(gp_Dir, GP3):
    def __init__(self, args):
        GP3.__init__(self)
        gp_Dir.__init__(self, args)

    def __repr__(self):
        return repr3(self)


class XYZ(gp_XYZ, GP3):
    def __init__(self, args):
        gp_XYZ.__init__(self, args)
        GP3.__init__(self)

    def __repr__(self):
        return repr3(self)


class Lin(gp_Lin, GP3):
    def __init__(self, *args):
        gp_Lin.__init__(self, *args)
        GP3.__init__(self)

    def __repr__(self):
        loc = self.Location()
        vec = self.Direction()
        return '<{}>p:({}, {}, {}) d:({}, {}, {})'.format(
            self.__class__.__name__, loc.X(), loc.Y(), loc.Z(),
            vec.X(), vec.Y(), vec.Z()
        )


class Pln(gp_Pln, GP3):
    def __init__(self, *args):
        gp_Pln.__init__(self, *args)
        GP3.__init__(self)

    def __repr__(self):
        loc = self.Location()
        vec = self.Axis().Direction()
        return '<{}>p:({}, {}, {}) d:({}, {}, {})'.format(
            self.__class__.__name__, loc.X(), loc.Y(), loc.Z(),
            vec.X(), vec.Y(), vec.Z()
        )

