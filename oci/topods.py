from OCC.TopoDS import (TopoDS_Wire,
                        TopoDS_Solid,
                        TopoDS_Edge,
                        TopoDS_Vertex,
                        TopoDS_Face,
                        TopoDS_Shape,
                        TopoDS_Builder,
                        topods,
                        TopoDS_Compound,
                        TopoDS_Shell,
                        TopoDS_CompSolid,
                        TopoDS_HShape,
                        TopoDS_Iterator,
                        TopoDS_ListIteratorOfListOfShape,
                        TopoDS_ListNodeOfListOfShape,
                        TopoDS_ListOfShape,
                        TopoDS_TCompound,
                        TopoDS_TCompSolid,
                        TopoDS_TWire,
                        TopoDS_TEdge, TopoDS_TFace, TopoDS_TShape

                        )
from OCC.TopExp import TopExp_Explorer, topexp, topexp_MapShapes
from OCC.TopTools import (
    TopTools_IndexedDataMapOfShapeAddress,
    TopTools_Array1OfListOfShape,
    TopTools_Array1OfShape,
    TopTools_Array2OfShape,
    TopTools_ShapeMapHasher,
    # BASIC
    TopTools_ShapeSet,

    # LIST
    TopTools_ListNodeOfListOfShape,
    TopTools_ListOfShape,
    TopTools_ListIteratorOfListOfShape,

    # MISC
    TopTools_MapOfOrientedShape,
    TopTools_IndexedMapOfOrientedShape,

    # SEQUENCE
    TopTools_SequenceOfShape,

    # SHAPE
    TopTools_MapOfShape,
    TopTools_MapIteratorOfMapOfShape,

    # SHAPE -- INDEX
    TopTools_IndexedMapOfShape,

    # SHAPE -- SHAPE
    TopTools_IndexedDataMapOfShapeShape,
    TopTools_SequenceNodeOfSequenceOfShape,
    TopTools_IndexedDataMapOfShapeListOfShape,

    # SHAPE <---> INTEGER(DATA)
    TopTools_DataMapOfShapeInteger,
    TopTools_DataMapOfShapeListOfInteger,
    TopTools_DataMapOfIntegerListOfShape,
    TopTools_DataMapNodeOfDataMapOfIntegerListOfShape,


    TopTools_IndexedMapNodeOfIndexedMapOfOrientedShape,
    TopTools_OrientedShapeMapHasher,

    TopTools_DataMapIteratorOfDataMapOfOrientedShapeInteger,
    TopTools_DataMapIteratorOfDataMapOfOrientedShapeShape,
    TopTools_DataMapIteratorOfDataMapOfShapeInteger,
    TopTools_DataMapIteratorOfDataMapOfShapeListOfInteger,
    TopTools_DataMapIteratorOfDataMapOfShapeSequenceOfShape,
    TopTools_DataMapIteratorOfDataMapOfShapeListOfShape,
    TopTools_DataMapIteratorOfDataMapOfShapeReal,
    TopTools_DataMapIteratorOfDataMapOfIntegerListOfShape,
    TopTools_DataMapIteratorOfDataMapOfIntegerShape,

    TopTools_IndexedMapNodeOfIndexedMapOfShape,
    TopTools_DataMapIteratorOfDataMapOfShapeShape,
)

from typing import Union

TShape = Union[TopoDS_TCompSolid, TopoDS_TCompound, TopoDS_TWire,
               TopoDS_TEdge, TopoDS_TFace, TopoDS_TShape]


StructShpShp = Union[

    TopTools_IndexedDataMapOfShapeShape,
    TopTools_IndexedDataMapOfShapeListOfShape,

]

StructNode = Union[
    TopTools_SequenceNodeOfSequenceOfShape,
    TopTools_IndexedMapNodeOfIndexedMapOfShape,
    TopTools_DataMapNodeOfDataMapOfIntegerListOfShape,
    TopTools_IndexedMapNodeOfIndexedMapOfOrientedShape

]


StructIterator = Union[
    TopTools_ListIteratorOfListOfShape,
    TopTools_DataMapIteratorOfDataMapOfShapeShape,
    TopTools_DataMapIteratorOfDataMapOfOrientedShapeInteger,
    TopTools_DataMapIteratorOfDataMapOfOrientedShapeShape,
    TopTools_DataMapIteratorOfDataMapOfShapeInteger,
    TopTools_DataMapIteratorOfDataMapOfShapeListOfInteger,
    TopTools_DataMapIteratorOfDataMapOfShapeSequenceOfShape,
    TopTools_DataMapIteratorOfDataMapOfShapeListOfShape,
    TopTools_DataMapIteratorOfDataMapOfShapeReal,
    TopTools_DataMapIteratorOfDataMapOfIntegerListOfShape,
    TopTools_DataMapIteratorOfDataMapOfIntegerShape,
]

