from OCC.TopLoc import (
    TopLoc_IndexedMapOfLocation,
    TopLoc_Datum3D,
    TopLoc_ItemLocation,
    TopLoc_Location,
    TopLoc_MapIteratorOfMapOfLocation,
    TopLoc_IndexedMapNodeOfIndexedMapOfLocation,
    TopLoc_MapLocationHasher,
    TopLoc_MapOfLocation,
    TopLoc_SListOfItemLocation,
    TopLoc_SListNodeOfItemLocation,
    TopLoc_StdMapNodeOfMapOfLocation,
)

from OCC.LocOpe import (
    LocOpe_SplitShape,
    LocOpe_BuildShape,

)

# https://www.opencascade.com/doc/occt-7.0.0/overview/html/occt_user_guides__boolean_operations.html#occt_algorithms_4
from OCC.BOPDS import (
    BOPDS_DS,
    BOPDS_CommonBlock,
    BOPDS_Curve,
    BOPDS_Tools
)