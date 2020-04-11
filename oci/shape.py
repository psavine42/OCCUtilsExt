from OCC.ShapeExtend import (
    ShapeExtend_WireData,
    ShapeExtend_DataMapOfShapeListOfMsg,
    ShapeExtend_ComplexCurve,
    ShapeExtend_Explorer,
    ShapeExtend_CompositeSurface,

    ShapeExtend_MsgRegistrator,
    ShapeExtend_BasicMsgRegistrator,
    ShapeExtend_DataMapIteratorOfDataMapOfShapeListOfMsg,
    ShapeExtend_DataMapIteratorOfDataMapOfTransientListOfMsg,
    ShapeExtend_DataMapNodeOfDataMapOfShapeListOfMsg,
    ShapeExtend_DataMapOfTransientListOfMsg,
    ShapeExtend_DataMapNodeOfDataMapOfTransientListOfMsg,
)

from OCC.ShapeFix import (
    ShapeFix_ComposeShell,

    ShapeFix_DataMapIteratorOfDataMapOfShapeBox2d,
    ShapeFix_DataMapNodeOfDataMapOfShapeBox2d,
    ShapeFix_DataMapOfShapeBox2d,
    ShapeFix_EdgeConnect,
    ShapeFix_Edge,
    ShapeFix_Face,
    ShapeFix_Shape,

)
# done

from OCC.ShapeUpgrade import (
    ShapeUpgrade_UnifySameDomain,

    ShapeUpgrade_ShapeConvertToBezier,
    ShapeUpgrade_ConvertCurve2dToBezier,
    ShapeUpgrade_ConvertCurve3dToBezier,
    ShapeUpgrade_ConvertSurfaceToBezierBasis,
    ShapeUpgrade_FixSmallBezierCurves,

    ShapeUpgrade_EdgeDivide,
    ShapeUpgrade_FixSmallCurves,
    ShapeUpgrade_SplitCurve,

    ShapeUpgrade_ClosedEdgeDivide,
    ShapeUpgrade_ClosedFaceDivide,
    ShapeUpgrade_FaceDivide,
    ShapeUpgrade_FaceDivideArea,

    ShapeUpgrade_RemoveInternalWires,
    ShapeUpgrade_RemoveLocations,

    ShapeUpgrade_ShapeDivide,
    ShapeUpgrade_ShapeDivideArea,
    ShapeUpgrade_ShapeDivideAngle,
    ShapeUpgrade_ShapeDivideClosedEdges,
    ShapeUpgrade_ShapeDivideClosed,

    ShapeUpgrade_SplitSurface,

    ShapeUpgrade_SplitCurve2d,
    ShapeUpgrade_SplitCurve2dContinuity,

    ShapeUpgrade_ShellSewing,
    ShapeUpgrade_ShapeDivideContinuity,
    ShapeUpgrade_SplitSurfaceArea,
    ShapeUpgrade_SplitSurfaceAngle,
    ShapeUpgrade_SplitCurve3d,

)

from OCC.ShapeAnalysis import (
    shapeanalysis,
    ShapeAnalysis_CheckSmallFace,
    ShapeAnalysis_Edge,
    ShapeAnalysis_Curve,
    ShapeAnalysis_Wire,
    ShapeAnalysis_Surface,
    ShapeAnalysis_FreeBounds,
    ShapeAnalysis_WireOrder
)

from OCC.ShapeAlgo import (
    ShapeAlgo_ToolContainer
)

from OCC.ShapeBuild import (
    shapebuild,
    ShapeBuild_Edge,
    ShapeBuild_ReShape,
    ShapeBuild_Vertex
)

from OCC.ShapeConstruct import (
    shapeconstruct,
    ShapeConstruct_Curve,

    ShapeConstruct_CompBezierCurves2dToBSplineCurve2d,
    ShapeConstruct_CompBezierCurvesToBSplineCurve,
    ShapeConstruct_MakeTriangulation,
    ShapeConstruct_ProjectCurveOnSurface
)

from OCC.ShapeCustom import (
    ShapeCustom_Curve,
    ShapeCustom_Curve2d,
    ShapeCustom_Surface,
    ShapeCustom_ConvertToBSpline,
    ShapeCustom_DirectModification,
    shapecustom
)