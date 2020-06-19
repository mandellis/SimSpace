//! ----------------
//! custom includes
//! ----------------
#include "graphicstools.h"
#include "mydefines.h"

//! ----
//! OCC
//! ----
#include <AIS_ColoredShape.hxx>

graphicsTools::graphicsTools()
{
    ;
}

//! -------------------------------
//! function: getModelFeatureColor
//! details:
//! -------------------------------
Quantity_NameOfColor graphicsTools::getModelFeatureColor(SimulationNodeClass::nodeType theType)
{
    Quantity_NameOfColor color;
    switch(theType)
    {
    case SimulationNodeClass::nodeType_namedSelectionElement:
    case SimulationNodeClass::nodeType_namedSelectionGeometry:
    case SimulationNodeClass::nodeType_namedSelectionNodal:
        color = NAMED_SELECTION_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
        color = FIXED_SUPPORT_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
        color = CYLINDRICAL_SUPPORT_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
        color = FRICTIONLESS_SUPPORT_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
        color = DISPLACEMENT_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
        color = REMOTE_DISPLACEMENT_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        color = FORCE_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        color = REMOTE_FORCE_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        color = MOMENT_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        color = PRESSURE_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        color = THERMALCONDITION_COLOR;
        break;
    case SimulationNodeClass::nodeType_meshEdgeSize:
    case SimulationNodeClass::nodeType_meshFaceSize:
    case SimulationNodeClass::nodeType_meshGrading:
    case SimulationNodeClass::nodeType_meshBodyMeshControl:
    case SimulationNodeClass::nodeType_meshBodyMeshMethod:
    case SimulationNodeClass::nodeType_meshMethod:
    case SimulationNodeClass::nodeType_meshMaxElementSize:
    case SimulationNodeClass::nodeType_meshMinElementSize:
    case SimulationNodeClass::nodeType_meshPrismaticLayer:
        color = MESH_CONTROL_COLOR;
        break;
    case SimulationNodeClass::nodeType_connectionGroup:
        color = CONNECTION_GROUP_COLOR;
        break;
    case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
    case SimulationNodeClass::nodeType_thermalAnalysisConvection:
    case SimulationNodeClass::nodeType_thermalAnalysisRadiation:
      case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
        color = THERMAL_BOUNDARY_COLOR;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
        color = COMPRESSION_ONLY_SUPPORT_COLOR;
        break;
    default:
        color = Quantity_NOC_YELLOW;
        break;
    }
    return color;
}

//! -------------------------
//! function: createColorBox
//! details:
//! -------------------------
void graphicsTools::createColorBox(double min, double max, int Nintervals, occHandle(AIS_ColorScaleExtended) &aColorBox)
{
    //! ---------------------------------------------------------
    //! https://tracker.dev.opencascade.org/view.php?id=25785
    //!
    //! GB: setColor sets the color of the text and of the frame
    //! the code of the class has been copied from OCC source
    //! ---------------------------------------------------------

    aColorBox = new AIS_ColorScaleExtended();
    aColorBox->SetRange(min,max);
    aColorBox->SetNumberOfIntervals(Nintervals);

    aColorBox->SetHeight(350);
    aColorBox->SetBreadth(75);  // it does not work
    Quantity_Color color(0,0,0,Quantity_TOC_RGB);
    aColorBox->SetColor(color);

    aColorBox->SetTextHeight(16);
    aColorBox->SetLabelAtBorder(Standard_True);

    occHandle(Graphic3d_TransformPers) trs = new Graphic3d_TransformPers(Graphic3d_TMF_2d);
    trs->SetCorner2d(Aspect_TOTP_LEFT_UPPER);

    Graphic3d_Vec2i offset(20,450);

    trs->SetOffset2d(offset);
    aColorBox->SetZLayer(Graphic3d_ZLayerId_TopOSD);
    aColorBox->SetTransformPersistence(trs);
    aColorBox->SetToUpdate();
}
