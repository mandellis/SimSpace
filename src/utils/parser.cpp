//! custom includes
#include "parser.h"
#include "tools.h"
#include "src/main/simulationmanager.h"

//! OCC
#include <Geom_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

parser::parser()
{

}

//! -------------------------------------------------------------
//! function: parseItem
//! details:
//! -------------------------------------------------------------
bool parser::parseItem(QExtendedStandardItem* item)
{
    if(item==NULL)
    {
        cerr<<"___error: NULL item____"<<endl;
        return false;
    }

    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    switch(node->getType())
    {
    case SimulationNodeClass::nodeType_namedSelectionGeometry:
        return parser::NamedSelection(item);
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
        return parser::FixedSupport(item);
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
        return parser::CylindricalSupport(item);
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
        return parser::FrictionlessSupport(item);
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
        return parser::Displacement(item);
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        return parser::Force(item);
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        return parser::Moment(item);
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        cout<<"parsing pressure"<<endl;
        return parser::Pressure(item);
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        return parser::ThermalCondition(item);
        break;
    case SimulationNodeClass::nodeType_meshBodyMeshMethod:
        return parser::MeshControl(item);
        break;
    case SimulationNodeClass::nodeType_meshBodyMeshControl:
        return parser::MeshControl(item);
        break;
    case SimulationNodeClass::nodeType_meshEdgeSize:
        return parser::MeshControl(item);
        break;
    case SimulationNodeClass::nodeType_meshVertexSize:
        return parser::MeshControl(item);
        break;
    case SimulationNodeClass::nodeType_connectionPair:
        return parser::ContactPair(item);
        break;
    case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
    case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
    case SimulationNodeClass::nodeType_solutionStructuralStress:
    case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
    case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
    case SimulationNodeClass::nodeType_solutionCFDK:
    case SimulationNodeClass::nodeType_solutionCFDpressure:
    case SimulationNodeClass::nodeType_solutionCFDvelocity:

        return parser::PostProcessingItem(item);
    break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
        return parser::BodyLoads(item);
        break;
    }
    return true;
}

//! -------------------------------------------------------------
//! function: NamedSelection
//! details:  parse a named selection
//! -------------------------------------------------------------
bool parser::NamedSelection(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::NamedSelection()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;

    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    itemToCheck = node->getPropertyItem("Tags");
    indexToCheck = node->getPropertyItem("Geometry")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0) invalidEntries.append(indexToCheck);

    parser::setItemBackground(node, invalidEntries);
    if(invalidEntries.length()==0)
    {
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
}

//! -------------------------------------------------------------
//! function: BodyLoads
//! details:  parse a body load
//! -------------------------------------------------------------
bool parser::BodyLoads(QExtendedStandardItem *item)
{
    {
        SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        cout<<"parser::BodyLoads()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;
        QExtendedStandardItem *itemToCheck;
        QModelIndex indexToCheck;
        QList<QModelIndex> invalidEntries;

        //! ----------------
        //! check the scope
        //! ----------------
        itemToCheck = node->getPropertyItem("Tags");
        if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
        else indexToCheck = node->getPropertyItem("Named selection")->index();

        std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
        if(vecLoc.size()==0)
        {
            //! the scope is empty => the item is not valid
            //! "Geometry" is an invalid entry
            invalidEntries.append(indexToCheck);
            parser::setItemBackground(node, invalidEntries);
        }
        else
        {
            for(std::vector<GeometryTag>::iterator itLoc = vecLoc.begin(); itLoc!=vecLoc.end(); ++itLoc)
            {
                GeometryTag loc = *itLoc;
                TopAbs_ShapeEnum type = loc.subShapeType;
                if(type != TopAbs_SOLID)
                {
                    invalidEntries.append(indexToCheck);
                    break;
                }
            }
        }
        parser::setItemBackground(node, invalidEntries);
        if(invalidEntries.length()==0)
        {
            item->setCheckState(Qt::Checked);
            node->setValid(true);
            return true;
        }
        else
        {
            item->setCheckState(Qt::Unchecked);
            node->setValid(false);
            return false;
        }
    }
}

//! ----------------------------------------------------------------------------------
//! function: FixedSupport
//! details:  parse a fixed support
//! ----------------------------------------------------------------------------------
bool parser::FixedSupport(QExtendedStandardItem *item)
{
    cout<<"parser::FixedSupport()->____parsing \"Fixed support\"____"<<endl;

    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToHighlight;
    QList<QModelIndex> invalidEntries;

    //! -------------------------------------------------------------------------
    //! check the scope: it can be defined through a direct "Geometry" selection
    //! or through a "Named selection"
    //! --------------------------------------------------------------------------
    itemToCheck = node->getPropertyItem("Tags");

    if(node->getPropertyItem("Geometry")!=NULL) indexToHighlight = node->getPropertyItem("Geometry")->index();
    else indexToHighlight = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        cerr<<"____the scope is empty____"<<endl;
        //! --------------------------------------------------------
        //! the scope is empty => the item is not valid
        //! "Geometry" is an invalid entry => use yellow background
        //! --------------------------------------------------------
        invalidEntries.append(indexToHighlight);
    }
    else
    {
        //! ----------------------------------------------------------
        //! if the "Tags" list is not empty, check what it contains:
        //! for "Fixed support" SOLIDS, COMPSOLIDS, and COMPOUNDS are
        //! note allowed
        //! ----------------------------------------------------------
        for(std::vector<GeometryTag>::iterator itLoc = vecLoc.begin(); itLoc!=vecLoc.end(); ++itLoc)
        {
            GeometryTag loc = *itLoc;
            TopAbs_ShapeEnum type = loc.subShapeType;
            if(type == TopAbs_SOLID || type == TopAbs_COMPOUND || type == TopAbs_COMPSOLID)
            {
                cerr<<"____the scope is not valid____"<<endl;
                invalidEntries.append(indexToHighlight);
                break;
            }
        }
        if(invalidEntries.length()==0)
        {
            cerr<<"____the scope is valid____"<<endl;
        }
    }

    parser::setItemBackground(node, invalidEntries);

    if(invalidEntries.length()==0)
    {
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
}

//! -------------------
//! function: Pressure
//! details:
//! -------------------
bool parser::Pressure(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::Pressure()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        //! the scope is empty => the item is not valid
        //! "Geometry" is an invalid entry
        invalidEntries.append(indexToCheck);
        parser::setItemBackground(node, invalidEntries);
    }
    else
    {
        for(std::vector<GeometryTag>::iterator itLoc = vecLoc.begin(); itLoc!=vecLoc.end(); ++itLoc)
        {
            GeometryTag loc = *itLoc;
            TopAbs_ShapeEnum type = loc.subShapeType;
            if(type != TopAbs_FACE)
            {
                invalidEntries.append(indexToCheck);
                break;
            }
        }
    }
    parser::setItemBackground(node, invalidEntries);
    //if(!item->isCheckable()) item->setCheckable(true);
    if(invalidEntries.length()==0)
    {
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
}

//! ----------------------------------------------------------------------------------
//! function: Force
//! details:  parse a force
//! ----------------------------------------------------------------------------------
bool parser::Force(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::Force()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        //! the scope is empty => the item is not valid
        //! "Geometry" is an invalid entry
        invalidEntries.append(indexToCheck);
        parser::setItemBackground(node, invalidEntries);
    }
    else
    {
        for(std::vector<GeometryTag>::iterator itLoc = vecLoc.begin(); itLoc!=vecLoc.end(); ++itLoc)
        {
            GeometryTag loc = *itLoc;
            TopAbs_ShapeEnum type = loc.subShapeType;
            if(type == TopAbs_SOLID || type == TopAbs_COMPOUND || type == TopAbs_COMPSOLID)
            {
                invalidEntries.append(indexToCheck);
                break;
            }
        }
    }
    parser::setItemBackground(node, invalidEntries);

    //! -----------------------------
    //! check the presence of values
    //! -----------------------------
    //! nothing to do

    if(invalidEntries.length()==0)
    {
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
}

//! ----------------------------------------------------------------------------------
//! function: ThermalCondition
//! details:  parse a thermal condition
//! ----------------------------------------------------------------------------------
bool parser::ThermalCondition(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::ThermalCondition()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        //! the scope is empty => the item is not valid
        //! "Geometry" is an invalid entry
        invalidEntries.append(indexToCheck);
        parser::setItemBackground(node, invalidEntries);
    }
    else
    {
        for(std::vector<GeometryTag>::iterator itLoc = vecLoc.begin(); itLoc!=vecLoc.end(); ++itLoc)
        {
            GeometryTag loc = *itLoc;
            TopAbs_ShapeEnum type = loc.subShapeType;
            if(type != TopAbs_SOLID && type != TopAbs_COMPSOLID)
            {
                invalidEntries.append(indexToCheck);
                break;
            }
        }
    }
    //if(!item->isCheckable()) item->setCheckable(true);
    if(invalidEntries.length()==0)
    {
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
}

//! ----------------------------------------------------------------------------------
//! function: CylindricalSupport
//! details:  parse a cylindrical support
//! ----------------------------------------------------------------------------------
bool parser::CylindricalSupport(QExtendedStandardItem *item)
{
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    geometryDataBase *gDB = sm->getDataBase();

    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    //cout<<"parser::CylindricalSupport()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;

    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        //! the scope is empty => the item is not valid
        invalidEntries.append(indexToCheck);
        parser::setItemBackground(node, invalidEntries);
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
    for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
    {
        GeometryTag curLoc = *it;

        //! check if the current topology is a face: if not, return false
        if(curLoc.subShapeType!=TopAbs_FACE)
        {
            //! one of the shapes within the scope is not a FACE
            //! the item is not valid
            invalidEntries.append(indexToCheck);
            parser::setItemBackground(node, invalidEntries);
            item->setCheckState(Qt::Unchecked);
            node->setValid(false);
            return false;
        }
        //! check if the current topology is a cylindrical face: if not return false
        int bodyIndex = curLoc.parentShapeNr;
        int faceNr = curLoc.subTopNr;
        const TopoDS_Face &curFace = TopoDS::Face(gDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceNr));
        occHandle(Geom_Surface) surface = BRep_Tool::Surface(curFace);
        GeomAdaptor_Surface adapt(surface);
        if(adapt.GetType()!=GeomAbs_Cylinder)
        {
            //cout<<"parser::CylindricalSupport()->____not a cylindrical surface____"<<endl;
            //! one of the shapes within the scope is not cylindrical
            invalidEntries.append(indexToCheck);
            parser::setItemBackground(node, invalidEntries);
            node->setValid(false);
            item->setCheckState(Qt::Unchecked);
            return false;
        }
    }

    if(invalidEntries.length()==0)
    {
        parser::setItemBackground(node, invalidEntries);
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        parser::setItemBackground(node, invalidEntries);
        item->setCheckState(Qt::Checked);
        node->setValid(false);
        return false;
    }
}

//! -------------------------------------------------------------
//! function: MeshControl
//! details:  parse a mesh control item
//! -------------------------------------------------------------
bool parser::MeshControl(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::MeshControl()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;

    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        invalidEntries.append(indexToCheck);
        parser::setItemBackground(node, invalidEntries);
        node->setValid(false);
        item->setCheckState(Qt::Unchecked);
        return false;
    }

    TopAbs_ShapeEnum shapeTypeReference;
    switch(node->getType())
    {
    case SimulationNodeClass::nodeType_meshVertexSize:
        cout<<"parser::MeshControl()->____parsing vertex mesh control____"<<endl;
        shapeTypeReference = TopAbs_VERTEX;
        break;
    case SimulationNodeClass::nodeType_meshEdgeSize:
        cout<<"parser::MeshControl()->____parsing edge mesh control____"<<endl;
        shapeTypeReference = TopAbs_EDGE;
        break;
    case SimulationNodeClass::nodeType_meshFaceSize:
        cout<<"parser::MeshControl()->____parsing face mesh control____"<<endl;
        shapeTypeReference = TopAbs_FACE;
        break;
    case SimulationNodeClass::nodeType_meshBodyMeshControl:
    case SimulationNodeClass::nodeType_meshBodyMeshMethod:
        cout<<"parser::MeshControl()->____parsing body mesh control____"<<endl;
        shapeTypeReference = TopAbs_SOLID;
        break;
    }

    for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
    {
        GeometryTag loc = *it;
        TopAbs_ShapeEnum shapeType = loc.subShapeType;
        if(shapeType != shapeTypeReference)
        {
            node->setValid(false);
            invalidEntries.append(indexToCheck);
            parser::setItemBackground(node, invalidEntries);
            item->setCheckState(Qt::Unchecked);
            return false;
        }
    }

    if(invalidEntries.length()==0)
    {
        parser::setItemBackground(node, invalidEntries);
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        parser::setItemBackground(node, invalidEntries);
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
}

//! -------------------------------------------------------------
//! function: ContactPair
//! details:  parse a contact pair
//! -------------------------------------------------------------
bool parser::ContactPair(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------------------------------
    //! check the presence of the "Tags master"
    //! ----------------------------------------
    std::vector<GeometryTag> vecLoc_master;
    itemToCheck = node->getPropertyItem("Tags master");
    indexToCheck = node->getPropertyItem("Master")->index();

    vecLoc_master= itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc_master.size()==0) invalidEntries.append(indexToCheck);
    parser::setItemBackground(node, invalidEntries);

    //! ----------------------------------------
    //! check the presence of the "Tags slave"
    //! ----------------------------------------
    std::vector<GeometryTag> vecLoc_slave;
    itemToCheck = node->getPropertyItem("Tags slave");
    indexToCheck = node->getPropertyItem("Slave")->index();

    vecLoc_slave = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc_slave.size()==0) invalidEntries.append(indexToCheck);
    parser::setItemBackground(node, invalidEntries);
    if(invalidEntries.length()==0)
    {
        node->setValid(true);
        item->setCheckState(Qt::Checked);
        return true;
    }
    else
    {
        node->setValid(false);
        item->setCheckState(Qt::Unchecked);
        return false;
    }
}

//! -------------------------------------------------------------
//! function: Displacement
//! details:  parse a Displacement boundary condition
//! -------------------------------------------------------------
bool parser::Displacement(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::MeshControl()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------------------------------------------------------------------------------------
    //! check the scope: it should not be empty, and it should contain only vertexes, edges, or faces
    //! ----------------------------------------------------------------------------------------------
    //itemToCheck = node->getPropertyItem("Geometry");
    //indexToCheck = itemToCheck->index();

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        //! the scope is empty => the item is not valid
        invalidEntries.append(indexToCheck);
    }
    else
    {
        for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
        {
            GeometryTag loc = *it;
            TopAbs_ShapeEnum type = loc.subShapeType;
            if(type == TopAbs_SOLID || type == TopAbs_COMPOUND || type == TopAbs_COMPSOLID)
            {
                invalidEntries.append(indexToCheck);
                break;
            }
        }
        parser::setItemBackground(node, invalidEntries);
    }

    //! ------------------------------------------------------------------
    //! check the components: at least one component should not be "free"
    //! ------------------------------------------------------------------
    QExtendedStandardItem *item_componentX = node->getPropertyItem("X component");
    QExtendedStandardItem *item_componentY = node->getPropertyItem("Y component");
    QExtendedStandardItem *item_componentZ = node->getPropertyItem("Z component");

    Property::loadDefinition load_componentX = item_componentX->data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
    Property::loadDefinition load_componentY = item_componentY->data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
    Property::loadDefinition load_componentZ = item_componentZ->data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();

    if(load_componentX==Property::loadDefinition_free && load_componentY==Property::loadDefinition_free && load_componentZ==Property::loadDefinition_free)
    {
        QModelIndex index_Xcomponent = item_componentX->index();
        QModelIndex index_Ycomponent = item_componentY->index();
        QModelIndex index_Zcomponent = item_componentZ->index();
        invalidEntries.append(index_Xcomponent);
        invalidEntries.append(index_Ycomponent);
        invalidEntries.append(index_Zcomponent);
    }

    if(!invalidEntries.isEmpty())
    {
        //cout<<"Number of invalid entries: "<<invalidEntries.size()<<endl;
        parser::setItemBackground(node, invalidEntries);
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
    else
    {
        parser::setItemBackground(node, invalidEntries);
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
}

//! ----------------------------------------------------------------------------------
//! function: FrictionlessSupport
//! details:  parse a frictionless support
//! ----------------------------------------------------------------------------------
bool parser::FrictionlessSupport(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::FrictionlessSupport()->____parsing: \"Frictionless support\"____"<<endl;
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        //! the scope is empty => the item is not valid
        //! "Geometry" is an invalid entry
        invalidEntries.append(indexToCheck);
        parser::setItemBackground(node, invalidEntries);
    }
    else
    {
        for(std::vector<GeometryTag>::iterator itLoc = vecLoc.begin(); itLoc!=vecLoc.end(); ++itLoc)
        {
            GeometryTag loc = *itLoc;
            TopAbs_ShapeEnum type = loc.subShapeType;
            if(type != TopAbs_FACE)
            {
                invalidEntries.append(indexToCheck);
                break;
            }
        }
    }
    parser::setItemBackground(node, invalidEntries);
    if(invalidEntries.length()==0)
    {
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
}

//! ----------------------------------------------------------------------------------
//! function: postProcessingItem
//! details:  parse a post processing item
//! ----------------------------------------------------------------------------------
bool parser::PostProcessingItem(QExtendedStandardItem *item)
{
    /*
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::FrictionlessSupport()->____parsing \"post processing item\"____"<<endl;
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        //! the scope is empty => the item is not valid
        //! "Geometry" is an invalid entry
        invalidEntries.append(indexToCheck);
        parser::setItemBackground(node, invalidEntries);
    }
    parser::setItemBackground(node, invalidEntries);
    if(invalidEntries.length()==0)
    {
        //! ----------------------------------------------------------------------
        //! check if the item has a post object and the post object contains data
        //! ----------------------------------------------------------------------
        QStandardItem* itemPostObject = node->getPropertyItem("Post object");
        if(itemPostObject==NULL)
        {
            //! no error in item definition, but the data lack
            item->setCheckState(Qt::PartiallyChecked);
        }
        else
        {
            postObject aPostObject = itemPostObject->data(Qt::UserRole).value<Property>().getData().value<postObject>();
            if(aPostObject.getData().isEmpty()) item->setCheckState(Qt::PartiallyChecked);
            else item->setCheckState(Qt::Checked);
        }
        node->setValid(true);
        return true;
    }
    else
    {
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
    */
    return true;
}

//! -------------------------------------------------------------
//! function: setItemBackground
//! details:
//! -------------------------------------------------------------
void parser::setItemBackground(SimulationNodeClass *node, QList<QModelIndex> invalidEntries)
{
    if(invalidEntries.length()>0)
    {
        for(int i=0; i<invalidEntries.length();i++)
        {
            QStandardItem *anInvalidItem = node->getModel()->itemFromIndex(invalidEntries.at(i));
            //! use solid yellow for an invalid item
            QBrush brush(Qt::yellow, Qt::SolidPattern);
            anInvalidItem->setBackground(brush);
        }
    }
    else
    {
        QStandardItemModel *nodeModel = node->getModel();
        QBrush brush;
        for(int i=0; i<nodeModel->rowCount(); i++)
        {
            QModelIndex separator_index =nodeModel->index(i,0);
            QStandardItem *separator = nodeModel->itemFromIndex(separator_index);
            int NbRowSeparator = separator->rowCount();
            //! this removes the brush
            for(int j=0; j<NbRowSeparator; j++) separator->child(j,1)->setBackground(brush);
        }
    }
}



//! ----------------------------------------------------------------------------------
//! function: Moment
//! details:  parse a Moment
//! ----------------------------------------------------------------------------------
bool parser::Moment(QExtendedStandardItem *item)
{
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    cout<<"parser::Force()->____parsing: "<<node->getName().toStdString()<<"____"<<endl;
    QExtendedStandardItem *itemToCheck;
    QModelIndex indexToCheck;
    QList<QModelIndex> invalidEntries;

    //! ----------------
    //! check the scope
    //! ----------------
    itemToCheck = node->getPropertyItem("Tags");
    if(node->getPropertyItem("Geometry")!=NULL) indexToCheck = node->getPropertyItem("Geometry")->index();
    else indexToCheck = node->getPropertyItem("Named selection")->index();

    std::vector<GeometryTag> vecLoc = itemToCheck->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
    if(vecLoc.size()==0)
    {
        //! the scope is empty => the item is not valid
        //! "Geometry" is an invalid entry
        invalidEntries.append(indexToCheck);
        parser::setItemBackground(node, invalidEntries);
    }
    else
    {
        for(std::vector<GeometryTag>::iterator itLoc = vecLoc.begin(); itLoc!=vecLoc.end(); ++itLoc)
        {
            GeometryTag loc = *itLoc;
            TopAbs_ShapeEnum type = loc.subShapeType;
            if(type == TopAbs_SOLID || type == TopAbs_COMPOUND || type == TopAbs_COMPSOLID)
            {
                invalidEntries.append(indexToCheck);
                break;
            }
        }
    }
    parser::setItemBackground(node, invalidEntries);

    //! -----------------------------
    //! check the presence of values
    //! -----------------------------
    //! nothing to do

    if(invalidEntries.length()==0)
    {
        item->setCheckState(Qt::Checked);
        node->setValid(true);
        return true;
    }
    else
    {
        item->setCheckState(Qt::Unchecked);
        node->setValid(false);
        return false;
    }
}
