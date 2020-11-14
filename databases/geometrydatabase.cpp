//! ----------------
//! custom includes
//! ----------------
#include <geometrydatabase.h>
#include "simulationnodeclass.h"
#include "property.h"
#include "qbackgroundevent.h"
#include "tools.h"
#include "ccout.h"
#include "topologytools.h"
#include "topods_shape_reg.h"
#include "src/geometry/geometryhealing.h"
#include "OCCface.h"

//! ----
//! OCC
//! ----
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_CompSolid.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Builder.hxx>
#include <TColStd_ListIteratorOfListOfAsciiString.hxx>
#include <TColStd_Array1OfAsciiString.hxx>
#include <GProp_GProps.hxx>
#include <GProp_PrincipalProps.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <ShapeAnalysis_ShapeContents.hxx>

//! ---
//! Qt
//! ---
#include <QStandardItem>
#include <QDebug>
#include <QApplication>

//! --------------------------------------
//! function: compute geometry properties
//! details:
//! --------------------------------------
void computeGeometryProperties(const TopoDS_Shape &aShape, geometryProperties_ &gprops)
{
    if(!aShape.IsNull())
    {
        GProp_GProps prop;

        switch (aShape.ShapeType())
        {
        case TopAbs_SOLID:
        case TopAbs_COMPOUND:
            BRepGProp::VolumeProperties(aShape,prop);
            gprops.volume = prop.Mass();
            break;

        case TopAbs_SHELL:
        case TopAbs_FACE:
            BRepGProp::SurfaceProperties(aShape,prop);
            gprops.area = prop.Mass();
            // a 2D body has not a volume
            gprops.volume = 0.0;
            break;
        }

        //! center of mass
        gp_Pnt CM = prop.CentreOfMass();
        gprops.xB = CM.X();
        gprops.yB = CM.Y();
        gprops.zB = CM.Z();

        //! principal moments of inertia
        Standard_Real Ixx, Iyy, Izz;
        prop.PrincipalProperties().Moments(Ixx,Iyy,Izz);
        gprops.Ixx = Ixx;
        gprops.Iyy = Iyy;
        gprops.Izz = Izz;

        //! volume/area ...
    }
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
geometryDataBase::~geometryDataBase()
{
    cout<<"geometryDataBase::~geometryDataBase()->____DESTRUCTOR CALLED____"<<endl;
}

/*
//! ------------------------------
//! function: default constructor
//! details:
//! ------------------------------
geometryDataBase::geometryDataBase(QObject *parent):
    QObject(parent),
    fullSourceFileName(""),
    myN3D(0),
    myN2D(0),
    myN1D(0),
    myTopoDS_Shape(TopoDS_Shape())
{
    ccout("geometryDataBase::geometryDataBase()->____DEFAULT CONSTRUCTOR CALLED____");
    cerr<<"geometryDataBase::geometryDataBase()->____default constructor called____"<<endl;

    //! ------------------------
    //! the standard item model
    //! ------------------------
    myModel = new QStandardItemModel(this);

    //! ----------------------------------------
    //! create the QStandardItemModel root item
    //! ----------------------------------------
    myInvisibleRootItem = static_cast<QExtendedStandardItem*>(myModel->invisibleRootItem());

    //! -----------------------------
    //! create the "Model" root item
    //! -----------------------------
    this->createModelRootItem();

    //! -------------------------
    //! create the "Import" item
    //! -------------------------
    this->createImportItem();

    //! --------------------------------
    //! create the "Geometry" root item
    //! --------------------------------
    this->createGeometryRoot();

    //! -------------------------------------
    //! create the "Coordinate systems" item
    //! -------------------------------------
    this->createCoordinateSystemsRoot();
}
*/
//! ------------------------------
//! function: constructor II
//! details:  with initialization
//! ------------------------------
geometryDataBase::geometryDataBase(const TopoDS_Shape &theShape, const QString &theFilePath, QObject *parent):
    QObject(parent),
    myTopoDS_Shape(theShape),
    fullSourceFileName(theFilePath)
{
    cout<<"geometryDataBase::geometryDataBase()->____constructor II called____"<<endl;

    //! ------------------------
    //! the standard item model
    //! ------------------------
    myModel = new QStandardItemModel(this);

    //! ----------------------------------------
    //! create the QStandardItemModel root item
    //! ----------------------------------------
    myInvisibleRootItem = static_cast<QExtendedStandardItem*>(myModel->invisibleRootItem());

    //! -----------------------------
    //! create the "Model" root item
    //! -----------------------------
    this->createModelRootItem();

    //! -------------------------
    //! create the "Import" item
    //! -------------------------
    this->createImportItem();

    //! ------------------------------
    //! create the Geometry root item
    //! ------------------------------
    this->createGeometryRoot();

    //! --------------------------
    //! create the geometry nodes
    //! --------------------------
    if(theShape.IsNull())
    {
        cerr<<"geometryDataBase::geometryDataBase()->____cannot create the geometry items: the input shape is null____"<<endl;
    }
    else this->createGeometryNodes();

    //! -------------------------------------
    //! create the "Coordinate systems" item
    //! -------------------------------------
    this->createCoordinateSystemsRoot();

    //! ---------------------------------------------------
    //! create the topology maps and initialize the arrays
    //! ---------------------------------------------------
    if(theShape.IsNull())
    {
        cerr<<"geometryDataBase::geometryDataBase()->____cannot init shapes and topology maps: the input shape is null____"<<endl;
    }
    else this->buildMaps(theShape);

    //! ----------------
    //! print a summary
    //! ----------------
    this->printSummary();
}

//! ------------------------------
//! function: createStandardModel
//! details:
//! ------------------------------
void geometryDataBase::createStandardModel()
{
    cout<<"geometryDataBase::createStandardModel()->____function called____"<<endl;

    //! -----------------------------
    //! create the "Model" root item
    //! -----------------------------
    this->createModelRootItem();

    //! -------------------------
    //! create the "Import" item
    //! -------------------------
    this->createImportItem();

    //! ------------------------------
    //! create the Geometry root item
    //! ------------------------------
    this->createGeometryRoot();

    //! -------------------------------------
    //! create the "Coordinate systems" item
    //! -------------------------------------
    this->createCoordinateSystemsRoot();
}

//! -----------------------
//! function: printSummary
//! details:
//! -----------------------
void geometryDataBase::printSummary()
{
    cout<<"/-----------------------------------/"<<endl;
    cout<<"/- SUMMARY OF IMPORT PROCESS        /"<<endl;
    cout<<"/-----------------------------------/"<<endl;
    int NbCompositeSolids = 0;
    int NbSolids = 0;
    int NbShells = 0;
    int NbFaces = 0;
    int NbWires = 0;
    int NbEdges = 0;
    for(QMap<int,TopoDS_Shape>::iterator it=bodyMap.begin(); it!=bodyMap.end(); ++it)
    {
        const TopoDS_Shape &shape = it.value();
        if(shape.ShapeType()==TopAbs_COMPSOLID) NbCompositeSolids++;
        if(shape.ShapeType()==TopAbs_SOLID) NbSolids++;
        if(shape.ShapeType()==TopAbs_SHELL) NbShells++;
        if(shape.ShapeType()==TopAbs_FACE) NbFaces++;
        if(shape.ShapeType()==TopAbs_WIRE) NbWires++;
        if(shape.ShapeType()==TopAbs_EDGE) NbEdges++;
    }
    cout<<"____composite solids: "<<NbCompositeSolids<<"____"<<endl;
    cout<<"____free solids: "<<NbSolids<<"____"<<endl;
    cout<<"____free shells: "<<NbShells<<"____"<<endl;
    cout<<"____free faces: "<<NbFaces<<"____"<<endl;
    cout<<"____free wires: "<<NbWires<<"____"<<endl;
    cout<<"____free edges: "<<NbEdges<<"____"<<endl;
}

//!-------------------------
//! function: getSubShapeNr
//! details:
//! -------------------------
bool geometryDataBase::getSubShapeNr(const TopoDS_Shape &aSubShape, int &mainShapeIndex, int &subShapeIndex, TopAbs_ShapeEnum &subShapeType) const
{
    //cout<<"geometryDataBase::getSubShapeNr()->____function called____"<<endl;
    TopAbs_ShapeEnum type = aSubShape.ShapeType();
    bool isTheShapeContained = false;

    switch(type)
    {
    case TopAbs_COMPSOLID:
    {
        ;
    }
        break;

    case TopAbs_SOLID:
    {
        mainShapeIndex = bodyMap.key(aSubShape);
        subShapeIndex = mainShapeIndex;
        subShapeType = TopAbs_SOLID;
        isTheShapeContained = true;
    }
        break;

    case TopAbs_FACE:
    case TopAbs_EDGE:
    case TopAbs_VERTEX:
    {
        for(QMap<int,TopoDS_Shape>::const_iterator it = bodyMap.cbegin(); it!=bodyMap.cend(); ++it)
        {
            int k = it.key();
            switch(type)
            {
            case TopAbs_FACE: isTheShapeContained = MapOfBodyTopologyMap.value(k).faceMap.Contains(aSubShape); break;
            case TopAbs_EDGE: isTheShapeContained = MapOfBodyTopologyMap.value(k).edgeMap.Contains(aSubShape); break;
            case TopAbs_VERTEX: isTheShapeContained = MapOfBodyTopologyMap.value(k).vertexMap.Contains(aSubShape); break;
            }

            if(isTheShapeContained==true)
            {
                subShapeType = aSubShape.ShapeType();
                mainShapeIndex = k;
                switch(type)
                {
                case TopAbs_FACE: subShapeIndex = MapOfBodyTopologyMap.value(k).faceMap.FindIndex(aSubShape); break;
                case TopAbs_EDGE: subShapeIndex = MapOfBodyTopologyMap.value(k).edgeMap.FindIndex(aSubShape); break;
                case TopAbs_VERTEX: subShapeIndex = MapOfBodyTopologyMap.value(k).vertexMap.FindIndex(aSubShape); break;
                }
            }
        }
    }
        break;
    }
    return isTheShapeContained;
}

//! --------------------------
//! function: getTopologyMaps
//! details:
//! --------------------------
TopologyMap geometryDataBase::getTopologyMaps(const TopoDS_Shape &aShape)
{
    TopologyMap aTopologyMap;
    TopExp::MapShapes(aShape,TopAbs_SHELL,aTopologyMap.shellMap);
    TopExp::MapShapes(aShape,TopAbs_FACE,aTopologyMap.faceMap);
    TopExp::MapShapes(aShape,TopAbs_WIRE,aTopologyMap.wireMap);
    TopExp::MapShapes(aShape,TopAbs_EDGE,aTopologyMap.edgeMap);
    TopExp::MapShapes(aShape,TopAbs_VERTEX,aTopologyMap.vertexMap);
    return aTopologyMap;
}

//!-------------------------------
//! function: createGeometryNodes
//! details:
//! ------------------------------
void geometryDataBase::createGeometryNodes()
{
    cout<<"geometryDataBase::createGeometryNodes()->____creating the geometry items____"<<endl;

    //! -------------
    //! generic data
    //! -------------
    QVariant data;

    //! -----------------------------------------------
    //! create the other nodes (1D, 2D, 3D, 4D bodies)
    //! -----------------------------------------------
    for(QMap<int,TopoDS_Shape>::iterator bodyIt = bodyMap.begin(); bodyIt!=bodyMap.end(); ++bodyIt)
    {
        TopoDS_Shape bodyShape = bodyIt.value();
        int bodyIndex = bodyIt.key();

        //! ---------------------
        //! vector of properties
        //! ---------------------
        QVector<Property> props;

        //! -------------------------------------------------------
        //! set the properties
        //! warning: at this stage MapOfBodyNames.value(bodyIndex)
        //!          will return an empty string
        //! -------------------------------------------------------
        QString name = QString("%1").arg(MapOfBodyNames.value(bodyIndex));
        data.setValue(name);
        Property prop_name("Name",data,Property::PropertyGroup_Definition);
        props.push_back(prop_name);

        Property::elementControl theElementControl = Property::elementControl_programControlled;    //! default value
        Property::integrationScheme theIntegrationScheme;
        if(theElementControl==Property::elementControl_programControlled)
        {
            ;   //! don't add
        }
        else
        {
            theIntegrationScheme = Property::integrationScheme_full;
            data.setValue(theIntegrationScheme);
            Property property_integrationScheme("Integration scheme",data,Property::PropertyGroup_Definition);
            props.push_back(property_integrationScheme);
        }

        //! ------
        //! scope
        //! ------
        int mainShapeIndex, subShapeIndex;
        TopAbs_ShapeEnum subShapeType;
        this->getSubShapeNr(bodyShape,mainShapeIndex,subShapeIndex,subShapeType);
        GeometryTag loc;
        loc.isParent = true;
        loc.parentShapeNr=mainShapeIndex;
        loc.subTopNr=subShapeIndex;
        loc.subShapeType=subShapeType;
        std::vector<GeometryTag> vecLoc;
        vecLoc.push_back(loc);
        data.setValue(vecLoc);
        Property property_geometry("Geometry",data,Property::PropertyGroup_Scope);
        props.push_back(property_geometry);
        Property prop_mapIndex("Map index",QVariant(bodyIndex),Property::PropertyGroup_Scope);
        props.push_back(prop_mapIndex);

        //! -------------------------------
        //! record the shape into the item
        //! -------------------------------
        const TopoDS_Shape &curShape = bodyIt.value();
        data.setValue(curShape);
        Property prop_shape("Shape",data,Property::PropertyGroup_Scope);
        props.push_back(prop_shape);

        //! -----------------------
        //! store the topology map
        //! -----------------------
        //TopologyMap aTopologyMap= this->getTopologyMaps(curShape);
        //data.setValue(aTopologyMap);
        //Property prop_topologyMap("Topology map",data,Property::PropertyGroup_Scope);
        //props.push_back(prop_topologyMap);

        //! ---------------------
        //! geometric properties
        //! ---------------------
        geometryProperties_ myProp;
        computeGeometryProperties(bodyShape,myProp);
        Property prop_xB("Center of mass x",myProp.xB,Property::PropertyGroup_Properties);
        Property prop_yB("Center of mass y",myProp.yB,Property::PropertyGroup_Properties);
        Property prop_zB("Center of mass z",myProp.zB,Property::PropertyGroup_Properties);
        Property prop_Ixx("Moment of inertia Ixx",myProp.Ixx,Property::PropertyGroup_Properties);
        Property prop_Iyy("Moment of inertia Iyy",myProp.Iyy,Property::PropertyGroup_Properties);
        Property prop_Izz("Moment of inertia Izz",myProp.Izz,Property::PropertyGroup_Properties);

        props.push_back(prop_xB);
        props.push_back(prop_yB);
        props.push_back(prop_zB);
        props.push_back(prop_Ixx);
        props.push_back(prop_Iyy);
        props.push_back(prop_Izz);

        //! -------------------
        //! graphic properties
        //! -------------------
        bool isVisible = true;
        Property property_visible("Visible",isVisible,Property::PropertyGroup_GraphicProperties);
        props.push_back(property_visible);

        //! ------------------
        //! bounding box size
        //! ------------------
        double LX, LY, LZ;
        getBoundingBox(bodyShape,LX,LY,LZ);
        Property BB_sizeX("Length X",LX,Property::PropertyGroup_BoundingBox);
        Property BB_sizeY("Length Y",LY,Property::PropertyGroup_BoundingBox);
        Property BB_sizeZ("Length Z",LZ,Property::PropertyGroup_BoundingBox);
        props.push_back(BB_sizeX);
        props.push_back(BB_sizeY);
        props.push_back(BB_sizeZ);

        //! ------------------
        //! property material
        //! ------------------
        int value = 0;
        data.setValue(value);
        Property prop_material("Assignment",data,Property::PropertyGroup_Material);
        props.push_back(prop_material);

        //! ---------------------
        //! property suppression
        //! ---------------------
        Property::SuppressionStatus ss = Property::SuppressionStatus_Active;
        data.setValue(ss);
        Property prop_suppression("Suppressed",data,Property::PropertyGroup_Definition);
        props.push_back(prop_suppression);

        //! -----
        //! node
        //! -----
        SimulationNodeClass *nodeGeometry = new SimulationNodeClass(name,SimulationNodeClass::nodeType_geometryBody,props,this);

        //! ---------
        //! time tag
        //! ---------
        nodeGeometry->addTimeTag();

        QExtendedStandardItem *GeometryItem = new QExtendedStandardItem();
        data.setValue(nodeGeometry);
        GeometryItem->setData(name,Qt::DisplayRole);
        GeometryItem->setData(data,Qt::UserRole);
        GeometryItem->setEditable(false);
        GeometryItem->setCheckState(Qt::Checked);

        //! ------------------------------------------
        //! append the item to the geometry root item
        //! ------------------------------------------
        GeometryRootItem->appendRow(GeometryItem);
        GeometryRootItem->setEditable(false);

        //! ------------------------
        //! add the parent time tag
        //! ------------------------
        QString parentTimeTag = GeometryRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
        data.setValue(parentTimeTag);
        Property property_parentTimeTag("Parent time tag",data,Property::PropertyGroup_Identifier);
        nodeGeometry->addProperty(property_parentTimeTag);
    }
}

//!-----------------------------------------------------------------------------
//! function: transferNodesNames
//! details:  transfer the names from the array of names to the items and nodes
//!-----------------------------------------------------------------------------
void geometryDataBase::transferNames()
{
    //!cout<<"geometryDataBase::transferNames()->____function called____"<<endl;
    if(GeometryRootItem->hasChildren())
    {
        //! ------------------------------------------------------------------
        //! change the name of the node, using the name of the CAD components
        //! ------------------------------------------------------------------
        int N = GeometryRootItem->rowCount();
        for(int row=0; row<N; row++)
        {
            QStandardItem *item = GeometryRootItem->child(row,0);
            SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(node->getType()==SimulationNodeClass::nodeType_pointMass) continue;
            int bodyIndex = node->getPropertyValue<int>("Map index");

            QString updatedName = QString("%1").arg(MapOfBodyNames.value(bodyIndex));

            //! -----------------------------------------------------------
            //! change the label of the item in the tree (Qt::DisplayRole)
            //! -----------------------------------------------------------
            item->setData(updatedName,Qt::DisplayRole);

            //! ----------------------------
            //! change the name of the node
            //! ----------------------------
            QVariant data;
            data.setValue(updatedName);
            //SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
            node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
            node->setName(updatedName);

            Property prop_name("Name",data,Property::PropertyGroup_Definition);
            node->getModel()->blockSignals(true);
            node->replaceProperty("Name",prop_name);
            node->getModel()->blockSignals(false);
        }
    }
}

//! -------------------------
//! function: getBoundingBox
//! details:
//! -------------------------
void geometryDataBase::getBoundingBox(const TopoDS_Shape &shape, double &L1, double &L2, double &L3)
{
    Bnd_Box boundingBox;
    BRepBndLib::Add(shape, boundingBox);
    Standard_Real Xmin,Ymin,Zmin,Xmax,Ymax,Zmax;
    boundingBox.Get(Xmin,Ymin,Zmin,Xmax,Ymax,Zmax);
    L1 = fabs(Xmax-Xmin);
    L2 = fabs(Ymax-Ymin);
    L3 = fabs(Zmax-Zmin);
}


//!-------------------------------------------
//! function: create "Coordinate system root"
//! details:
//!-------------------------------------------
#include "handle_ais_trihedron_reg.h"
#include "markers.h"
void geometryDataBase::createCoordinateSystemsRoot()
{
    cout<<"geometryDataBase::createCoordinateSystemsRoot()->____function called____"<<endl;

    //! -----------------------
    //! Coordinate system root
    //! -----------------------
    QVariant data;
    QVector<Property> props;

    SimulationNodeClass *node = new SimulationNodeClass("Coordinate systems", SimulationNodeClass::nodeType_coordinateSystems, props, this);
    data.setValue(node);

    //! ---------------------------
    //! create and append the item
    //! ---------------------------
    CoordinateSystemsRootitem = new QExtendedStandardItem();
    CoordinateSystemsRootitem->setData(data,Qt::UserRole);
    CoordinateSystemsRootitem->setData("Coordinate systems",Qt::DisplayRole);
    CoordinateSystemsRootitem->setEditable(false);
    myRootItem->appendRow(CoordinateSystemsRootitem);

    //! -----------------------------------------
    //! add the time tag and the parent time tag
    //! -----------------------------------------
    node->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    node->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! -------------------------------------------
    //! create the "Global coordinate system" item
    //! -------------------------------------------
    Property prop_originX("Origin X",0,Property::PropertyGroup_Origin);
    Property prop_originY("Origin Y",0,Property::PropertyGroup_Origin);
    Property prop_originZ("Origin Z",0,Property::PropertyGroup_Origin);

    props.push_back(prop_originX);
    props.push_back(prop_originY);
    props.push_back(prop_originZ);

    //! -----------------
    //! directional data
    //! -----------------
    QVariant datax,datay,dataz;
    QVector<double> dirX, dirY, dirZ;
    dirX.append(1); dirX.append(0); dirX.append(0);
    dirY.append(0); dirY.append(1); dirY.append(0);
    dirZ.append(0); dirZ.append(0); dirZ.append(1);
    QVector<QVector<double>> directionalData;
    directionalData<<dirX<<dirY<<dirZ;

    datax.setValue(dirX);
    datay.setValue(dirY);
    dataz.setValue(dirZ);

    Property property_X_axisData("X axis data",datax,Property::PropertyGroup_DirectionalData);
    Property property_Y_axisData("Y axis data",datay,Property::PropertyGroup_DirectionalData);
    Property property_Z_axisData("Z axis data",dataz,Property::PropertyGroup_DirectionalData);
    props.push_back(property_X_axisData);
    props.push_back(property_Y_axisData);
    props.push_back(property_Z_axisData);

    //! --------------------------------------
    //! add the trihedron as a graphic object
    //! --------------------------------------
    QVector<double> origin;
    origin.push_back(0.0);
    origin.push_back(0.0);
    origin.push_back(0.0);

    AIS_Trihedron_handle_reg theTrihedron = markers::buildTrihedron(origin,directionalData,20);
    data.setValue(theTrihedron);
    Property prop_trihedron("Graphic object",data,Property::PropertyGroup_GraphicObjects);
    props.push_back(prop_trihedron);

    SimulationNodeClass *node_globalCS = new SimulationNodeClass("Global coordinate system", SimulationNodeClass::nodeType_coordinateSystem_global, props, this);

    //! -----------------------------------------
    //! create the global coordinate system item
    //! -----------------------------------------
    QExtendedStandardItem *item_globalCS = new QExtendedStandardItem();
    data.setValue(node_globalCS);
    item_globalCS->setData(data, Qt::UserRole);
    item_globalCS->setData("Global coordinate system", Qt::DisplayRole);
    item_globalCS->setEditable(false);

    //! -----------------------------------------
    //! add the time tag and the parent time tag
    //! -----------------------------------------
    node_globalCS->addTimeTag();
    timeTag = node->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    node_globalCS->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! -----------------------------------------
    //! append the global coordinate system item
    //! -----------------------------------------
    CoordinateSystemsRootitem->appendRow(item_globalCS);
}

//!-------------------------------------------------------
//! function: toleranceForHealing
//! details:  return tolerance value for geometry healing
//! ------------------------------------------------------
double geometryDataBase::toleranceForHealing(const TopoDS_Shape &shape, int typeOfTolerance)
{
    double toleranceValue;
    double level;
    switch(typeOfTolerance)
    {
    case 1: level = -4.0; toleranceValue = boundingBox(shape)*pow(10.0,level); break; //! "Normal"
    case 2: level = -3.0; toleranceValue = boundingBox(shape)*pow(10.0,level); break; //! "Loose"
    default: level = -4.0; toleranceValue = boundingBox(shape)*pow(10.0,level); break; //! "Default" - use "Normal"
    }
    return toleranceValue;
}

//!-----------------------------------------------
//! function: buildMaps
//! details:  build the maps and setup the arrays
//! ----------------------------------------------
#include "geomtoolsclass.h"
void geometryDataBase::buildMaps(const TopoDS_Shape &shape)
{
    ccout("geometryDataBase::buildMaps()->____function called: building maps____");

    TopoDS_Compound aComp = TopoDS::Compound(shape);
    QList<TopoDS_Shape> csolids,solids,shells,faces,wires,edges;
    const TopoDS_Shape &compound = TopologyTools::exploreCompound(aComp,csolids,solids,shells,faces,wires,edges,true,true,true);
    Q_UNUSED(compound)

    //! -----------------------
    //! Total number of bodies
    //! -----------------------
    myN3D = solids.length();
    myN2D = shells.length()+faces.length();
    myN1D = wires.length()+edges.length();
    int Nt = myN3D+myN2D+myN1D;

    cout<<"\n@____Summary of loaded shape types____@"<<endl;
    cout<<"____Nb solids: "<<myN3D<<"____"<<endl;
    cout<<"____Nb surface bodies: "<<myN2D<<"____"<<endl;
    cout<<"____Nb line bodies: "<<myN1D<<"____"<<endl;
    cout<<"____Total number of bodies: "<<Nt<<"____"<<endl<<endl;

    //! -----------------------------
    //! access the importing options
    //! and the healing parameters
    //! -----------------------------
    bool importSolidBodies = true;
    bool importSurfaceBodies = true;
    bool importLineBodies = true;
    bool healGeometry = false;

    SimulationNodeClass *importNode = myInvisibleRootItem->child(0,0)->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();

    cout<<"____check if the \"Import\" node is present____"<<endl;
    if(importNode->getType()==SimulationNodeClass::nodeType_import) cout<<"____\"Import\" node not present____\n"<<endl;
    else cout<<"____\"Import\" node present____\n"<<endl;

    if(importNode->getType()==SimulationNodeClass::nodeType_import)
    {
        importSolidBodies = importNode->getPropertyValue<bool>("Solid bodies");
        importSurfaceBodies = importNode->getPropertyValue<bool>("Surface bodies");
        importLineBodies = importNode->getPropertyValue<bool>("Line bodies");
        healGeometry = importNode->getPropertyValue<bool>("Geometry healing");
    }
    else
    {
        cout<<"geometryDataBase::buildMaps()->____\"Import\" item not found: use default options____"<<endl;
    }

    int bodyIndex = 0;
    //! --------------------------------
    //! solids (solid parts): 3D bodies
    //! --------------------------------
    if(importSolidBodies)
    {
        for(int i=0; i<solids.length(); i++)
        {
            cout<<"____building the topology map for a \"Solid\"____"<<endl;
            bodyIndex++;
            const TopoDS_Shape &curBody = solids.at(i);

            if(healGeometry)
            {
                //! -----------------------------
                //! retrieve the tolerance value
                //! -----------------------------
                double toleranceValue;
                int toleranceType = importNode->getPropertyValue<double>("Tolerance");
                if(toleranceType!=0)
                {
                    //! auto-computed tolerance value
                    toleranceValue = this->toleranceForHealing(curBody,toleranceType);
                }
                else
                {
                    //! user defined tolerance value
                    toleranceValue = importNode->getPropertyValue<double>("Tolerance value");
                }

                cout<<"____geometry healing is active: tolerance for body nr: "<<bodyIndex<<" tolerance = "<<toleranceValue<<"____"<<endl;

                //! ---------------------
                //! geometry defeaturing
                //! ---------------------
                TopoDS_Shape curBody = solids.at(i);
                geometryHealing aGeometryHealing(curBody);

                aGeometryHealing.perform(toleranceValue);
                curBody = aGeometryHealing.getShape();
            }

            bodyMap.insert(bodyIndex,curBody);
            MapOfIsActive.insert(bodyIndex,true);
            TopologyMap map;
            TopExp::MapShapes(curBody,TopAbs_SHELL,map.shellMap);
            TopExp::MapShapes(curBody,TopAbs_FACE,map.faceMap);
            TopExp::MapShapes(curBody,TopAbs_WIRE,map.wireMap);
            TopExp::MapShapes(curBody,TopAbs_EDGE,map.edgeMap);
            TopExp::MapShapes(curBody,TopAbs_VERTEX,map.vertexMap);
            MapOfBodyTopologyMap.insert(bodyIndex,map);

            //! experimental
            //bool isSweepable = GeomToolsClass::checkIfSweepable(curBody);
            //if(isSweepable)
            //    cout<<"____the body is sweepable____"<<endl;
        }
    }
    //! -----------------------
    //! shells: surface bodies
    //! -----------------------
    if(importSurfaceBodies)
    {
        for(int i=0; i<shells.length(); i++)
        {
            cout<<"____building the topology map for a shell____"<<endl;
            bodyIndex++;
            const TopoDS_Shape &curBody = shells.at(i);
            bodyMap.insert(bodyIndex,curBody);
            MapOfIsActive.insert(bodyIndex,true);
            TopologyMap map;
            TopExp::MapShapes(curBody,TopAbs_FACE,map.faceMap);
            TopExp::MapShapes(curBody,TopAbs_WIRE,map.wireMap);
            TopExp::MapShapes(curBody,TopAbs_EDGE,map.edgeMap);
            TopExp::MapShapes(curBody,TopAbs_VERTEX,map.vertexMap);
            MapOfBodyTopologyMap.insert(bodyIndex,map);
        }
    }

    //! ----------------------
    //! faces: surface bodies
    //! ----------------------
    if(importSurfaceBodies)
    {
        for(int i=0; i<faces.length(); i++)
        {
            cout<<"____building the topology map for a face____"<<endl;
            bodyIndex++;
            const TopoDS_Shape &curBody = faces.at(i);
            bodyMap.insert(bodyIndex,curBody);
            MapOfIsActive.insert(bodyIndex,true);
            TopologyMap map;
            TopExp::MapShapes(curBody,TopAbs_WIRE,map.wireMap);
            TopExp::MapShapes(curBody,TopAbs_EDGE,map.edgeMap);
            TopExp::MapShapes(curBody,TopAbs_VERTEX,map.vertexMap);
            MapOfBodyTopologyMap.insert(bodyIndex,map);
        }
    }
    //! -------------------
    //! wires: line bodies
    //! -------------------
    if(importLineBodies)
    {
        for(int i=0; i<wires.length(); i++)
        {
            cout<<"____building the topology map for a composite wire____"<<endl;
            bodyIndex++;
            const TopoDS_Shape &curBody = wires.at(i);
            bodyMap.insert(bodyIndex,curBody);
            MapOfIsActive.insert(bodyIndex,true);
            TopologyMap map;
            TopExp::MapShapes(curBody,TopAbs_WIRE,map.wireMap);
            TopExp::MapShapes(curBody,TopAbs_EDGE,map.edgeMap);
            TopExp::MapShapes(curBody,TopAbs_VERTEX,map.vertexMap);
            MapOfBodyTopologyMap.insert(bodyIndex,map);
        }
    }

    //! -------------------
    //! edges: line bodies
    //! -------------------
    if(importLineBodies)
    {
        for(int i=0; i<edges.length(); i++)
        {
            cout<<"____building the topology map for an edge____"<<endl;
            bodyIndex++;
            const TopoDS_Shape &curBody = edges.at(i);
            bodyMap.insert(bodyIndex,curBody);
            MapOfIsActive.insert(bodyIndex,true);
            TopologyMap map;
            TopExp::MapShapes(curBody,TopAbs_VERTEX,map.vertexMap);
            MapOfBodyTopologyMap.insert(bodyIndex,map);
        }
    }

    //! --------------------------------------------------------------------------
    //! activate the bodies
    //! at the begininning, all the shapes are active
    //! 28/09/2020 - please change it. The activation status is driven by the GUI
    //! --------------------------------------------------------------------------
    cout<<endl;
    for(QMap<int,TopoDS_Shape>::iterator it = bodyMap.begin(); it!=bodyMap.end(); ++it)
    {
        cout<<"____activating body: "<<it.key()<<"____"<<endl;
        MapOfIsActive.insert(it.key(),true);
    }
    cout<<endl;

    //for(QMap<int,TopologyMap>::iterator it = MapOfBodyTopologyMap.begin(); it!=MapOfBodyTopologyMap.end(); it++)
    //{
    //    TopologyMap aTopologyMap = it.value();
    //    cout<<"____body index "<<it.key()<<"____"<<endl;
    //    int NbFaces = aTopologyMap.faceMap.Extent();
    //    for(int n=1; n<=NbFaces; n++) cout<<"____the face nr: "<<n<<(aTopologyMap.faceMap.FindKey(n).IsNull()? " is not valid":" is valid")<<endl;
    //}
}

//! ------------------------------
//! function: constructor
//! details:  from a list of node
//! ------------------------------
geometryDataBase::geometryDataBase(const QList<SimulationNodeClass*> listOfNodes, const QString &archiveFileName, QObject *parent):QObject(parent)
{
    cout<<"geometryDataBase::geometryDataBase()->____constructor from a list of nodes____"<<endl;

    //! -----------------
    //! source file name
    //! -----------------
    fullSourceFileName = archiveFileName;

    //! -----------
    //! tree model
    //! -----------
    myModel = new QStandardItemModel(this);

    //! ---------------
    //! invisible root
    //! ---------------
    myInvisibleRootItem =  static_cast<QExtendedStandardItem*>(myModel->invisibleRootItem());

    //! -------------
    //! generic data
    //! -------------
    QVariant data;

    //! -----------------------------
    //! search for the geometry root
    //! -----------------------------
    for(int i=0; i<listOfNodes.length(); i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        if(curNode->getType()==SimulationNodeClass::nodeType_root)
        {
            cout<<"geometryDataBase::geometryDataBase()->____model root found____"<<endl;

            myRootItem = new QExtendedStandardItem();
            data.setValue(curNode);
            myRootItem->setData(data,Qt::UserRole);
            myRootItem->setData(curNode->getName(),Qt::DisplayRole);
            myRootItem->setEditable(false);
            myInvisibleRootItem->appendRow(myRootItem);

            //! --------------------------------------------------
            //! immediately apply backgound colors: post an event
            //! --------------------------------------------------
            QWidget *w = tools::getWidgetByName("maingwindow");
            if(w!=NULL)
            {
                int gradient = curNode->getPropertyItem("Gradient")->data(Qt::UserRole).value<Property>().getData().toInt();

                QVector<int> rgb1 = curNode->getPropertyItem("First color")->data(Qt::UserRole).value<Property>().getData().value<QVector<int>>();
                QColor firstColor(rgb1.at(0),rgb1.at(1),rgb1.at(2),255);

                QVector<int> rgb2 = curNode->getPropertyItem("Second color")->data(Qt::UserRole).value<Property>().getData().value<QVector<int>>();
                QColor secondColor(rgb2.at(0),rgb2.at(1),rgb2.at(2),255);

                ccout(QString("geometryDataBase::geometryDataBase()->____1st(%1, %2, %3)_2nd(%4, %5, %6)_type: %7____")
                      .arg(firstColor.red()).arg(firstColor.green()).arg(firstColor.blue())
                      .arg(secondColor.red()).arg(secondColor.green()).arg(secondColor.blue()).arg(gradient));

                QBackgroundEvent *event = new QBackgroundEvent(gradient, firstColor, secondColor);
                QApplication::postEvent(w,event);
            }
            break;
        }
    }

    //! ------------------------------------------------------
    //! "Import" item
    //! 22/12/2019 - commented: do not reload the import item
    //! ------------------------------------------------------
    for(int i=0;i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        if(curNode->getType()==SimulationNodeClass::nodeType_import)
        {
            cout<<"geometryDataBase::geometryDataBase()->____\"Import\" found____"<<endl;
            QExtendedStandardItem *importItem = new QExtendedStandardItem();
            data.setValue(curNode);
            importItem->setData(data,Qt::UserRole);
            importItem->setData(curNode->getName(),Qt::DisplayRole);
            importItem->setEditable(false);
            myRootItem->appendRow(importItem);
            break;
        }
    }

    //! -----------------
    //! "Geometry" items
    //! -----------------
    for(int i=0;i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        if(curNode->getType()==SimulationNodeClass::nodeType_geometry)
        {
            cout<<"geometryDataBase::geometryDataBase()->____Geometry root found____"<<endl;

            GeometryRootItem = new QExtendedStandardItem();
            data.setValue(curNode);

            GeometryRootItem->setData(data,Qt::UserRole);
            GeometryRootItem->setData(curNode->getName(),Qt::DisplayRole);
            GeometryRootItem->setEditable(false);
            myRootItem->appendRow(GeometryRootItem);

            //! -------------------------------------------------------------
            //! extract the TopoDS_Shape (compound) from the "Geometry" node
            //! the geometry is stored into the node as a TopoDS_shape
            //! -------------------------------------------------------------
            if(curNode->getPropertyItem("Source geometry")==NULL)
                cout<<"geometryDataBase::geometryDataBase()->____error in retrieving \"Source geometry\"____"<<endl;
            else cout<<"geometryDataBase::geometryDataBase()->____\"Source geometry\" found____"<<endl;

            cout<<"geometryDataBase::geometryDataBase()->____"<<curNode->getPropertyValue<QString>("Source geometry").toStdString()<<"_____"<<endl;

            if(!curNode->getPropertyItem("Source geometry")->data(Qt::UserRole).value<Property>().getData().canConvert<TopoDS_Shape>())
            {
                cout<<"geometryDataBase::geometryDataBase()->____bad data format in \"Source geometry\"____"<<endl;
            }

            TopoDS_Shape sourceGeometry = curNode->getPropertyValue<TopoDS_Shape>("Source geometry");
            if(!sourceGeometry.IsNull())
            {
                cout<<"geometryDataBase::geometryDataBase()->____shape valid____"<<endl;
            }
            else
            {
                cout<<"geometryDataBase::geometryDataBase()->____error in loading the shape____"<<endl;
            }

            //! ------------------------------
            //! initialize the private member
            //! ------------------------------
            myTopoDS_Shape = sourceGeometry;

            //! -----------------------
            //! build the topology map
            //! -----------------------
            this->buildMaps(myTopoDS_Shape);
            break;
        }
    }

    //! --------------------------
    //! attach the geometry items
    //! --------------------------
    for(int i=0; i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        if(curNode->getType()==SimulationNodeClass::nodeType_geometryBody)
        {
            int bodyIndex = curNode->getPropertyValue<int>("Map index");
            const TopoDS_Shape &theShape = bodyMap.value(bodyIndex);

            int mainShapeIndex, subShapeIndex;
            TopAbs_ShapeEnum subShapeType;
            this->getSubShapeNr(theShape,mainShapeIndex,subShapeIndex,subShapeType);
            GeometryTag loc;
            loc.isParent = false;
            loc.parentShapeNr = mainShapeIndex;
            loc.subTopNr = subShapeIndex;
            loc.subShapeType = subShapeType;
            std::vector<GeometryTag> vecLoc;
            vecLoc.push_back(loc);

            QVariant data;
            data.setValue(vecLoc);
            Property prop("Geometry",data,Property::PropertyGroup_Scope);
            curNode->replaceProperty("Geometry",prop);
            QExtendedStandardItem *geomItem = new QExtendedStandardItem();
            data.setValue(curNode);
            geomItem->setData(data,Qt::UserRole);
            geomItem->setData(curNode->getName(),Qt::DisplayRole);
            GeometryRootItem->appendRow(geomItem);
        }
    }

    //! ---------------------
    //! "Point masses" items
    //! ---------------------
    for(int i=0; i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        if(curNode->getType()==SimulationNodeClass::nodeType_pointMass)
        {
            cout<<"____ADDING POINT MASS____"<<endl;
            QExtendedStandardItem *pointMassItem= new QExtendedStandardItem();
            data.setValue(curNode);
            pointMassItem->setData(data,Qt::UserRole);
            pointMassItem->setData(curNode->getName(),Qt::DisplayRole);
            GeometryRootItem->appendRow(pointMassItem);
        }
    }

    //! --------------------------
    //! "Coordinate systems" root
    //! --------------------------
    for(int i=0;i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        if(curNode->getType()==SimulationNodeClass::nodeType_coordinateSystems)
        {
            CoordinateSystemsRootitem = new QExtendedStandardItem();
            data.setValue(curNode);

            CoordinateSystemsRootitem->setData(data,Qt::UserRole);
            CoordinateSystemsRootitem->setData(curNode->getName(),Qt::DisplayRole);
            myRootItem->appendRow(CoordinateSystemsRootitem);
            break;
        }
    }

    //! -------------------------------------
    //! attach the global coordinate systems
    //! -------------------------------------
    for(int i=0;i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        if(curNode->getType()==SimulationNodeClass::nodeType_coordinateSystem_global)
        {
            QExtendedStandardItem *item = new QExtendedStandardItem();
            data.setValue(curNode);

            item->setData(data,Qt::UserRole);
            item->setData(curNode->getName(),Qt::DisplayRole);
            CoordinateSystemsRootitem->appendRow(item);
            break;
        }
    }

    //! ----------------------------------------
    //! attach the remaining coordinate systems
    //! ----------------------------------------
    for(int i=0;i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        if(curNode->getType()==SimulationNodeClass::nodeType_coordinateSystem)
        {
            QExtendedStandardItem *item = new QExtendedStandardItem();
            data.setValue(curNode);

            item->setData(data,Qt::UserRole);
            item->setData(curNode->getName(),Qt::DisplayRole);
            CoordinateSystemsRootitem->appendRow(item);
        }
    }
}

//! -------------------------
//! function: getSubShapeNr1
//! details:  experimental
//! -------------------------
bool geometryDataBase::getSubShapeNr1(const TopoDS_Shape &aSubShape, int &mainShapeIndex, int &subShapeIndex, TopAbs_ShapeEnum &subShapeType) const
{
    for(QMap<int,TopoDS_Shape>::const_iterator it = bodyMap.cbegin(); it!=bodyMap.cend(); ++it)
    {
        int bodyIndex =it.key();
        const TopoDS_Shape &curShape = it.value();
        if(MapOfBodyTopologyMap.value(bodyIndex).solidMap.Contains(curShape))
        {
            mainShapeIndex = bodyIndex;
            subShapeIndex = MapOfBodyTopologyMap.value(bodyIndex).solidMap.FindIndex(aSubShape);
            subShapeType = aSubShape.ShapeType();
            return true;
        }
        if(MapOfBodyTopologyMap.value(bodyIndex).shellMap.Contains(curShape))
        {
            mainShapeIndex = bodyIndex;
            subShapeIndex = MapOfBodyTopologyMap.value(bodyIndex).shellMap.FindIndex(aSubShape);
            subShapeType = aSubShape.ShapeType();
            return true;
        }
        if(MapOfBodyTopologyMap.value(bodyIndex).faceMap.Contains(curShape))
        {
            mainShapeIndex = bodyIndex;
            subShapeIndex = MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindIndex(aSubShape);
            subShapeType = aSubShape.ShapeType();
            return true;
        }
        if(MapOfBodyTopologyMap.value(bodyIndex).wireMap.Contains(curShape))
        {
            mainShapeIndex = bodyIndex;
            subShapeIndex = MapOfBodyTopologyMap.value(bodyIndex).wireMap.FindIndex(aSubShape);
            subShapeType = aSubShape.ShapeType();
            return true;
        }
        if(MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Contains(curShape))
        {
            mainShapeIndex = bodyIndex;
            subShapeIndex = MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindIndex(aSubShape);
            subShapeType = aSubShape.ShapeType();
            return true;
        }
        if(MapOfBodyTopologyMap.value(bodyIndex).vertexMap.Contains(curShape))
        {
            mainShapeIndex = bodyIndex;
            subShapeIndex = MapOfBodyTopologyMap.value(bodyIndex).vertexMap.FindIndex(aSubShape);
            subShapeType = aSubShape.ShapeType();
            return true;
        }
    }
    return false;
}

//! -----------------------------
//! function: createGeometryRoot
//! details:
//! -----------------------------
void geometryDataBase::createGeometryRoot()
{
    cerr<<"geometryDataBase::createGeometryRoot()->____function called____"<<endl;

    //! -------------
    //! generic data
    //! -------------
    QVariant data;

    //! --------------------------------
    //! create the "Geometry" root node
    //! --------------------------------
    QVector<Property> props;

    //! ----------------------------
    //! property group "Definition"
    //! ----------------------------
    QString dataSourceLocation = fullSourceFileName;
    data.setValue(dataSourceLocation);
    Property property_source("Source",data,Property::PropertyGroup_Definition);
    props.push_back(property_source);

    Property::elementControl theElementControl = Property::elementControl_programControlled;    //! default value
    data.setValue(theElementControl);
    Property property_elementControl("Element control",data,Property::PropertyGroup_Definition);
    props.push_back(property_elementControl);

    //! -----------------------------------
    //! property "Scope" - the whole model
    //! -----------------------------------
    data.setValue(myTopoDS_Shape);
    Property property_sourceGeometry("Source geometry",data,Property::PropertyGroup_Scope);
    props.push_back(property_sourceGeometry);
    if(!myTopoDS_Shape.IsNull())
    {
        //! --------------------
        //! geometry properties
        //! --------------------
        geometryProperties_ myProp;
        computeGeometryProperties(myTopoDS_Shape,myProp);
        Property xB("Center of mass x",myProp.xB,Property::PropertyGroup_Properties);
        Property yB("Center of mass y",myProp.yB,Property::PropertyGroup_Properties);
        Property zB("Center of mass z",myProp.zB,Property::PropertyGroup_Properties);
        Property Ixx("Moment of inertia Ixx",myProp.Ixx,Property::PropertyGroup_Properties);
        Property Iyy("Moment of inertia Iyy",myProp.Iyy,Property::PropertyGroup_Properties);
        Property Izz("Moment of inertia Izz",myProp.Izz,Property::PropertyGroup_Properties);
        props.push_back(xB);
        props.push_back(yB);
        props.push_back(zB);
        props.push_back(Ixx);
        props.push_back(Iyy);
        props.push_back(Izz);

        //! ------------------
        //! bounding box size
        //! ------------------
        double LX, LY, LZ;
        getBoundingBox(myTopoDS_Shape,LX,LY,LZ);
        Property BB_sizeX("Length X",LX,Property::PropertyGroup_BoundingBox);
        Property BB_sizeY("Length Y",LY,Property::PropertyGroup_BoundingBox);
        Property BB_sizeZ("Length Z",LZ,Property::PropertyGroup_BoundingBox);
        props.push_back(BB_sizeX);
        props.push_back(BB_sizeY);
        props.push_back(BB_sizeZ);
    }

    //! -----
    //! node
    //! -----
    SimulationNodeClass *nodeGeometryRoot = new SimulationNodeClass("Geometry",SimulationNodeClass::nodeType_geometry,props,this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeGeometryRoot->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeGeometryRoot->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! --------------------------------
    //! create the "Geometry" root item
    //! --------------------------------
    GeometryRootItem = new QExtendedStandardItem();
    GeometryRootItem->setData("Geometry",Qt::DisplayRole);
    data.setValue(nodeGeometryRoot);
    GeometryRootItem->setData(data, Qt::UserRole);
    GeometryRootItem->setCheckState(Qt::Checked);

    //! -----------------------------------------------------
    //! append the geometry root item to the model root item
    //! -----------------------------------------------------
    myRootItem->appendRow(GeometryRootItem);
}

//! ------------------------------
//! function: createModelRootItem
//! details:
//! ------------------------------
void geometryDataBase::createModelRootItem()
{
    //! ------------------------------------------------
    //! create the QStandardItemModel "Model" root item
    //! ------------------------------------------------
    QVariant data;

    //! ---------------------
    //! vector of properties
    //! ---------------------
    QVector<Property> vecProps;

    //! -------------------------------------------
    //! property group "Background"
    //! 0 => uniform 1 => horizontal 2 => vertical
    //! -------------------------------------------
    Property property_gradient("Gradient",int(0),Property::PropertyGroup_Background);

    QColor firstColor(Qt::white);
    QVector<double> rgb;
    rgb.push_back(firstColor.redF());
    rgb.push_back(firstColor.greenF());
    rgb.push_back(firstColor.blueF());
    data.setValue(rgb);
    Property property_firstColor("First color",data,Property::PropertyGroup_Background);

    QColor secondColor(Qt::white);
    rgb.clear();
    rgb.push_back(secondColor.redF());
    rgb.push_back(secondColor.greenF());
    rgb.push_back(secondColor.blueF());
    data.setValue(rgb);
    Property property_secondColor("Second color",data,Property::PropertyGroup_Background);

    vecProps.push_back(property_gradient);
    vecProps.push_back(property_firstColor);
    vecProps.push_back(property_secondColor);

    //! --------------------------
    //! property group "Lighting"
    //! --------------------------
    Property property_ambient("Ambient",0.6,Property::PropertyGroup_Lighting);
    Property property_diffuse("Diffuse",0.1,Property::PropertyGroup_Lighting);
    Property property_specular("Specular",1.0,Property::PropertyGroup_Lighting);

    vecProps.push_back(property_ambient);
    vecProps.push_back(property_diffuse);
    vecProps.push_back(property_specular);

    SimulationNodeClass *nodeProjectRoot = new SimulationNodeClass("Model",SimulationNodeClass::nodeType_root,vecProps,this);
    data.setValue(nodeProjectRoot);

    //! ---------
    //! time tag
    //! ---------
    nodeProjectRoot->addTimeTag();

    myRootItem = new QExtendedStandardItem();
    myRootItem->setData("Model",Qt::DisplayRole);
    myRootItem->setData(data,Qt::UserRole);

    myInvisibleRootItem->appendRow(myRootItem);
}

//! ---------------------------
//! function: createImportItem
//! details:
//! ---------------------------
void geometryDataBase::createImportItem()
{
    QVariant data;

    //! vector of properties
    QVector<Property> vecProps;

    //! Source file path
    QString sourceFilePath;
    data.setValue(sourceFilePath);
    Property prop_sourceFilePath("Source file path",data,Property::PropertyGroup_Definition);
    vecProps.push_back(prop_sourceFilePath);

    //! Import solid bodies: true
    data.setValue(true);
    Property prop_importSolidBodies("Solid bodies",data,Property::PropertyGroup_Definition);
    vecProps.push_back(prop_importSolidBodies);

    //! Import surface bodies: true
    data.setValue(true);
    Property prop_importSurfaceBodies("Surface bodies",data,Property::PropertyGroup_Definition);
    vecProps.push_back(prop_importSurfaceBodies);

    //! Import line bodies: false
    data.setValue(false);
    Property prop_importLineBodies("Line bodies",data,Property::PropertyGroup_Definition);
    vecProps.push_back(prop_importLineBodies);

    //! geometry healing: true
    data.setValue(true);
    Property prop_geometryHealing("Geometry healing",data,Property::PropertyGroup_Definition);
    vecProps.push_back(prop_geometryHealing);

    //! ------------------------------------
    //! tolerance type: "Normal" by default
    //! 0 => "User defined"
    //! 1 => "Normal"
    //! 2 =" "Loose"
    //! ------------------------------------
    int toleranceType = 1;
    data.setValue(toleranceType);
    Property prop_toleranceType("Tolerance",data,Property::PropertyGroup_Definition);
    vecProps.push_back(prop_toleranceType);

    //! ----------------
    //! create the node
    //! ----------------
    SimulationNodeClass *nodeImport = new SimulationNodeClass("Import",SimulationNodeClass::nodeType_import,vecProps,this);
    data.setValue(nodeImport);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeImport->addTimeTag();
    //QString parentTimeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    //data.setValue(parentTimeTag);
    //nodeImport->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ----------------
    //! create the item
    //! ----------------
    myImportNode = new QExtendedStandardItem();
    myImportNode->setData(data,Qt::UserRole);
    data.setValue(QString("Import"));
    myImportNode->setData(data,Qt::DisplayRole);
    myImportNode->setEditable(false);

    //! -------------------------
    //! append the import otions
    //! -------------------------
    //GeometryRootItem->appendRow(item_import);
    myRootItem->appendRow(myImportNode);

}

//! ------------------
//! function: update
//! details:
//! ------------------
void geometryDataBase::update(const TopoDS_Shape &aShape)
{
    if(aShape.IsNull())
    {
        cerr<<"geometryDataBase::update()->____cannot update/init the geometry data base: the input shape is NULL____"<<endl;
        return;
    }
    cout<<"geometryDataBase::update()->____function called____"<<endl;

    //! ------------------------
    //! init the private member
    //! ------------------------
    myTopoDS_Shape = aShape;

    //! ----------------------------
    //! put the shape into the tree
    //! ----------------------------
    cout<<"geometryDataBase::update()->____updating the \"Geometry source\"____"<<endl;
    cout<<"geometryDataBase::update()->____source geometry status: "<<(!aShape.IsNull()? "OK":"not valid")<<"____"<<endl;

    SimulationNodeClass *geometryRoot = GeometryRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    data.setValue(aShape);
    Property prop("Source geometry",data,Property::PropertyGroup_Scope);
    geometryRoot->replaceProperty("Source geometry",prop);

    //! -------------------------
    //! create the topology maps
    //! -------------------------
    cout<<"geometryDataBase::update()->____build maps____"<<endl;
    this->buildMaps(aShape);

    //! --------------------------
    //! create the geometry items
    //! --------------------------
    cout<<"geometryDataBase::update()->____create geometry nodes____"<<endl;
    this->createGeometryNodes();
}

//! ------------------------
//! function: resetDataBase
//! details:
//! ------------------------
void geometryDataBase::resetDataBase()
{
    cout<<"geometryDataBase::resetDataBase()->____resetting the geometry database____"<<endl;

    //! ------------------------
    //! reset the internal maps
    //! ------------------------
    bodyMap.clear();
    MapOfBodyTopologyMap.clear();
    MapOfIsActive.clear();
    MapOfBodyNames.clear();

    myN1D = 0;
    myN2D = 0;
    myN3D = 0;

    bodyTwinFaces.clear();
    bodyTwinEdges.clear();

    //! --------------------------
    //! remove the geometry nodes
    //! --------------------------
    for(int n = GeometryRootItem->rowCount()-1; n>=0; n--) GeometryRootItem->removeRow(n);
}

//! ---------------------------------
//! function: removeBodyFromDataBase
//! details:  experimental
//! ---------------------------------
void geometryDataBase::removeBodyFromDataBase(int bodyIndex)
{
    cout<<"geometryDataBase::removeBodyFromDataBase()->____function called____"<<endl;

    //! -----------------------------
    //! update the geometry database
    //! -----------------------------
    bodyMap.remove(bodyIndex);
    MapOfBodyTopologyMap.remove(bodyIndex);
    MapOfBodyNames.remove(bodyIndex);
    MapOfIsActive.remove(bodyIndex);

    //! --------------------------
    //! - work not finished yet -
    //! --------------------------
    bodyTwinFaces.remove(bodyIndex);
    bodyTwinEdges.remove(bodyIndex);
    cout<<"geometryDataBase::removeBodyFromDataBase()->____final number of bodies: "<<bodyMap.size()<<"____"<<endl;

    //! ----------------------------------------
    //! rebuild and replace the geometry source
    //! ----------------------------------------
    cout<<"geometryDataBase::removeBodyFromDataBase()->____rebuilding the source geometry____"<<endl;
    TopoDS_Builder aBuilder;
    TopoDS_Compound aCompound;
    aBuilder.MakeCompound(aCompound);
    for(QMap<int,TopoDS_Shape>::iterator it=bodyMap.begin(); it!=bodyMap.end(); it++)
    {
        cout<<"geometryDataBase::removeBodyFromDataBase()->____adding shape____"<<endl;
        const TopoDS_Shape &curShape = it.value();
        aBuilder.Add(aCompound,curShape);
    }
    QVariant data;
    TopoDS_Shape aShape = aCompound;
    myTopoDS_Shape = aShape;
    data.setValue(aShape);
    Property prop_sourceGeometry("Source geometry",data,Property::PropertyGroup_Scope);

    SimulationNodeClass *node = GeometryRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    node->replaceProperty("Source geometry",prop_sourceGeometry);
    myTopoDS_Shape = aCompound;
}
