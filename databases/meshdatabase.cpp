//! ----------------
//! custom includes
//! ----------------
#include <meshdatabase.h>
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasourceface.h>
#include <edgemeshdatasourcebuildesclass.h>
#include "nodefactory.h"
#include "property.h"
#include "mydefines.h"
#include <meshtools.h>

//! ----
//! OCC
//! ----
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <MeshVS_DataSource.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_PGProps.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <MeshVS_MeshPrsBuilder.hxx>
#include <BRep_Tool.hxx>

//! ---
//! Qt
//! ---
#include <QDirIterator>
//! -----------------------------------------------
//! function - static: some characteristic lengths
//! details: 1) diagonal of the BB
//!          2) the shortest edge
//! -----------------------------------------------
void meshDataBase::shapeSize(const TopoDS_Shape &aShape, double &diag, double &minlength)
{
    if(!aShape.IsNull())
    {
        //! diagonal of the bounding box
        Bnd_Box shapeBB;
        BRepBndLib::Add(aShape, shapeBB);
        diag = sqrt(shapeBB.SquareExtent());

        //! minimum edge lenght
        GProp_GProps props;
        minlength = 1e99;
        for(TopExp_Explorer exp(aShape,TopAbs_EDGE);exp.More();exp.Next())
        {
            if(BRep_Tool::Degenerated(TopoDS::Edge(exp.Current()))) continue;
            BRepGProp::LinearProperties(exp.Current(),props);
            if(props.Mass()<minlength)minlength=props.Mass();
        }
    }
    else
    {
        diag = 0;
        minlength = 0;
    }
}
/*
//! --------------------------------
//! function: constructor - default
//! details:
//! --------------------------------
meshDataBase::meshDataBase(QObject *parent): geometryDataBase(parent)
{
    cerr<<"meshDataBase::meshDataBase()->____default constructor called____"<<endl;

    //! create the mesh root item, and append it to the root item
    this->createMeshRootNode();

    //! configure the "array-style" part of the database
    //this->initializeMeshArrays();
}
*/
//! ---------------------
//! function: destructor
//! details:
//! ---------------------
meshDataBase::~meshDataBase()
{
    cout<<"meshDataBase::~meshDataBase()->____DESTRUCTOR CALLED____"<<endl;
}

/*
//! --------------------------------
//! function: constructor - default
//! details:
//! --------------------------------
meshDataBase::meshDataBase(QObject *parent):geometryDataBase(parent)
{
    cerr<<"meshDataBase::meshDataBase()->____default constructor called____"<<endl;

    //! create the mesh root item, and append it to the "Model" root item
    this->createMeshRootNode();
}
*/

//! --------------------------
//! function: constructor III
//! details:
//! --------------------------
meshDataBase::meshDataBase(const QList<SimulationNodeClass*> listOfNodes, const QString &archiveFileName, QObject *parent):
    geometryDataBase(listOfNodes, archiveFileName, parent)
{
    cout<<"meshDataBase::meshDataBase()->____constructor III called____"<<endl;

    //! ----------------------------
    //! retrieve the mesh root node
    //! ----------------------------
    QVariant data;
    for(int i=0; i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        SimulationNodeClass::nodeType theType, theFamily;
        theType = curNode->getType();
        theFamily = curNode->getFamily();
        if(theType==SimulationNodeClass::nodeType_meshControl)
        {
            //! the "Mesh" root node has been found
            MeshRootItem = new QExtendedStandardItem();
            data.setValue(curNode);
            MeshRootItem->setData(data,Qt::UserRole);
            MeshRootItem->setData(curNode->getName(),Qt::DisplayRole);
            myRootItem->appendRow(MeshRootItem);
            cout<<"meshDataBase::meshDataBase()->____mesh root found____"<<endl;
            break;
        }
    }
    for(int i=0; i<listOfNodes.length();i++)
    {
        SimulationNodeClass *curNode = listOfNodes.at(i);
        SimulationNodeClass::nodeType theType, theFamily;
        theType = curNode->getType();
        theFamily = curNode->getFamily();
        if(theFamily==SimulationNodeClass::nodeType_meshControl && theType!=SimulationNodeClass::nodeType_meshControl)
        {
            //! a mesh control node has been found
            QExtendedStandardItem *item = new QExtendedStandardItem();
            data.setValue(curNode);
            item->setData(data, Qt::UserRole);
            item->setData(curNode->getName(),Qt::DisplayRole);
            MeshRootItem->appendRow(item);
        }
    }

    //! configure the "array-style" part of the database
    this->initializeMeshArrays();

    //! ----------------------------------------
    //! access the folder containing the meshes
    //! ----------------------------------------
    QString project_filesDir = fullSourceFileName;
    project_filesDir.chop(4);
    project_filesDir.append("_files");

    QDir curDir(project_filesDir);
    curDir.cd("MDB");

    cout<<"---------------------------------------------------------"<<endl;
    cout<<"- begin reading mesh database"<<endl;
    cout<<"- "<<curDir.absolutePath().toStdString()<<endl;
    cout<<"---------------------------------------------------------"<<endl;

    //! --------------------------
    //! Read the main body meshes
    //! --------------------------
    curDir.cd("Volume");
    QDirIterator *dirIterator = new QDirIterator(curDir.absolutePath(),QDirIterator::NoIteratorFlags);
    while(dirIterator->hasNext())
    {
        QString meshFileName = dirIterator->next();
        QFileInfo info(meshFileName);
        if(!info.isDir())
        {
            cout<<"Reading volume mesh file: "<<meshFileName.toStdString()<<endl;

            int bodyIndex;
            QString item = meshFileName.split("/").last();
            sscanf(item.toStdString().c_str(),"%d",&bodyIndex);

            //! rebuild the mesh data source by calling the constructor
            occHandle(Ng_MeshVS_DataSource3D) meshDS = new Ng_MeshVS_DataSource3D(meshFileName.toStdString());
            if(meshDS.IsNull())
            {
                cout<<"Volume mesh nr: "<<bodyIndex<<" is null"<<endl;
                ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, true);
            }
            ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, false);
            ArrayOfMeshDS.insert(bodyIndex, meshDS);
            //{
            //    ArrayOfMeshDS.insert(bodyIndex, meshDS);
            //   ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, false);
            //}
            //else
            //{
            //    cout<<"Volume mesh nr: "<<bodyIndex<<" is null"<<endl;
            //occHandle(MeshVS_DataSource) nullMeshHandle;
            //ArrayOfMeshDS.insert(bodyIndex, occHandle(MeshVS_DataSource)());
            //ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, true);
            //}
        }
        //emit operationReadingFromDiskDone();
    }
    delete dirIterator;

    //! ----------------------------------
    //! Read the main surface body meshes
    //! ----------------------------------
    curDir.cdUp();
    curDir.cd("Surface");
    QDirIterator *dirIterator1 = new QDirIterator(curDir.absolutePath(),QDirIterator::NoIteratorFlags);
    while(dirIterator1->hasNext())
    {
        QString meshFileName = dirIterator1->next();
        QFileInfo info(meshFileName);
        if(!info.isDir())
        {
            cout<<"Reading surface mesh file: "<<meshFileName.toStdString()<<"<-"<<endl;
            QString item = meshFileName.split("/").last();
            int bodyIndex;
            sscanf(item.toStdString().c_str(),"%d_2D",&bodyIndex);
            //! rebuild the mesh data source by calling the constructor
            occHandle(Ng_MeshVS_DataSource2D) meshDS = new Ng_MeshVS_DataSource2D(meshFileName.toStdString());

            if(!meshDS.IsNull())
            {
                ArrayOfMeshDS2D.insert(bodyIndex, meshDS);
            }
            else
            {
                occHandle(MeshVS_DataSource) nullMeshHandle;
                ArrayOfMeshDS2D.insert(bodyIndex,nullMeshHandle);
            }
        }
        //emit operationReadingFromDiskDone();
    }
    delete dirIterator1;

    //! ----------------------------------------
    //! enter the folder containing the meshes
    //! ----------------------------------------

    curDir.cdUp();
    curDir.cd("Faces");

    //! ---------------------------------------
    //! Reading the subtopology meshes - faces
    //! ---------------------------------------
    QDirIterator *dirIterator2 = new QDirIterator(curDir.absolutePath(),QDirIterator::NoIteratorFlags);
    while(dirIterator2->hasNext())
    {
        QString meshFileName = dirIterator2->next();
        QFileInfo info(meshFileName);
        if(!info.isDir())
        {
            cout<<"meshDataBase::meshDataBase()->____Reading face mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
            QList<QString> List;
            List = meshFileName.split("/");
            QString item = List.last();

            List.clear();
            List = item.split("_");
            QString item1 = List.first();
            int bodyIndex = item1.toInt();
            QString item2 = List.last();
            int subTopologyNr = item2.toInt();

            //! --------------------------------------------------------
            //! rebuild the mesh data source by calling the constructor
            //! --------------------------------------------------------
            occHandle(Ng_MeshVS_DataSourceFace) meshDS = new Ng_MeshVS_DataSourceFace(meshFileName.toStdString());
            ArrayOfMeshDSOnFaces.setValue(bodyIndex, subTopologyNr, meshDS);
        }
    }
    delete dirIterator2;

    //!cout<<"meshDataBase::meshDataBase()->____current dir: "<<curDir.path().toStdString()<<"____"<<endl;
    curDir.cdUp();
    //!cout<<"meshDataBase::meshDataBase()->____current dir: "<<curDir.path().toStdString()<<"____"<<endl;
    curDir.cdUp();
    //!cout<<"meshDataBase::meshDataBase()->____current dir: "<<curDir.path().toStdString()<<"____"<<endl;
    curDir.cdUp();
    //!cout<<"meshDataBase::meshDataBase()->____current dir: "<<curDir.path().toStdString()<<"____"<<endl;

    //! ------------------------------
    //! Rebuild the 1D meshes - edges
    //! ------------------------------
    //edgeMeshDataSourceBuildesClass *the1DMeshBuilder = new edgeMeshDataSourceBuildesClass(this,this);
    //the1DMeshBuilder->perform();
}

//! -------------------------
//! function: constructor II
//! details:
//! -------------------------
meshDataBase::meshDataBase(const TopoDS_Shape &shape, const QString& theFilePath, QObject *parent):
    geometryDataBase(shape,theFilePath,parent)
{
    cout<<"meshDataBase::meshDataBase()->____CONSTRUCTOR II CALLED____"<<endl;

    //! ------------------------------------------------------------------
    //! create the mesh root item, and append it to the "Model" root item
    //! ------------------------------------------------------------------
    this->createMeshRootNode();

    //! -------------------------------------------------
    //! configure the "array-style" part of the database
    //! -------------------------------------------------
    if(shape.IsNull())
    {
        cout<<"meshDataBase::meshDataBase()->____cannot initialize the mesh data base: shape is NULL____"<<endl;
    }
    else this->initializeMeshArrays();
}

//! ------------------------------
//! function: createStandardModel
//! details:
//! ------------------------------
void meshDataBase::createStandardModel()
{
    geometryDataBase::createStandardModel();
    this->createMeshRootNode();
}

//! -----------------
//! function: update
//! details:
//! -----------------
void meshDataBase::update(const TopoDS_Shape &aShape)
{
    if(aShape.IsNull())
    {
        cerr<<"meshDataBase::meshDataBase()->____cannot update/init the mesh data base: the input shape is NULL____"<<endl;
        return;
    }
    //! ------------------------------------------
    //! update the unredrlying geometry data base
    //! ------------------------------------------
    geometryDataBase::update(aShape);

    //! -------------------------------------------------
    //! configure the "array-style" part of the database
    //! -------------------------------------------------
    this->initializeMeshArrays();
}

//! --------------------------------------------
//! function: calculateDefaultMeshingParameters
//! details:
//! --------------------------------------------
void meshDataBase::calculateDefaultMeshingParameters()
{
    //! ------------------------------------------
    //! calculate the default meshing paramenters
    //! ------------------------------------------
    for(QMap<int,TopoDS_Shape>::iterator it = bodyMap.begin(); it!=bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();
        const TopoDS_Shape &theShape = it.value();

        //! retrieve the typical dimension of the current body
        double diag, minlength;
        this->shapeSize(theShape, diag, minlength);

        //! -------------------------------
        //! the default meshing parameters
        //! -------------------------------
        default2DMeshEngine.insert(bodyIndex, Property::meshEngine2D_Netgen);
        default3DMeshEngine.insert(bodyIndex, Property::meshEngine3D_Netgen);
        defaultMeshOrder.insert(bodyIndex, Property::meshOrder_First);
        defaultVolumeGrading.insert(bodyIndex,0.2);
        defaultVolumeMinElementSize.insert(bodyIndex,0.1);
        defaultVolumeMaxElementSize.insert(bodyIndex,diag/20);
        defaultSmoothingSteps.insert(bodyIndex,0);              //! smoothing steps
        defaultIsElementStraight.insert(bodyIndex,false);       //! init isoparametric elements
        defaultIsBodyHealingOn.insert(bodyIndex,false);         //! init no mesh-based healing for all the bodies
        defaultHasPrismaticFaces.insert(bodyIndex,false);       //! by default a body has not inflation
        ArrayOfDefaultMeshType.insert(bodyIndex,0);             //! full tri-tet mesh by default

        //! ---------------------------------------
        //! feature preserving and mesh correction
        //! ---------------------------------------
        defaultMapOfFeaturePreserving.insert(bodyIndex,false);
        defaultMapOfGeometryCorrection.insert(bodyIndex,false);
        featuredTags.clear();

        //! -------------------------------
        //! the default meshing parameters
        //! -------------------------------
        defaultMapOfUseBRep.insert(bodyIndex,true);
        default2DMeshEngine.insert(bodyIndex, Property::meshEngine2D_Netgen);
        default3DMeshEngine.insert(bodyIndex, Property::meshEngine3D_Netgen);
        defaultMeshOrder.insert(bodyIndex, Property::meshOrder_First);
        defaultVolumeGrading.insert(bodyIndex,0.2);
        defaultVolumeMinElementSize.insert(bodyIndex,0.1);
        defaultVolumeMaxElementSize.insert(bodyIndex,diag/20);
        defaultSmoothingSteps.insert(bodyIndex,0);
        defaultIsElementStraight.insert(bodyIndex,false);
        defaultIsBodyHealingOn.insert(bodyIndex,false);
        defaultIsBodyDefeaturingOn.insert(bodyIndex,false);
        defaultIsMeshingRunningInMemory.insert(bodyIndex,true);

        defaultMapOfDiscretizer.insert(bodyIndex,Property::meshEngine2D_OCC_STL);
        defaultMapOfAngularDeflection.insert(bodyIndex,0.01);
        defaultMapOfLinearDeflection.insert(bodyIndex,0.1);
        defaultMinFaceSize.insert(bodyIndex,0.1);
        defaultMaxFaceSize.insert(bodyIndex,1000);

        //! ---------------------------
        //! default TetWild parameters
        //! ---------------------------
        defaultMapOfEnvelopeSizingType.insert(bodyIndex,0);                 //! "Relative" by default
        defaultMapOfRelativeEnvelopeSize.insert(bodyIndex,0.0025);        //! "0.0025" by default
        defaultMapOfAbsoluteEnvelopeSize.insert(bodyIndex,1.0);           //! "1.0" 1 mm by default

        defaultMapOfIdealLengthSizingType.insert(bodyIndex,0);              //! "Relative" by default
        defaultMapOfIdealLengthRelativeSize.insert(bodyIndex,0.10);         //! "0.10" 10% by default
        defaultMapOfIdealLengthAbsoluteSize.insert(bodyIndex,100.0);        //! "100.0" 100 mm by default
    }
}

//! ------------------------------------------------
//! function: apply the default meshing paramenters
//! details:
//! ------------------------------------------------
void meshDataBase::applyDefaultMeshingParameters()
{
    //! ----------------------------------
    //! Total number of bodies (3D+2D+1D)
    //! ----------------------------------
    int Nt = N3D()+N2D()+N1D();

    //! ----------------------------------
    //! initialize the meshing parameters
    //! ----------------------------------
    for(QMap<int,TopoDS_Shape>::iterator it = bodyMap.begin(); it!=bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();

        //! ----------------------------
        //! init the status of the mesh
        //! ----------------------------
        ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, Standard_True);

        //! ----------------------------------
        //! initialize the meshing parameters
        //! ----------------------------------
        mapOfUseBRep.insert(bodyIndex, defaultMapOfUseBRep.value(bodyIndex));
        ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, true);
        ArrayOfMesh2DEngine.insert(bodyIndex,default2DMeshEngine.value(bodyIndex));
        ArrayOfMesh3DEngine.insert(bodyIndex,default3DMeshEngine.value(bodyIndex));
        ArrayOfMeshOrder.insert(bodyIndex,defaultMeshOrder.value(bodyIndex));
        ArrayOfGradingValue.insert(bodyIndex,defaultVolumeGrading.value(bodyIndex));
        ArrayOfMinBodyElementSize.insert(bodyIndex,defaultVolumeMinElementSize.value(bodyIndex));
        ArrayOfMaxBodyElementSize.insert(bodyIndex,defaultVolumeMaxElementSize.value(bodyIndex));
        ArrayOfSmoothingSteps.insert(bodyIndex,defaultSmoothingSteps.value(bodyIndex));
        mapOfIsElementStraight.insert(bodyIndex,defaultIsElementStraight.value(bodyIndex));
        ArrayOfMeshType.insert(bodyIndex,ArrayOfDefaultMeshType.value(bodyIndex));

        MapOfIsBodyDefeaturingOn.insert(bodyIndex,defaultIsBodyDefeaturingOn.value(bodyIndex));
        MapOfIsBodyHealingOn.insert(bodyIndex,defaultIsBodyHealingOn.value(bodyIndex));
        MapOfIsMeshingRunningInMemory.insert(bodyIndex,defaultIsMeshingRunningInMemory.value(bodyIndex));
        MapOfIsMeshSimplificationOn.insert(bodyIndex,defaultIsMeshSimplificationOn.value(bodyIndex));
        MapOfMeshSimplificationBy.insert(bodyIndex,defaultMeshDefeaturingBy.value(bodyIndex));
        MapOfBodyDefeaturingParameterValue.insert(bodyIndex,defaultMeshDefeaturingParameterValue.value(bodyIndex));

        mapOfDiscretizer.insert(bodyIndex,defaultMapOfDiscretizer.value(bodyIndex));
        mapOfAngularDeflection.insert(bodyIndex,defaultMapOfAngularDeflection.value(bodyIndex));
        mapOfLinearDeflection.insert(bodyIndex,defaultMapOfLinearDeflection.value(bodyIndex));
        mapOfMinFaceSize.insert(bodyIndex,defaultMinFaceSize.value(bodyIndex));
        mapOfMaxFaceSize.insert(bodyIndex,defaultMaxFaceSize.value(bodyIndex));

        //! ----------------------------------------------------------------------------------------------
        //! feeature preserving flag
        //! geometry correction flag (mesh points projection onto geometry for a patch independent method
        //! ----------------------------------------------------------------------------------------------
        mapOfFeaturePreserving.insert(bodyIndex,defaultMapOfFeaturePreserving.value(bodyIndex));
        mapOfGeometryCorrection.insert(bodyIndex,defaultMapOfGeometryCorrection.value(bodyIndex));
    }

    //! -----------------------------------------------------------------------------------
    //! each face mesh sizing value is initialized with a size and "not modified" flag [*]
    //! -----------------------------------------------------------------------------------
    int NMaxFace =0;
    int NMaxEdge =0;
    int NMaxVertex = 0;
    this->NmaxEntity(NMaxFace, NMaxEdge, NMaxVertex);
    for(int bodyIndex=1; bodyIndex<=Nt; bodyIndex++)
    {
        double elementSizeOnFace = 1;   //! see the note [*]
        for(int faceNr=1; faceNr<=NMaxFace; faceNr++)
        {
            //! the element size on face has no effect if myArrayOfIsFaceModified->Value(faceNr) == Standard_False
            MapOfIsFaceModified.setValue(bodyIndex,faceNr,false);
            MapOfElementSizeOnFace.setValue(bodyIndex,faceNr,elementSizeOnFace);
        }
        double elementSizeOnEdge = 1;   //! see the note below [**]
        for(int edgeNr=1; edgeNr<=NMaxEdge; edgeNr++)
        {
            //! the element size on edge has no effect if myArrayOfIsEdgeModified->Value(edgeNr) == Standard_False
            MapOfIsEdgeModified.setValue(bodyIndex,edgeNr,false);
            MapOfElementSizeOnEdge.setValue(bodyIndex,edgeNr,elementSizeOnEdge);
            MapOfSizingTypeOnEdge.setValue(bodyIndex,edgeNr,Standard_True);
        }
        double elementSizeOnVertex = 1; //! see the note below [***]
        for(int vertexNr=1; vertexNr<=NMaxVertex; vertexNr++)
        {
            //! the element size around the vertex has no effect if myArrayOfIsVertexModified->Value(edgeNr) == Standard_False
            MapOfIsEdgeModified.setValue(bodyIndex,vertexNr,false);
            MapOfElementSizeOnVertex.setValue(bodyIndex,vertexNr,elementSizeOnVertex);
            MapOfVertexPinball.setValue(bodyIndex,vertexNr,elementSizeOnVertex*10.0);
        }
    }

    //! --------------------------------------------------
    //! the same of the previous for the new meshdatabase
    //! --------------------------------------------------
    for(QMap<int,TopoDS_Shape>::const_iterator it = bodyMap.cbegin(); it!=bodyMap.cend(); ++it)
    {
        int bodyIndex = it.key();
        int NbFaces = MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
        int NbEdges = MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent();
        int NbVertices= MapOfBodyTopologyMap.value(bodyIndex).vertexMap.Extent();
        for(int faceNr=1; faceNr<=NbFaces; faceNr++)
        {
            MapOfIsFaceModified.setValue(bodyIndex,faceNr,false);
            MapOfElementSizeOnFace.setValue(bodyIndex,faceNr,1.0);
        }

        for(int edgeNr=1; edgeNr<=NbEdges; edgeNr++)
        {
            MapOfIsEdgeModified.setValue(bodyIndex,edgeNr,false);
            MapOfSizingTypeOnEdge.setValue(bodyIndex,edgeNr,true);
            MapOfElementSizeOnEdge.setValue(bodyIndex,edgeNr,1.0);
        }

        for(int vertexNr=1; vertexNr<=NbVertices; vertexNr++)
        {
            MapOfIsVertexModified.setValue(bodyIndex,vertexNr,false);
            MapOfVertexPinball.setValue(bodyIndex,vertexNr,1.0);
            MapOfElementSizeOnVertex.setValue(bodyIndex,vertexNr,1.0);
        }
    }

    //! -----------------------------------------------
    //! set the default settings for inflation on body
    //! By default no body has inflation, the list of
    //! faces having inflation is empty for each body
    //! -----------------------------------------------
    for(QMap<int,TopoDS_Shape>::iterator it = bodyMap.begin(); it!=bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();
        int NbFacesOnBody = MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
        HasPrismaticFaces.insert(bodyIndex,defaultHasPrismaticFaces.value(bodyIndex));
        //!cout<<"____body nr: "<<it.key()<<" "<<(defaultHasPrismaticFaces.value(bodyIndex)==true? "has":"has not")<<" inflation____"<<endl;

        //! ------------------------------------------
        //! empty, uninitialized list of face numbers
        //! ------------------------------------------
        QList<int> listOfFacesHavingInflation = QList<int>();
        for(int faceNr = 0; faceNr<=NbFacesOnBody; faceNr++) prismaticFaces.insert(bodyIndex,listOfFacesHavingInflation);
    }

    //! ----------------------------------------------
    //! default parameters for the new TetWild mesher
    //! ----------------------------------------------
    for(QMap<int,TopoDS_Shape>::iterator it = bodyMap.begin(); it!=bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();
        mapOfEnvelopeSizingType.insert(bodyIndex,defaultMapOfEnvelopeSizingType.value(bodyIndex));
        mapOfRelativeEnvelopeSize.insert(bodyIndex,defaultMapOfRelativeEnvelopeSize.value(bodyIndex));
        mapOfAbsoluteEnvelopeSize.insert(bodyIndex,defaultMapOfAbsoluteEnvelopeSize.value(bodyIndex));

        mapOfIdealLengthSizingType.insert(bodyIndex,defaultMapOfIdealLengthSizingType.value(bodyIndex));
        mapOfIdealLengthRelativeSize.insert(bodyIndex,defaultMapOfIdealLengthRelativeSize.value(bodyIndex));
        mapOfIdealLengthAbsoluteSize.insert(bodyIndex,defaultMapOfIdealLengthAbsoluteSize.value(bodyIndex));
    }
}

//! ------------------------------------------------------------
//! function: max number of subtopologies
//! details:  find the absolute maximum number of faces, edges,
//!           among the bodies
//! ------------------------------------------------------------
void meshDataBase::NmaxEntity(int &NMaxFace, int &NMaxEdge, int &NMaxVertex)
{
    for(int k=1; k<=MapOfBodyTopologyMap.size(); k++)
    {
        int NFace = MapOfBodyTopologyMap.value(k).faceMap.Extent();
        int NEdge = MapOfBodyTopologyMap.value(k).edgeMap.Extent();
        int NVertex = MapOfBodyTopologyMap.value(k).vertexMap.Extent();

        if(NMaxFace<=NFace) NMaxFace=NFace;
        if(NMaxVertex<=NVertex) NMaxVertex=NVertex;
        if(NMaxEdge<=NEdge) NMaxEdge=NEdge;
    }
}

//! -----------------------------
//! function: createMeshRootNode
//! details:
//! -----------------------------
void meshDataBase::createMeshRootNode()
{
    //! --------------------------
    //! create the mesh root item
    //! --------------------------
    QVariant data;
    data.setValue(int(0));
    Property property_NN("Number of nodes",data,Property::PropertyGroup_Statistics);

    data.setValue(int(0));
    Property property_NE("Number of elements",data,Property::PropertyGroup_Statistics);

    data.setValue(false);
    Property property_showMeshNodes("Show mesh nodes",data,Property::PropertyGroup_MeshViewOptions);

    //! initial size seed
    data.setValue(0);
    Property property_initialSizeSeed("Initial size seed",data,Property::PropertyGroup_Sizing);

    //! relevance
    int relevance = 0;
    data.setValue(relevance);
    Property property_relevance("Relevance",data,Property::PropertyGroup_Sizing);

    //! min edge length
    double diag,ml;
    this->shapeSize(this->shape(),diag,ml);
    data.setValue(ml);
    Property property_minedgelen("Minimum edge length",data,Property::PropertyGroup_Sizing);

    //! 0 => low = 1 step
    //! 1 => medium = 5 steps
    //! 2 => high = 10 steps
    data.setValue(int(0));
    Property property_smoothing("Smoothing",data,Property::PropertyGroup_Sizing);

    //! submeshes generation: true => deferred false => forward
    data.setValue(false);
    Property property_submeshes("Submeshes",data,Property::PropertyGroup_Advanced);

    //! curved elements
    data.setValue(false);
    Property property_curvedElements("Straight sided elements",data, Property::PropertyGroup_Advanced);
    QVector<Property> props;

    //! mesh order - element midside nodes
    //! 0 => program controlled - 1 => first order - 2 => second order
    data.setValue(0);
    Property property_elementMidsideNodes("Element midside nodes", data, Property::PropertyGroup_Sizing);

    //! -----------------------------------------
    //! run in memory option - default in memory
    //! -----------------------------------------
    //data.setValue(false);
    //Property property_runInMemory("Run in memory",data,Property::PropertyGroup_Advanced);

    props.push_back(property_NN);
    props.push_back(property_NE);
    props.push_back(property_showMeshNodes);
    props.push_back(property_relevance);
    props.push_back(property_initialSizeSeed);
    props.push_back(property_smoothing);
    props.push_back(property_minedgelen);
    props.push_back(property_submeshes);
    //props.push_back(property_surfaceMesher);
    props.push_back(property_elementMidsideNodes);
    props.push_back(property_curvedElements);
    //props.push_back(property_runInMemory);

    SimulationNodeClass *nodeMeshRoot = new SimulationNodeClass("Mesh",SimulationNodeClass::nodeType_meshControl,props,this);

    //! --------------------------
    //! create the mesh root node
    //! --------------------------
    MeshRootItem = new QExtendedStandardItem();
    MeshRootItem->setData("Mesh", Qt::DisplayRole);
    data.setValue(nodeMeshRoot);
    MeshRootItem->setData(data, Qt::UserRole);
    MeshRootItem->setEditable(false);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeMeshRoot->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeMeshRoot->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! -------------------------------------------
    //! append the mesh root item to the root item
    //! -------------------------------------------
    myRootItem->appendRow(MeshRootItem);
}

//!--------------------------------
//! function: initializeMeshArrays
//! details:
//! -------------------------------
void meshDataBase::initializeMeshArrays()
{
    //! -----------------------------------------
    //! allocate the memory for the data sources
    //! -----------------------------------------
    int Nt = N3D()+N2D()+N1D();
    for(int i=1;i<=Nt;i++)
    {
        ArrayOfMeshDS.insert(i,occHandle(MeshVS_DataSource)());
        ArrayOfMeshDS2D.insert(i,occHandle(MeshVS_DataSource)());
    }

    //! ----------------------------------------------------------
    //! max number of vertices, edges, faces among all the bodies
    //! ----------------------------------------------------------
    int NMaxVertex=0;
    int NMaxEdge=0;
    int NMaxFace=0;
    this->NmaxEntity(NMaxFace, NMaxEdge, NMaxVertex);

    //! -----------------------------------------------------------------------------------
    //! local mesh data sources
    //! faces: the stl repairing tool can add triangles to the surface mesh, in order
    //! to fill the holes. These triangles are labeled with "0", and a a MeshVS_DataSource
    //! is built: this is put into the array ArrayOfMeshDSOnFaces at position "0"
    //! -----------------------------------------------------------------------------------

    //! --------------------------------------------------
    //! calculate and apply the default meshing parameter
    //! --------------------------------------------------
    this->calculateDefaultMeshingParameters();
    this->applyDefaultMeshingParameters();
}

//! --------------------------------------------------
//! function: getMeshingParameters
//! details:  get the meshing parameters for the body
//! --------------------------------------------------
meshParam meshDataBase::getMeshingParameters(int bodyIndex)
{
    cout<<"meshDataBase::getMeshingParameters()->____function called____"<<endl;

    //! --------------------------------
    //! retrieve the meshing parameters
    //! --------------------------------
    meshParam mp;

    //! ---------------
    //! handling faces
    //! ---------------
    for(int faceNr=1; faceNr<=this->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent(); faceNr++)
    {
        if(this->MapOfIsFaceModified.getValue(bodyIndex,faceNr)==true)
        {
            double faceSize = this->MapOfElementSizeOnFace.getValue(bodyIndex,faceNr);
            mp.faceSize.insert(faceNr,faceSize);
        }
    }

    //! ---------------
    //! handling edges
    //! ---------------
    for(int edgeNr=1; edgeNr<=this->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent(); edgeNr++)
    {
        if(this->MapOfIsEdgeModified.getValue(bodyIndex,edgeNr)==true)
        {
            int sizingType = MapOfSizingTypeOnEdge.getValue(bodyIndex,edgeNr);
            double edgeElementSize;
            switch(sizingType)
            {
            case 0:
                //! --------------------------------------------------
                //! the edge element size has been directly specified
                //! --------------------------------------------------
                edgeElementSize = this->MapOfElementSizeOnEdge.getValue(bodyIndex,edgeNr);
                break;

            case 1:
                //! ---------------------------------------------------------------------
                //! the edge element size has been specified through number of divisions
                //! ---------------------------------------------------------------------
            {
                int numberOfDivisionsOnEdge = this->MapOfNumberOfDivisionOnEdge.getValue(bodyIndex,edgeNr);
                const TopoDS_Shape &theCurEdge = this->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(edgeNr);
                GProp_GProps LProps;
                BRepGProp::LinearProperties(theCurEdge,LProps);
                double L = LProps.Mass();
                edgeElementSize = L/numberOfDivisionsOnEdge;
            }
                break;
            }
            mp.edgeSize.insert(edgeNr,edgeElementSize);
        }
    }

    //! ------------------
    //! handling vertexes
    //! ------------------
    for(int vertexNr=1; vertexNr<=this->MapOfBodyTopologyMap.value(bodyIndex).vertexMap.Extent(); vertexNr++)
    {
        if(this->MapOfIsVertexModified.getValue(bodyIndex,vertexNr)==true)
        {
            double vertexSize = this->MapOfElementSizeOnVertex.getValue(bodyIndex,vertexNr);
            double pinball = this->MapOfVertexPinball.getValue(bodyIndex,vertexNr);
            mp.vertexSize.insert(vertexNr,vertexSize);
            mp.pinballSize.insert(vertexNr,pinball);
        }
    }

    //! ----------------------------------
    //! handling global (at "body" level)
    //! ----------------------------------
    mp.grading = this->ArrayOfGradingValue.value(bodyIndex);
    mp.minElementSize = this->ArrayOfMinBodyElementSize.value(bodyIndex);
    mp.maxElementSize = this->ArrayOfMaxBodyElementSize.value(bodyIndex);
    Property::meshOrder mo1 = this->ArrayOfMeshOrder.value(bodyIndex);
    mp.isSecondOrder = (mo1 == Property::meshOrder_Second)? true: false;
    mp.isElementStraight = this->mapOfIsElementStraight.value(bodyIndex);

    return mp;
}

//! ------------------------
//! function: resetDataBase
//! details:
//! ------------------------
void meshDataBase::resetDataBase()
{
    //! ----------------------------
    //! reset the geometry database
    //! ----------------------------
    geometryDataBase::resetDataBase();

    //! ------------------------------
    //! face/edge/vertex element size
    //! ------------------------------
    MapOfElementSizeOnFace.clear();
    MapOfElementSizeOnEdge.clear();
    MapOfNumberOfDivisionOnEdge.clear();
    MapOfSizingTypeOnEdge.clear();
    MapOfElementSizeOnVertex.clear();
    MapOfVertexPinball.clear();

    //! ----------------
    //! prismatic faces
    //! ----------------
    HasPrismaticFaces.clear();
    defaultHasPrismaticFaces.clear();

    //! ---------------------------------------------------
    //! key => bodyIndex
    //! value => QList<int> list of faces having inflation
    //! ---------------------------------------------------
    prismaticFaces.clear();

    //! ----------------------------------------------------------------------------------
    //! key => bodyIndex
    //! value => QList<prismaticLayerParameters> prismatic layers parameters on that body
    //! ----------------------------------------------------------------------------------
    prismaticMeshParameters.clear();

    //! ------------------------------------
    //! face/edge/vertex modification flags
    //! ------------------------------------
    MapOfIsFaceModified.clear();
    MapOfIsEdgeModified.clear();
    MapOfIsVertexModified.clear();

    //! ---------------------------------
    //! healing/defeaturing and defaults
    //! ---------------------------------
    mapOfDiscretizer.clear();
    mapOfAngularDeflection.clear();
    mapOfLinearDeflection.clear();
    MapOfIsBodyDefeaturingOn.clear();
    MapOfIsBodyHealingOn.clear();
    MapOfIsMeshSimplificationOn.clear();
    MapOfMeshSimplificationBy.clear();
    MapOfBodyDefeaturingParameterValue.clear();
    defaultIsBodyDefeaturingOn.clear();
    defaultIsBodyHealingOn.clear();
    defaultIsMeshSimplificationOn.clear();
    defaultMeshDefeaturingParameterValue.clear();
    defaultMeshDefeaturingBy.clear();

    //! ----------------------------------
    //! meshing process in memory/on disk
    //! ----------------------------------
    MapOfIsMeshingRunningInMemory.clear();
    defaultIsMeshingRunningInMemory.clear();

    //! ----------------------
    //! The mesh data sources
    //! ----------------------
    ArrayOfMeshDS.clear();
    ArrayOfMeshDS2D.clear();
    ArrayOfMeshDSOnFaces.clear();
    ArrayOfMeshDSOnEdges.clear();
    //! -----------------
    //!  mesh parameters
    //! -----------------
    mapOfUseBRep.clear();
    ArrayOfMesh2DEngine.clear();
    ArrayOfMesh3DEngine.clear();
    ArrayOfMeshOrder.clear();
    ArrayOfGradingValue.clear();
    ArrayOfMinBodyElementSize.clear();
    ArrayOfMaxBodyElementSize.clear();
    mapOfIsElementStraight.clear();
    ArrayOfSmoothingSteps.clear();

    //! -------------------------------
    //! the default meshing parameters
    //! -------------------------------
    defaultMapOfUseBRep.clear();
    default2DMeshEngine.clear();
    default3DMeshEngine.clear();
    defaultMeshOrder.clear();
    defaultVolumeGrading.clear();
    defaultVolumeMinElementSize.clear();
    defaultVolumeMaxElementSize.clear();
    defaultIsElementStraight.clear();
    defaultSmoothingSteps.clear();
    defaultMapOfDiscretizer.clear();
    defaultMapOfAngularDeflection.clear();
    defaultMapOfLinearDeflection.clear();

    //! -------------------------------------
    //! a flag on the mesh: must be updated?
    //! -------------------------------------
    ArrayOfMeshIsToBeUdpdated.clear();

    //! ------------------------------------------------------------
    //! the type of mesh for each body - surface mesh & volume mesh
    //! ------------------------------------------------------------
    ArrayOfDefaultMeshType.clear();
}
