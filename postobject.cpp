//! ----------------
//! custom includes
//! ----------------
#include "postobject.h"
#include "graphicstools.h"
#include "meshtools.h"
#include <meshdatabase.h>
#include "ccout.h"
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_deformeddatasource2d.h>

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include <TopAbs_ShapeEnum.hxx>

//! ---
//! Qt
//! ---
#include <QMap>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! -------
//! system
//! -------
#include <Windows.h>

/*
___postObject();
___postObject(const QMap<int, QMap<int, double> > &resMap);
___postObject(const QMap<int,QMap<int,double>> &resMap, const QVector<GeometryTag> &aVecLoc);
___postObject(const QMap<int,QMap<int,double>> &resMap, const QVector<GeometryTag> &aVecLoc, const QString &aName);
*/

//! ----------------------
//! function: constructor
//! details:  default
//! ----------------------
postObject::postObject()
{
    //! empty
    theData = QMap<GeometryTag,QList<QMap<int,double>>>();
    vecLoc = QVector<GeometryTag>();
    name = QString();
    theMeshes = QMap<GeometryTag,opencascade::handle<MeshVS_Mesh>>();
    mySolutionDataComponent = 0;
    myShowSolidMeshAsSurface = true;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myShowSolidMeshAsSurface = true;
    myScale = 1.0;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
postObject::postObject(const QMap<GeometryTag, QList<QMap<int, double>>> &resMap)
{
    theData = resMap;

    //! empty
    vecLoc = QVector<GeometryTag>();
    name = QString();
    theMeshes = QMap<GeometryTag,opencascade::handle<MeshVS_Mesh>>();
    mySolutionDataComponent = 0;
    myShowSolidMeshAsSurface = true;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myShowSolidMeshAsSurface = true;
    myScale = 1.0;
}

//! --------------------------------
//! function: constructor
//! details:  unused for the moment
//! --------------------------------
postObject::postObject(const QMap<GeometryTag,QList<QMap<int,double>>> &resMap, const QVector<GeometryTag> &aVecLoc)
{
    theData = resMap;
    vecLoc = aVecLoc;

    //! empty
    name = QString();
    theMeshes = QMap<GeometryTag,opencascade::handle<MeshVS_Mesh>>();
    mySolutionDataComponent = 0;
    myShowSolidMeshAsSurface = true;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myShowSolidMeshAsSurface = true;
    myScale = 1.0;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
postObject::postObject(const QMap<GeometryTag,QList<QMap<int,double>>> &resMap,
                       const QVector<GeometryTag> &aVecLoc,
                       const QString& aName)
{
    theData = resMap;
    vecLoc = aVecLoc;
    name = aName;
    theMeshes = QMap<GeometryTag,opencascade::handle<MeshVS_Mesh>>();
    mySolutionDataComponent = 0;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myShowSolidMeshAsSurface = true;
    myScale = 1.0;
}

//! --------------------------------
//! function: constructor
//! details:  this is actually used
//! --------------------------------
postObject::postObject(const QMap<GeometryTag,QList<QMap<int,double>>> &resMap,
                       const QVector<GeometryTag> &aVecLoc,
                       const QMap<GeometryTag,QMap<int,gp_Vec>> mapOfMapOfNodalDiplacements,
                       const QString& aName,
                       bool aShowSolidMeshAsSurface)
{
    theData = resMap;
    vecLoc = aVecLoc;
    name = aName;
    myMapOfNodalDisplacements = mapOfMapOfNodalDiplacements;
    myScale = 1.0;  // initial scale "true scale"
    theMeshes = QMap<GeometryTag,opencascade::handle<MeshVS_Mesh>>();
    mySolutionDataComponent = 0;
    myShowSolidMeshAsSurface = aShowSolidMeshAsSurface;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myShowSolidMeshAsSurface = true;
    myScale = 1.0;
}

//! ----------------
//! function: write
//! details:
//! ----------------
void postObject::write(ofstream &file)
{
    cout<<"postObject::write()->____write post object____"<<endl;

    //! ---------------
    //! write the name
    //! ---------------
    QList<QString> nameLines = name.split("\n");
    int NbNameLines = nameLines.length();
    //cout<<"postObject::write()->____number of lines: "<<NbNameLines<<"____"<<endl;

    //! ------------------------------------------
    //! write the name: write the number of lines
    //! ------------------------------------------
    file<<NbNameLines<<endl;

    //! ---------------------------------------
    //! write the lines of the name one by one
    //! ---------------------------------------
    for(int i=0; i<NbNameLines; i++)
    {
        //cout<<"postObject::write()->____line "<<i+1<<"th: "<<nameLines.at(i).toStdString()<<"____"<<endl;
        file<<nameLines.at(i).toStdString()<<endl;
    }

    //! -----------------------------------------
    //! write the scope: the number of locations
    //! -----------------------------------------
    file<<vecLoc.size()<<endl;
    //cout<<"postObject::write()->____number of locations: "<<vecLoc.size()<<"____"<<endl;

    //! -----------------------------------------
    //! write the scope: write the geometry tags
    //! -----------------------------------------
    for(QVector<GeometryTag>::iterator itLoc = vecLoc.begin();itLoc!=vecLoc.end(); ++itLoc)
    {
        const GeometryTag &aLoc = *itLoc;
        file<<aLoc.isParent<<endl;
        file<<aLoc.parentShapeNr<<endl;
        file<<aLoc.subTopNr<<endl;
        file<<aLoc.subShapeType<<endl;
    }

    //! --------------------------------------------------
    //! write the data: write the number of data blocks
    //! QMap<GeometryTag,QList<QMap<int,double>>> &resMap
    //! --------------------------------------------------
    file<<theData.keys().length()<<endl;
    //cout<<"postObject::write()->____number of keys in data: "<<theData.keys().length()<<"____"<<endl;

    for(QMap<GeometryTag,QList<QMap<int,double>>>::iterator it=theData.begin(); it!=theData.end(); ++it)
    {
        const GeometryTag &tag = it.key();
        const QList<QMap<int,double>> &l = it.value();

        //! --------------------------------------------------------
        //! write the number of components the result is defined by
        //! --------------------------------------------------------
        file<<l.size()<<endl;
        //cout<<"postObject::write()->____number of components: "<<l.size()<<"____"<<endl;

        //! --------------
        //! write the tag
        //! --------------
        file<<tag.isParent<<endl;
        file<<tag.parentShapeNr<<endl;
        file<<tag.subTopNr<<endl;
        file<<tag.subShapeType<<endl;

        for(int k=0; k<l.size(); k++)
        {
            const QMap<int,double> &res = l.at(k);

            //! ----------------------------------------------------
            //! write the size of the map for the current component
            //! ----------------------------------------------------
            int L = res.size();
            file<<L<<endl;
            //cout<<"postObject::write()->____writing the size of the result map____"<<L<<"____"<<endl;

            //! ---------------
            //! write the data
            //! ---------------
            for(QMap<int,double>::const_iterator itRes = res.begin(); itRes != res.cend(); ++itRes)
            {
                int nodeID = itRes.key();
                double val = itRes.value();
                file<<nodeID<<"\t"<<val<<endl;
            }
        }
    }

    //! -------------------------------------------
    //! write the mesh: write the number of meshes
    //! -------------------------------------------
    int NbMeshes = theMeshes.size();
    //cout<<"postObject::write()->____number of mesh objects: "<<NbMeshes<<"____"<<endl;

    for(QMap<GeometryTag,occHandle(MeshVS_Mesh)>::iterator it = theMeshes.begin(); it != theMeshes.end(); ++it)
    {
        const occHandle(MeshVS_Mesh) aMeshVS = it.value();
        const occHandle(MeshVS_DataSource) &aMeshDS = aMeshVS->GetDataSource();
        if(aMeshDS.IsNull())
        {
            //cout<<"postObject::write()->____cannot write the mesh____"<<endl;
            NbMeshes=0;
            break;
        }
    }

    file<<NbMeshes<<endl;

    for(QMap<GeometryTag,occHandle(MeshVS_Mesh)>::iterator it = theMeshes.begin(); it != theMeshes.end(); ++it)
    {
        const GeometryTag &aTag = it.key();
        const occHandle(MeshVS_Mesh) aMeshVS = it.value();

        //! --------------
        //! write the tag
        //! --------------
        file<<aTag.isParent<<endl;
        file<<aTag.parentShapeNr<<endl;
        file<<aTag.subTopNr<<endl;
        file<<aTag.subShapeType<<endl;

        //! ---------------------------
        //! write the mesh data source
        //! ---------------------------
        const occHandle(MeshVS_DataSource) &aMeshDS = aMeshVS->GetDataSource();

        switch(aTag.subShapeType)
        {
        case TopAbs_SOLID:
        {
            //cout<<"/---------------------------/"<<endl;
            //cout<<"/   writing a volume mesh   /"<<endl;
            //cout<<"/---------------------------/"<<endl;
            this->writeIntoStream(file,aMeshDS);
        }
            break;

        case TopAbs_FACE:
        {
            //cout<<"/---------------------------/"<<endl;
            //cout<<"/   writing a face mesh     /"<<endl;
            //cout<<"/---------------------------/"<<endl;
            this->writeIntoStream(file,aMeshDS);
        }
            break;

        case TopAbs_EDGE:
        {
            //! to be implemented
        }
            break;

        case TopAbs_VERTEX:
        {
            //! to be implemented
        }
            break;
        }
    }
}

//! --------------------------------
//! function: constructor from file
//! details:
//! --------------------------------
postObject::postObject(ifstream &file)
{
    cout<<"postObject::read()->____start reading post object____"<<endl;

    //! --------------
    //! read the name
    //! --------------
    int NbNameLines;
    file>>NbNameLines;
    cout<<"postObject::read()->____reading number of lines: "<<NbNameLines<<"____"<<endl;

    for(int i=0; i<NbNameLines; i++)
    {
        std::string lineName;
        std::getline(file,lineName);
        cout<<"postObject::read()->____"<<lineName<<"____"<<endl;
        name.append(QString::fromStdString(lineName)+"\n");
    }

    //! ----------------------------------------
    //! read the scope: the number of locations
    //! ----------------------------------------
    int NbLocs;
    file>>NbLocs;
    cout<<"postObject::read()->____number of locations: "<<NbLocs<<"____"<<endl;

    //! ---------------------------------------
    //! read the scope: read the geometry tags
    //! ---------------------------------------
    for(int i=0; i<NbLocs; i++)
    {
        //! --------------------
        //! build a GeometryTag
        //! --------------------
        GeometryTag loc;

        int isParentInt;
        bool isParent;
        file>>isParentInt;
        if(isParentInt) isParent = true; else isParent = false;
        loc.isParent = isParent;

        file>>loc.parentShapeNr;
        file>>loc.subTopNr;
        int shapeTypeInt;
        file>>shapeTypeInt;
        TopAbs_ShapeEnum type = static_cast<TopAbs_ShapeEnum>(shapeTypeInt);
        loc.subShapeType = type;

        //! -------
        //! append
        //! -------
        vecLoc.push_back(loc);
    }

    //! --------------------------------------------------
    //! read the data: read the number of data blocks
    //! QMap<GeometryTag,QList<QMap<int,double>>> &resMap
    //! --------------------------------------------------
    int NbBlocks;
    file>>NbBlocks;
    //cout<<"postObject::read()->____number of outer keys: "<<NbBlocks<<"____"<<endl;

    for(int i=0; i<NbBlocks; i++)
    {
        //! -------------------------------------------------------
        //! read the number of components the result is defined by
        //! -------------------------------------------------------
        int NbComponents;
        file>>NbComponents;
        //cout<<"postObject::read()->____number of components within the data file: "<<NbComponents<<"____"<<endl;

        //! -------------
        //! read the tag
        //! -------------
        GeometryTag aTag;
        int val;
        file>>val;
        if(val==0) aTag.isParent = false; else aTag.isParent = true;
        file>>aTag.parentShapeNr;
        file>>aTag.subTopNr;
        int type;
        file>>type;
        TopAbs_ShapeEnum shapeType = static_cast<TopAbs_ShapeEnum>(type);
        aTag.subShapeType = shapeType;

        QList<QMap<int,double>> listRes;
        for(int k=0; k<NbComponents; k++)
        {
            //! ---------------------------------------------------
            //! read the size of the map for the current component
            //! ---------------------------------------------------
            int L;
            file>>L;

            //! --------------
            //! read the data
            //! --------------
            int nodeID;
            double val;
            QMap<int,double> res;
            for(int j=0; j<L; j++)
            {
                file>>nodeID>>val;
                res.insert(nodeID,val);
            }
            listRes<<res;
        }
        theData.insert(aTag,listRes);
    }

    //! -----------------------------------------
    //! read the mesh: read the number of meshes
    //! -----------------------------------------
    int NbMeshes;
    file>>NbMeshes;
    //cout<<"postObject::read()->____number of meshes "<<NbMeshes<<"____"<<endl;

    //! --------------
    //! read the mesh
    //! --------------
    QMap<GeometryTag, occHandle(MeshVS_DataSource)> aMapOfMeshDataSources;
    for(int i=0; i<NbMeshes; i++)
    {
        //! -------------
        //! read the tag
        //! -------------
        GeometryTag aTag;
        int val;
        file>>val;
        if(val==1) aTag.isParent = true; else aTag.isParent = false;
        file>>aTag.parentShapeNr;
        file>>aTag.subTopNr;
        int type;
        file>>type;
        aTag.subShapeType = static_cast<TopAbs_ShapeEnum>(type);

        //! --------------------------
        //! read the mesh data source
        //! --------------------------
        occHandle(MeshVS_DataSource) aMeshDS;
        this->readMeshFromStream(file,aMeshDS);

        //! ----------------------------------
        //! fill the map of mesh data sources
        //! ----------------------------------
        aMapOfMeshDataSources.insert(aTag,aMeshDS);
    }

    //! -----------------------------------------------------
    //! create the colored meshes using the autoscale option
    //! -----------------------------------------------------
    this->buildMeshIO(aMapOfMeshDataSources,-9999,9999,10,true,0);
}

//! ----------------------
//! funtion: resetMeshes
//! details: reset meshes
//! ----------------------
void postObject::resetMeshes()
{
    this->theMeshes.clear();
}

//! ---------------------------------------------------
//! function: buildMeshIO
//! details:  build the colored mesh and the color box
//! ---------------------------------------------------
void postObject::buildMeshIO(const mapOfMeshDataSources &aMapOfMeshDataSources,
                             double min, double max,
                             int Nlevels,
                             bool autoscale,
                             int component,
                             bool showMeshEdges)
{
    cout<<"postObject::buildMeshIO()->____function called____"<<endl;

    //! -----------------------------------
    //! this parameter should be connected
    //! -----------------------------------
    bool showSolidMeshAsSurface = false;

    //! -----------------
    //! the data content
    //! -----------------
    double theMin, theMax;
    QMap<GeometryTag,QList<QMap<int,double>>>::const_iterator it = theData.cbegin();
    for(mapOfMeshDataSources::const_iterator anIt = aMapOfMeshDataSources.begin(); anIt!= aMapOfMeshDataSources.cend() && it!= theData.cend(); ++anIt, ++it)
    {
        const GeometryTag &loc = anIt.key();
        occHandle(MeshVS_DataSource) curMeshDS;

        //! -------------------------------------------------------
        //! if the current shape on which results are defined
        //! is a solid, it means that a volume mesh for that solid
        //! is defined: generate the interactive object according
        //! to the value of the flag "showSolidMeshAsSurface"
        //! -------------------------------------------------------
        if(loc.subShapeType == TopAbs_SOLID && showSolidMeshAsSurface == true)
        {
            const occHandle(Ng_MeshVS_DataSource3D) &volumeMeshDS = occHandle(Ng_MeshVS_DataSource3D)::DownCast(anIt.value());

            //! ------------------------------------------------------------------
            //! generate the surface mesh topologically - to be changed
            //! (I mean access the mesh data base directly... to do, if possible,
            //! since this object does not have access to the mesh data base)
            //! ------------------------------------------------------------------
            if(volumeMeshDS->myFaceToElements.size()==0) volumeMeshDS->buildFaceToElementConnectivity();
            curMeshDS = new Ng_MeshVS_DataSource2D(volumeMeshDS);
        }
        else
        {
            curMeshDS = anIt.value();
        }
        const QList<QMap<int, double>> &listOfRes = theData.value(loc);
        const QMap<int, double> &res = listOfRes.at(component);

        //! -----------------------------
        //! min and max for the colorbox
        //! -----------------------------
        if(autoscale)
        {
            std::pair<double,double> minmax = this->getMinMax(component);
            theMin = minmax.first;
            theMax = minmax.second;
        }
        else
        {
            theMin = min;
            theMax = max;
        }

        //! ----------------------------
        //! the new mesh (colored) mesh
        //! ----------------------------
        occHandle(MeshVS_Mesh) aColoredMesh;
        //MeshTools::buildColoredMesh(curMeshDS,res,theMin,theMax,Nlevels,aColoredMesh,false);

        QMap<int,gp_Vec> displacementMap = myMapOfNodalDisplacements.value(loc);
        MeshTools::buildDeformedColoredMesh(curMeshDS,res,displacementMap,1.0,theMin,theMax,Nlevels,aColoredMesh,showMeshEdges);
        theMeshes.insert(loc,aColoredMesh);
    }

    //! ---------------------------------
    //! create the color box with labels
    //! ---------------------------------
    //cout<<"____"<<theMin<<"____"<<endl;
    //cout<<"____"<<theMax<<"____"<<endl;

    graphicsTools::createColorBox(theMin, theMax, Nlevels, AISColorScale);
    TCollection_ExtendedString title(name.toStdString().c_str());
    AISColorScale->SetTitle(title);

    //! ------------------------
    //! set the private members
    //! ------------------------
    mySolutionDataComponent = component;
    myIsAutoscale = autoscale;
    myMin = theMin;
    myMax = theMax;
    myNbLevels = Nlevels;
    myShowSolidMeshAsSurface = showSolidMeshAsSurface;
}

//! -------------------------------------------
//! function: getMinMax
//! details:  used for autoscaling data values
//! -------------------------------------------
std::pair<double,double> postObject::getMinMax(int component)
{
    cout<<"postObject::getMinMax()->____function called____"<<endl;

    //! ---------------------------------------------
    //! first element => min; second element => max;
    //! ---------------------------------------------
    std::pair<double,double> minmax;
    minmax.first = 1e80;
    minmax.second = -1e80;

    //! -------------------
    //! scan the locations
    //! -------------------
    if(theData.isEmpty())
    {
        cout<<"postObject::getMinMax()->___empty data____"<<endl;
    }
    QMap<GeometryTag,QList<QMap<int,double>>>::const_iterator rIt;
    for(rIt = theData.cbegin(); rIt!= theData.cend(); ++rIt)
    {
        QList<QMap<int,double>> lres = rIt.value();
        QMap<int,double> res = lres.at(component);          //! component
        QMap<int,double>::iterator rrIt;

        if(res.isEmpty())
        {
            cout<<"postObject::getMinMax()->___empty map data____"<<endl;
        }
        for(rrIt = res.begin(); rrIt!=res.end(); ++rrIt)
        {
            double aVal = rrIt.value();
            if(aVal<=minmax.first) minmax.first = aVal;
            if(aVal>=minmax.second) minmax.second = aVal;
        }
    }
    cout<<"postObject::getMinMax()->____Min = "<<minmax.first<<"____Max = "<<minmax.second<<"____"<<endl;
    return minmax;
}

//! ---------------------------------------------------
//! function: update
//! details:  retrieve the mesh from the database
//!           internally calls postObject::buildMeshIO
//! ---------------------------------------------------
void postObject::update(meshDataBase *mDB, int component)
{
    cout<<"postObject::update()->____update has been requested____"<<endl;

    if(theData.isEmpty())
    {
        cout<<"postObject::update()->____cannot generate the mesh object: no data____"<<endl;
        return;
    }

    int topNr;
    int parentShapeNr;
    occHandle(MeshVS_DataSource) meshDS;
    mapOfMeshDataSources aMapOfMeshDS;

    for(QVector<GeometryTag>::iterator it = vecLoc.begin(); it!= vecLoc.end(); ++it)
    {
        const GeometryTag &loc = *it;
        TopAbs_ShapeEnum type = static_cast<TopAbs_ShapeEnum>(loc.subShapeType);
        switch(type)
        {
        case TopAbs_SOLID:
        {
            cout<<"postObject::update()->____handling the \"SOLID\" case____"<<endl;
            topNr = loc.parentShapeNr;
            aMapOfMeshDS.insert(loc, mDB->ArrayOfMeshDS.value(topNr));
        }
            break;

        case TopAbs_FACE:
        {
            cout<<"postObject::update()->____handling the \"FACE\" case____"<<endl;
            parentShapeNr = loc.parentShapeNr;
            topNr = loc.subTopNr;
            occHandle(MeshVS_DataSource) theMeshToInsert = mDB->ArrayOfMeshDSOnFaces.getValue(parentShapeNr,topNr);
            aMapOfMeshDS.insert(loc, theMeshToInsert);
        }
            break;

        case TopAbs_EDGE:
        {
            cout<<"postObject::update()->____handling the \"EDGE\" case____"<<endl;
            topNr = loc.subTopNr;
            parentShapeNr = loc.parentShapeNr;
            aMapOfMeshDS.insert(loc, mDB->ArrayOfMeshDSOnEdges.getValue(parentShapeNr,topNr));
        }
            break;

        case TopAbs_VERTEX:
            cout<<"postObject::update()->____handling the \"VERTEX\" case____"<<endl;
            //! to be implemented. To do ...
            break;
        }
    }

    //! ------------------------------------------------
    //! create the post object using "Autoscale" option
    //! ------------------------------------------------
    bool autoscale = true;
    this->buildMeshIO(aMapOfMeshDS,0,0,INITIAL_NUMBER_OF_COLORBOX_LEVELS,autoscale,component);
}

//! ---------------------
//! function: updateView
//! details:
//! ---------------------
void postObject::updateView(bool showMeshEdges)
{
    cout<<"postObject::updateView()->____show elements: "<<(showMeshEdges? "yes, show elements":"no, hide elements")<<"____"<<endl;
    for(QMap<GeometryTag,occHandle(MeshVS_Mesh)>::iterator it = theMeshes.begin(); it!=theMeshes.end(); ++it)
    {
        occHandle(MeshVS_Mesh) curMeshVS_Mesh = *it;
        curMeshVS_Mesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges,showMeshEdges);
    }
}

//! -------------------------
//! function: update mapping
//! details:
//! -------------------------
void postObject::updateMapping(int mapping)
{
    cout<<"postObject::updateView()->____function called____"<<endl;
    for(QMap<GeometryTag,occHandle(MeshVS_Mesh)>::iterator it = theMeshes.begin(); it!=theMeshes.end(); ++it)
    {
        occHandle(MeshVS_Mesh) aMeshObject = *it;

    }
}

//! ---------------------------
//! function: updateScaledView
//! details:
//! ---------------------------
void postObject::updateScaledView()
{
    if(theMeshes.isEmpty()) return;
    for(QMap<GeometryTag,occHandle(MeshVS_Mesh)>::iterator it = theMeshes.begin(); it!= theMeshes.end(); ++it)
    {
        const occHandle(MeshVS_Mesh) &curMeshVS = it.value();
        const occHandle(Ng_MeshVS_DeformedDataSource2D)  &curMeshDeformedDS = occHandle(Ng_MeshVS_DeformedDataSource2D)::DownCast(curMeshVS->GetDataSource());
        if(myScale == 0.0)
        {
            cerr<<"postObject::updateScaledView()->____updating scale: undeformed____"<<endl;
            occHandle(Ng_MeshVS_DataSource2D) nonDeformedDS = occHandle(Ng_MeshVS_DataSource2D)::DownCast(curMeshDeformedDS->GetNonDeformedDataSource());
            curMeshVS->SetDataSource(nonDeformedDS);
        }
        else
        {
            cerr<<"postObject::updateScaledView()->____updating scale: "<<myScale<<"____"<<endl;
            curMeshDeformedDS->SetMagnify(myScale);
        }
    }
}


//! ---------------------------
//! function: writeMesh
//! details:  use a C++ stream
//! ---------------------------
void postObject::writeIntoStream(ofstream &os, const occHandle(MeshVS_DataSource) &aMeshDS)
{
    if(!os.is_open()) return;

    cout<<"postObject::writeIntoStream()->____function called____"<<endl;

    //! ----------------
    //! nodes, elements
    //! ----------------
    const TColStd_PackedMapOfInteger &aNodes = aMeshDS->GetAllNodes();
    const TColStd_PackedMapOfInteger &mapOfElements = aMeshDS->GetAllElements();

    if(aNodes.Extent()==0)
    {
        cout<<"Ng_MeshVS_DataSource3D::writeIntoStream()->____no node to write____"<<endl;
        return;
    }
    if(mapOfElements.Extent()==0)
    {
        cout<<"Ng_MeshVS_DataSource3D::writeIntoStream()->____no element to write____"<<endl;
        return;
    }

    //! ---------------
    //! mesh dimension
    //! ---------------
    os<<int(3)<<endl;
    os<<aNodes.Extent()<<endl;

    Standard_Real aCoordsBuf[3];
    TColStd_Array1OfReal aCoords(*aCoordsBuf,1,3);
    Standard_Integer nbNodes;
    MeshVS_EntityType aType;
    TColStd_MapIteratorOfPackedMapOfInteger anIter;

    //! ------------------------------------
    //! write the nodID and the coordinates
    //! ------------------------------------
    int k=1;
    for(anIter.Initialize(aNodes);anIter.More();anIter.Next())
    {
        int aKey = anIter.Key();
        if(aMeshDS->GetGeom(k,Standard_False,aCoords,nbNodes,aType)==true)
        {
            os<<aKey<<"\t"<<aCoords(1)<<"\t"<<aCoords(2)<<"\t"<<aCoords(3)<<endl;
            k++;
        }
    }

    //! -----------------------------
    //! write the number of elements
    //! -----------------------------
    int NbElements = mapOfElements.Extent();
    os<<NbElements<<endl;

    //! ------------------------------
    //! write the elements definition
    //! ------------------------------
    for(anIter.Initialize(mapOfElements);anIter.More();anIter.Next())
    {
        int aKey = anIter.Key();
        int NbNodes;
        int aNodesBuf[20];
        TColStd_Array1OfInteger nodeIDs(*aNodesBuf,1,20);
        aMeshDS->GetNodesByElement(aKey,nodeIDs,NbNodes);

        MeshVS_EntityType type;
        aMeshDS->GetGeomType(aKey,true,type);

        char header[16];
        if(type==MeshVS_ET_Face)
        {
            switch(NbNodes)
            {
            case 3: sprintf(header,"TRIG"); break;
            case 4: sprintf(header,"QUAD"); break;
            case 6: sprintf(header,"TRIG6"); break;
            case 8: sprintf(header,"QUAD8"); break;
            default: break;
            }
        }
        else if(type==MeshVS_ET_Volume)
        {
            switch(NbNodes)
            {
            case 4: sprintf(header,"TET"); break;
            case 8: sprintf(header,"HEXA"); break;
            case 10: sprintf(header,"TET10"); break;
            case 20: sprintf(header,"HEXA20"); break;
            case 6: sprintf(header,"PRISM"); break;
            case 5: sprintf(header,"PYRAM"); break;
            default: break;
            }
        }

        //! ----------------------------------------
        //! write the element key, the element type
        //! ----------------------------------------
        os<<aKey<<endl;
        os<<header<<endl;

        for(int i=1; i<NbNodes; i++) os<<nodeIDs(i)<<"\t";
        os<<nodeIDs(NbNodes)<<endl;
    }
}

//! -----------------------------
//! function: readMeshFromStream
//! details:
//! -----------------------------
#include <mesh.h>
bool postObject::readMeshFromStream(ifstream &stream, occHandle(MeshVS_DataSource) &aMeshDS)
{
    cout<<"postObject::readMeshFromStream()->____function called____"<<endl;

    std::vector<mesh::meshPoint> nodes;
    std::vector<mesh::meshElement> elements;

    if(!stream.is_open()) return false;

    //! -------------------------------------------
    //! read the dimension and the number of nodes
    //! -------------------------------------------
    int dim, NbNodes;
    stream>>dim;
    stream>>NbNodes;

    cout<<"postObject::readMeshFromStream()->____number of nodes: "<<NbNodes<<"____"<<endl;

    if(NbNodes==0) return false;
    if(dim == 3 && NbNodes<4) return false;
    if(dim ==2 && NbNodes<3) return false;

    //! ---------------------------
    //! read the nodes coordinates
    //! ---------------------------
    mesh::meshPoint aNode;
    for(int i=0; i<NbNodes;i++)
    {
        stream>>aNode.ID>>aNode.x>>aNode.y>>aNode.z;
        nodes.push_back(aNode);
        //cout<<"postObject::readMeshFromStream()->____reading node ID: "<<aNode.ID<<"("<<
        //      aNode.x<<", "<<aNode.y<<", "<<aNode.z<<")____"<<endl;
    }

    //! ----------------------------
    //! read the number of elements
    //! ----------------------------
    int NbElements;
    stream>>NbElements;
    cout<<"postObject::readMeshFromStream()->____number of elements: "<<NbElements<<"____"<<endl;
    if(NbElements==0) return false;

    //! -----------------------------
    //! read the element definitions
    //! -----------------------------
    for(int i=0; i<NbElements; i++)
    {
        mesh::meshElement anElement;
        char eType[16];

        stream>>anElement.ID;
        stream>>eType;

        if(strcmp(eType,"TRIG")==0)
        {
            anElement.type= TRIG;
            int V1,V2,V3;
            stream>>V1>>V2>>V3;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
        }
        if(strcmp(eType,"TRIG6")==0)
        {
            anElement.type= TRIG6;
            int V1,V2,V3,V4,V5,V6;
            stream>>V1>>V2>>V3>>V4>>V5>>V6;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
            anElement.theNodeIDs.push_back(V5);
            anElement.theNodeIDs.push_back(V6);
        }
        else if(strcmp(eType,"TET")==0)
        {
            anElement.type= TET;
            int V1,V2,V3,V4;
            stream>>V1>>V2>>V3>>V4;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
        }
        else if(strcmp(eType,"TET10")==0)
        {
            anElement.type= TET10;
            int V1,V2,V3,V4,V5,V6,V7,V8,V9,V10;
            stream>>V1>>V2>>V3>>V4>>V5>>V6>>V7>>V8>>V9>>V10;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
            anElement.theNodeIDs.push_back(V5);
            anElement.theNodeIDs.push_back(V6);
            anElement.theNodeIDs.push_back(V7);
            anElement.theNodeIDs.push_back(V8);
            anElement.theNodeIDs.push_back(V9);
            anElement.theNodeIDs.push_back(V10);
        }
        else if(strcmp(eType,"PYRAM")==0)
        {
            anElement.type= PYRAM;
            int V1,V2,V3,V4,V5;
            stream>>V1>>V2>>V3>>V4>>V5;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
            anElement.theNodeIDs.push_back(V5);
        }
        else if(strcmp(eType,"PRISM")==0)
        {
            anElement.type= PRISM;
            int V1,V2,V3,V4,V5,V6;
            stream>>V1>>V2>>V3>>V4>>V5>>V6;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
            anElement.theNodeIDs.push_back(V5);
            anElement.theNodeIDs.push_back(V6);
        }
        elements.push_back(anElement);
    }

    switch(dim)
    {
    case 2:
    {
        aMeshDS = new Ng_MeshVS_DataSourceFace(elements, nodes);
        if(aMeshDS.IsNull()) return false;
    }
        break;
    case 3:
    {
        aMeshDS = new Ng_MeshVS_DataSource3D(elements, nodes);
        if(aMeshDS.IsNull()) return false;
    }
        break;
    }
    return true;
}
