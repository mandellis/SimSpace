//! ----------------
//! custom includes
//! ----------------
#include "postobject.h"
#include "graphicstools.h"
#include "meshtools.h"
#include <mesh.h>
#include <meshdatabase.h>
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_deformeddatasource2d.h>
#include <isostripbuilder.h>
#include <isostrip.h>
#include <meshslicer.h>

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <MeshVS_MeshPrsBuilder.hxx>
#include <MeshVS_ElementalColorPrsBuilder.hxx>

//! ---
//! Qt
//! ---
#include <QMap>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:  default
//! ----------------------
postObject::postObject()
{
    theData = std::map<GeometryTag,std::vector<std::map<int,double>>>();

    myVecLoc = std::vector<GeometryTag>();
    name = QString();

    theMeshes = std::map<GeometryTag,occHandle(MeshVS_Mesh)>();
    theMeshDataSources = std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>();

    mySolutionDataComponent = 0;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myUseSurfaceMeshForVolumeResults = true;
    myScale = 1.0;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
postObject::postObject(const std::map<GeometryTag, std::vector<std::map<int, double>>> &resMap)
{
    theData = resMap;

    myVecLoc = std::vector<GeometryTag>();
    name = QString();
    theMeshes = std::map<GeometryTag,occHandle(MeshVS_Mesh)>();
    theMeshDataSources = std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>();

    mySolutionDataComponent = 0;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myUseSurfaceMeshForVolumeResults = true;
    myScale = 1.0;
}

//! --------------------------------
//! function: constructor
//! details:  unused for the moment
//! --------------------------------
postObject::postObject(const std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap, const std::vector<GeometryTag> &aVecLoc)
{
    theData = resMap;
    myVecLoc = aVecLoc;

    name = QString();
    theMeshes = std::map<GeometryTag,occHandle(MeshVS_Mesh)>();
    theMeshDataSources = std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>();
    theMeshDataSourcesForView = std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>();
    mySolutionDataComponent = 0;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myUseSurfaceMeshForVolumeResults = true;
    myScale = 1.0;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
postObject::postObject(const std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap,
                       const std::vector<GeometryTag> &aVecLoc,
                       const QString& aName)
{
    theData = resMap;
    myVecLoc = aVecLoc;
    name = aName;
    theMeshes = std::map<GeometryTag,occHandle(MeshVS_Mesh)>();
    theMeshDataSources = std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>();
    theMeshDataSourcesForView = std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>();
    mySolutionDataComponent = 0;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
    myUseSurfaceMeshForVolumeResults = true;
    myScale = 1.0;
}

//! --------------------------------
//! function: constructor
//! details:  this is actually used
//! --------------------------------
postObject::postObject(const std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap,
                       const std::vector<GeometryTag> &aVecLoc,
                       const std::map<GeometryTag,std::map<int,gp_Vec>> mapOfMapOfNodalDiplacements,
                       const QString& aName,
                       bool useSurfaceMeshForVolumeResults)

{
    theData = resMap;
    myVecLoc = aVecLoc;
    name = aName;
    myMapOfNodalDisplacements = mapOfMapOfNodalDiplacements;
    myScale = 1.0;  // initial scale "true scale"
    theMeshes = std::map<GeometryTag,occHandle(MeshVS_Mesh)>();
    theMeshDataSources = std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>();
    theMeshDataSourcesForView = std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>();
    mySolutionDataComponent = 0;
    myUseSurfaceMeshForVolumeResults = useSurfaceMeshForVolumeResults;
    myMax = 100;
    myMin = -100;
    myIsAutoscale = true;
    myNbLevels = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
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

    //! ------------------------------------------
    //! write the name: write the number of lines
    //! ------------------------------------------
    file<<NbNameLines<<endl;

    //! ---------------------------------------
    //! write the lines of the name one by one
    //! ---------------------------------------
    for(int i=0; i<NbNameLines; i++) file<<nameLines.at(i).toStdString()<<endl;

    //! -----------------------------------------
    //! write the scope: the number of locations
    //! -----------------------------------------
    file<<myVecLoc.size()<<endl;

    //! -----------------------------------------
    //! write the scope: write the geometry tags
    //! -----------------------------------------
    for(std::vector<GeometryTag>::iterator itLoc = myVecLoc.begin();itLoc!=myVecLoc.end(); ++itLoc)
    {
        const GeometryTag &aLoc = *itLoc;
        file<<aLoc.isParent<<endl;
        file<<aLoc.parentShapeNr<<endl;
        file<<aLoc.subTopNr<<endl;
        file<<aLoc.subShapeType<<endl;
    }

    //! ----------------------------------------------------------------
    //! write the data: write the number of data blocks
    //! std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap
    //! ----------------------------------------------------------------
    file<<theData.size()<<endl;

    for(std::map<GeometryTag,std::vector<std::map<int,double>>>::iterator it=theData.begin(); it!=theData.end(); ++it)
    {
        const GeometryTag &tag = it->first;
        const std::vector<std::map<int,double>> &l = it->second;

        //! --------------------------------------------------------
        //! write the number of components the result is defined by
        //! --------------------------------------------------------
        file<<l.size()<<endl;

        //! --------------
        //! write the tag
        //! --------------
        file<<tag.isParent<<endl;
        file<<tag.parentShapeNr<<endl;
        file<<tag.subTopNr<<endl;
        file<<tag.subShapeType<<endl;

        for(int k=0; k<l.size(); k++)
        {
            const std::map<int,double> &res = l.at(k);

            //! ----------------------------------------------------
            //! write the size of the map for the current component
            //! ----------------------------------------------------
            int L = (int)res.size();
            file<<L<<endl;

            //! ---------------
            //! write the data
            //! ---------------
            for(std::map<int,double>::const_iterator itRes = res.begin(); itRes != res.cend(); ++itRes)
            {
                int nodeID = itRes->first;
                double val = itRes->second;
                file<<nodeID<<"\t"<<val<<endl;
            }
        }
    }

    //! -------------------------------------------
    //! write the mesh: write the number of meshes
    //! -------------------------------------------
    int NbMeshes = (int)theMeshDataSources.size();
    for(std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>::iterator it = theMeshDataSources.begin(); it!=theMeshDataSources.end(); it++)
    {
        const occHandle(MeshVS_DeformedDataSource) &aMeshDS = it->second;
        if(aMeshDS.IsNull())
        {
            NbMeshes=0;
            break;
        }
    }
    file<<NbMeshes<<endl;

    for(std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>::iterator it = theMeshDataSources.begin(); it != theMeshDataSources.end(); ++it)
    {
        const GeometryTag &aTag = it->first;
        const occHandle(MeshVS_DeformedDataSource) &aDeformedMeshDS = it->second;

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
        const occHandle(MeshVS_DataSource) &aMeshDS = aDeformedMeshDS->GetNonDeformedDataSource();
        switch(aTag.subShapeType)
        {
        case TopAbs_SOLID: this->writeIntoStream(file,aMeshDS); break;
        case TopAbs_FACE: this->writeIntoStream(file,aMeshDS); break;
        case TopAbs_EDGE: /* to be implemented */ break;
        case TopAbs_VERTEX: /* to be implemented */ break;
        }
    }
    //! ----------------------------
    //! write generation paramaters
    //! ----------------------------
    file<<mySolutionDataComponent<<endl;
    file<<myMin<<endl;
    file<<myMax<<endl;
    file<<myNbLevels<<endl;
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

    for(int i=0; i<NbNameLines; i++)
    {
        std::string lineName;
        std::getline(file,lineName);
        name.append(QString::fromStdString(lineName)+"\n");
    }

    //! ----------------------------------------
    //! read the scope: the number of locations
    //! ----------------------------------------
    int NbLocs;
    file>>NbLocs;

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
        myVecLoc.push_back(loc);
    }

    //! ----------------------------------------------------------------
    //! read the data: read the number of data blocks
    //! std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap
    //! ----------------------------------------------------------------
    int NbBlocks;
    file>>NbBlocks;

    for(int i=0; i<NbBlocks; i++)
    {
        //! -------------------------------------------------------
        //! read the number of components the result is defined by
        //! -------------------------------------------------------
        int NbComponents;
        file>>NbComponents;

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

        std::vector<std::map<int,double>> listRes;
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
            std::map<int,double> res;
            for(int j=0; j<L; j++)
            {
                file>>nodeID>>val;
                res.insert(std::make_pair(nodeID,val));
            }
            listRes.push_back(res);
        }
        theData.insert(std::make_pair(aTag,listRes));
    }

    //! -----------------------------------------
    //! read the mesh: read the number of meshes
    //! -----------------------------------------
    int NbMeshes;
    file>>NbMeshes;

    //! --------------
    //! read the mesh
    //! --------------
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
        occHandle(MeshVS_DeformedDataSource) aDefMeshDS = new MeshVS_DeformedDataSource(aMeshDS,1.0);
        theMeshDataSources.insert(std::make_pair(aTag,aDefMeshDS));
    }
    //! -------------------------------
    //! read the generation parameters
    //! -------------------------------
    file>>mySolutionDataComponent;
    file>>myMin;
    file>>myMax;
    file>>myNbLevels;
}

//! ----------------------
//! funtion: resetMeshes
//! details: reset meshes
//! ----------------------
void postObject::resetMeshes()
{
    this->theMeshes.clear();
}

//! ----------------------
//! function: buildMeshIO
//! details:
//! ----------------------
void postObject::buildMeshIO(double min, double max, int Nlevels, bool autoscale, int component, double deformationScale)
{
    cout<<"postObject::buildMeshIO()->____function called____"<<endl;

    //! ------------------------
    //! set the private members
    //! ------------------------
    mySolutionDataComponent = component;
    myIsAutoscale = autoscale;
    myNbLevels = Nlevels;
    theMeshes.clear();

    //! -----------------
    //! the data content
    //! -----------------
    std::map<GeometryTag,std::vector<std::map<int,double>>>::const_iterator it = theData.cbegin();
    for(std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>::const_iterator anIt = theMeshDataSources.cbegin(); anIt!= theMeshDataSources.cend() && it!= theData.cend(); ++anIt, ++it)
    {
        const GeometryTag &loc = anIt->first;
        occHandle(MeshVS_DataSource) curMeshDS;
        if(loc.subShapeType == TopAbs_SOLID && myUseSurfaceMeshForVolumeResults == true)
        {
            cout<<"____USING THE EXTERIOR MESH____"<<endl;
            occHandle(Ng_MeshVS_DataSource3D) volumeMeshDS = new Ng_MeshVS_DataSource3D(anIt->second->GetNonDeformedDataSource());

            //! ----------------------------------------
            //! generate the surface mesh topologically
            //! ----------------------------------------
            if(volumeMeshDS->myFaceToElements.isEmpty()) volumeMeshDS->buildFaceToElementConnectivity();
            curMeshDS = new Ng_MeshVS_DataSource2D(volumeMeshDS);
        }
        else
        {
            curMeshDS = anIt->second->GetNonDeformedDataSource();
        }

        //! ----------------------------------
        //! build a deformed mesh data source
        //! ----------------------------------
        occHandle(MeshVS_DeformedDataSource) theDeformedDS = new MeshVS_DeformedDataSource(curMeshDS,deformationScale);
        const std::map<int,gp_Vec> &displacementMap = myMapOfNodalDisplacements.at(loc);
        theDeformedDS->SetNonDeformedDataSource(curMeshDS);

        for(TColStd_MapIteratorOfPackedMapOfInteger it(curMeshDS->GetAllNodes()); it.More(); it.Next())
        {
            int globalNodeID = it.Key();
            const gp_Vec &d = displacementMap.at(globalNodeID);
            theDeformedDS->SetVector(globalNodeID,d);
        }
        theDeformedDS->SetMagnify(deformationScale);
        theMeshDataSourcesForView[loc]=theDeformedDS;   //! abruptly replace - do not use "insert"

        const std::vector<std::map<int, double>> &listOfRes = theData.at(loc);
        const std::map<int, double> &res = listOfRes.at(component);

        //! -------------------------------------------
        //! min and max for the colorbox and isostrips
        //! -------------------------------------------
        if(autoscale)
        {
            std::pair<double,double> minmax = this->getMinMax(component);
            myMin = minmax.first;
            myMax = minmax.second;
        }
        else { myMin = min; myMax = max; }

        occHandle(MeshVS_Mesh) aColoredMesh;
        switch(Global::status().myResultPresentation.theTypeOfPresentation)
        {
        case resultPresentation::typeOfPresentation_isostrips:
            MeshTools::buildIsoStrip(theDeformedDS,res,myMin,myMax,myNbLevels,aColoredMesh);
            break;
        case resultPresentation::typeOfPresentation_nodalresults:
            MeshTools::buildColoredMesh(theDeformedDS,res,aColoredMesh,myMin,myMax,myNbLevels);
            break;
        case resultPresentation::typeOfPresentation_isosurfaces:
            MeshTools::buildIsoSurfaces(theDeformedDS,res,myMin,myMax,myNbLevels,aColoredMesh,false);
            break;
        case resultPresentation::typeOfPresentation_isolines:
            break;
        }

        theMeshes.insert(std::make_pair(loc,aColoredMesh));
    }

    //! ---------------------------------
    //! create the color box with labels
    //! ---------------------------------
    graphicsTools::createColorBox(myMin, myMax, myNbLevels, AISColorScale);
    TCollection_ExtendedString title(name.toStdString().c_str());
    AISColorScale->SetTitle(title);
    cout<<"postObject::buildMeshIO()->____exiting function____"<<endl;
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
    minmax.first = 1e20;
    minmax.second = -1e20;

    //! -------------------
    //! scan the locations
    //! -------------------
    if(theData.size()==0)
    {
        cout<<"postObject::getMinMax()->___empty data____"<<endl;
        return std::make_pair(0.0,0.0);
    }
    for(std::map<GeometryTag,std::vector<std::map<int,double>>>::const_iterator rIt = theData.cbegin(); rIt!= theData.cend(); ++rIt)
    {
        std::vector<std::map<int,double>> lres = rIt->second;
        std::map<int,double> res = lres.at(component);          //! component
        if(res.size()==0)
        {
            cout<<"postObject::getMinMax()->___empty map data____"<<endl;
        }
        for(std::map<int,double>::iterator rrIt = res.begin(); rrIt!=res.end(); ++rrIt)
        {
            double aVal = rrIt->second;
            if(aVal<=minmax.first) minmax.first = aVal;
            if(aVal>=minmax.second) minmax.second = aVal;
        }
    }
    //cout<<"postObject::getMinMax()->____min: "<<minmax.first<<", max: "<<minmax.second<<"____"<<endl;
    return minmax;
}

//! ----------------------------------------
//! function: init
//! details:  provide the mesh data sources
//! ----------------------------------------
void postObject::init(meshDataBase *mDB)
{
    cout<<"postObject::init()->____function called____"<<endl;

    if(theData.size()==0)
    {
        cout<<"postObject::update()->____cannot generate the mesh object: no data____"<<endl;
        return;
    }

    int topNr, parentShapeNr;
    occHandle(MeshVS_DataSource) meshDS;

    for(std::vector<GeometryTag>::iterator it = myVecLoc.begin(); it!= myVecLoc.end(); ++it)
    {
        const GeometryTag &loc = *it;
        TopAbs_ShapeEnum type = static_cast<TopAbs_ShapeEnum>(loc.subShapeType);
        switch(type)
        {
        case TopAbs_SOLID:
        {
            topNr = loc.parentShapeNr;
            const occHandle(MeshVS_DeformedDataSource) &aDeformedMeshDS = new MeshVS_DeformedDataSource(mDB->ArrayOfMeshDS.value(topNr),1.0);
            aDeformedMeshDS->SetNonDeformedDataSource(mDB->ArrayOfMeshDS.value(topNr));
            aDeformedMeshDS->SetMagnify(1.0);
            theMeshDataSources.insert(std::make_pair(loc, aDeformedMeshDS));
        }
            break;

        case TopAbs_FACE:
        {
            parentShapeNr = loc.parentShapeNr;
            topNr = loc.subTopNr;
            occHandle(MeshVS_DataSource) theMeshToInsert = mDB->ArrayOfMeshDSOnFaces.getValue(parentShapeNr,topNr);            
            const occHandle(MeshVS_DeformedDataSource) &aDeformedMeshDS = new MeshVS_DeformedDataSource(mDB->ArrayOfMeshDSOnFaces.getValue(parentShapeNr,topNr),1.0);
            theMeshDataSources.insert(std::make_pair(loc, aDeformedMeshDS));
        }
            break;

        case TopAbs_EDGE:
        {
            topNr = loc.subTopNr;
            parentShapeNr = loc.parentShapeNr;
            const occHandle(MeshVS_DeformedDataSource) &aDeformedMeshDS = new MeshVS_DeformedDataSource(mDB->ArrayOfMeshDSOnEdges.getValue(parentShapeNr,topNr),1.0);
            theMeshDataSources.insert(std::make_pair(loc, aDeformedMeshDS));
        }
            break;

        case TopAbs_VERTEX:
            //! to be implemented. To do ...
            break;
        }
    }
    //cout<<"postObject::init()->____exiting function____"<<endl;
}

//! ------------------------
//! function: updateMapping
//! details:
//! ------------------------
void postObject::updateMapping(int mapping)
{
    cout<<"postObject::updateView()->____function called____"<<endl;
    for(std::map<GeometryTag,occHandle(MeshVS_Mesh)>::iterator it = theMeshes.begin(); it!=theMeshes.end(); ++it)
    {
        occHandle(MeshVS_Mesh) aMeshObject = it->second;
    }
}

//! ---------------------------
//! function: updateScaledView
//! details:
//! ---------------------------
void postObject::updateScaledView()
{
    cout<<"postObject::updateScaledView()->____function called____"<<endl;
    if(theMeshDataSources.size()==0) return;
    buildMeshIO(myMin, myMax, myNbLevels, myIsAutoscale, mySolutionDataComponent, myScale);
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
    }

    //! ----------------------------
    //! read the number of elements
    //! ----------------------------
    int NbElements;
    stream>>NbElements;
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

//! --------------------------------
//! function: computeHiddenElements
//! details:
//! --------------------------------
void postObject::computeHiddenElements(const std::map<int,std::vector<double>> &mapOfClipPlanes)
{
    //std::map<GeometryTag,occHandle(TColStd_HPackedMapOfInteger)> hiddenElementsIDs;

    //! ---------------------------------------
    //! check the number of active clip planes
    //! ---------------------------------------
    int NbClipPlanes = (int)mapOfClipPlanes.size();
    if(NbClipPlanes==0)
    {
       TColStd_PackedMapOfInteger e;
       occHandle(TColStd_HPackedMapOfInteger) he = new TColStd_HPackedMapOfInteger;
       he->ChangeMap() = e;
       for(std::map<GeometryTag,occHandle(MeshVS_Mesh)>::iterator it = theMeshes.begin(); it!=theMeshes.end(); it++)
       {
           const GeometryTag &aTag = it->first;
           myMapOfHiddenElements.insert(std::make_pair(aTag,he));
       }
       return;
    }

    //! ------------
    //! mesh slicer
    //! ------------
    meshSlicer aMeshSlicer;

    //! ------------------------
    //! iterate over the meshes
    //! ------------------------
    for(std::map<GeometryTag,occHandle(MeshVS_Mesh)>::iterator it_ = theMeshes.begin(); it_ != theMeshes.end(); it_++)
    {
        const GeometryTag &aTag = it_->first;
        const occHandle(MeshVS_Mesh) &aMeshObject = it_->second;
        const occHandle(MeshVS_DataSource) &aMeshDS = aMeshObject->GetDataSource();
        aMeshSlicer.setMeshDataSource(aMeshDS);

        TColStd_PackedMapOfInteger hiddenElementsIDsForTag;

        //! ------------------------
        //! iterate over the planes
        //! ------------------------
        for(std::map<int,std::vector<double>>::const_iterator it = mapOfClipPlanes.cbegin(); it != mapOfClipPlanes.cend(); it++)
        {
            //int clipPlaneID = it->first;
            const std::vector<double> &coeffs = it->second;
            occHandle(TColStd_HPackedMapOfInteger) HMapOfIDs;
            bool isDone = aMeshSlicer.perform(coeffs[0],coeffs[1],coeffs[2],coeffs[3],HMapOfIDs);
            if(isDone==false) continue;
            const TColStd_PackedMapOfInteger &mapOfIDs = HMapOfIDs->Map();
            hiddenElementsIDsForTag.Unite(mapOfIDs);
        }

        occHandle(TColStd_HPackedMapOfInteger) HHiddenElementsIDsForTag = new TColStd_HPackedMapOfInteger;
        HHiddenElementsIDsForTag->ChangeMap() = hiddenElementsIDsForTag;
        myMapOfHiddenElements.insert(std::make_pair(aTag,HHiddenElementsIDsForTag));
    }
}
