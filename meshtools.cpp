//! ----------------
//! custom includes
//! ----------------
#include "meshtools.h"
#include "mydefines.h"
#include "stlapiwriter.h"
#include <meshelementbycoords.h>
#include <tetgen.h>
#include "ais_colorscaleextended.h"
#include "tools.h"
#include "exportingtools.h"
#include "extendedrwstl.h"
#include "ccout.h"
#include <ng_meshvs_datasource1d.h>
#include "qprogressindicator.h"
#include "qprogressevent.h"
#include <elementtypes.h>
#include <isostrip.h>
#include <isostripbuilder.h>
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource3d.h>

//! ---
//! Qt
//! ---
#include <QApplication>
#include <QHash>
#include <QList>
#include <QDir>
#include <QMessageBox>

//! ----
//! OCC
//! ----
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfMapOfInteger.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Array1OfTriangle.hxx>
#include <gp_Pnt.hxx>
#include <MeshVS_EntityType.hxx>
#include <BRepMesh_ShapeTool.hxx>
#include <BRepTools.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TColgp_SequenceOfXYZ.hxx>
#include <MeshVS_NodalColorPrsBuilder.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DisplayModeFlags.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include <RWStl.hxx>
#include <StlMesh_Mesh.hxx>
#include <Poly_Triangulation.hxx>
#include <BRep_Tool.hxx>
#include <TColgp_Array1OfDir.hxx>
#include <Poly_Connect.hxx>
#include <CSLib_DerivativeStatus.hxx>
#include <CSLib_NormalStatus.hxx>
#include <CSLib.hxx>
#include <StlAPI_ErrorStatus.hxx>
#include <Graphic3d_MaterialAspect.hxx>
#include <Quantity_TypeOfColor.hxx>
#include <Aspect_TypeOfColorScaleData.hxx>

#include <OSD_Path.hxx>
#include <Message_ProgressIndicator.hxx>
#include <MeshVS_DataSource.hxx>
#include <Standard_ErrorHandler.hxx>

#include <GeomLib_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>

#include <MeshVS_ElementalColorPrsBuilder.hxx>


//! ------
//! EMESH
//! ------
#include <OMFTools.hxx>
#include <OMFDS_Mesh.hxx>
#include <OMFVS_DataSource.hxx>

//! -------
//! global
//! -------
#include "global.h"

//! -------
//! helper
//! -------
int hueFromValue(int theValue,int theMin,int theMax)
{
    int aMinLimit (0), aMaxLimit (230);
    int aHue = aMaxLimit;
    if (theMin!=theMax) aHue = (int)(aMaxLimit -(aMaxLimit-aMinLimit)*(theValue-theMin)/(theMax-theMin));
    aHue = std::min (std::max (aMinLimit, aHue), aMaxLimit);
    return aHue;
}

#define DELTA_TOL 1e-10

//! --------------------------------------------------------------------
//! function: toPolyTriangulation
//! details:  translate the MeshVS_DataSource into a Poly_Triangulation
//!           in case of second order elements the midside nodes are
//!           discarded
//! --------------------------------------------------------------------
bool MeshTools::toPolyTriangulation(const occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS,
                                    occHandle(Poly_Triangulation) &aTriangulation)
{
    cout<<"MeshTools::toPolyTriangulation->____function called____"<<endl;

    //! -------------
    //! sanity check
    //! -------------
    if(faceMeshDS.IsNull()) return false;
    if(faceMeshDS->GetAllElements().Extent()<1) return false;
    if(faceMeshDS->GetAllNodes().Extent()<3) return false;

    int NN = faceMeshDS->GetAllNodes().Extent();
    int NE = faceMeshDS->GetAllElements().Extent();

    //! ----------------------------
    //! array of triagles and nodes
    //! ----------------------------
    Poly_Array1OfTriangle trigs(1,NE);
    TColgp_Array1OfPnt nodes(1,NN);
    TColStd_IndexedMapOfInteger map;

    //! --------------------
    //! scan the mesh nodes
    //! --------------------
    int i=1;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(faceMeshDS->GetAllNodes());it.More();it.Next())
    {
        int NbNodes;
        MeshVS_EntityType eType;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);

        int globalNodeID = it.Key();
        faceMeshDS->GetGeom(globalNodeID,Standard_False,coords,NbNodes,eType);
        gp_Pnt V(coords.Value(1),coords.Value(2),coords.Value(3));
        nodes.SetValue(i,V);
        i++;
        map.Add(globalNodeID);
    }

    //! --------------------------
    //! scan the surface elements
    //! --------------------------
    i=1;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(faceMeshDS->GetAllElements());it.More();it.Next())
    {
        int globalElementID = it.Key();
        int NbNodes, nbuf[8];
        TColStd_Array1OfInteger nodeIDs(*nbuf,1,8);

        faceMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        switch(NbNodes)
        {
        case 3:
        {
            int V1 = nodeIDs.Value(1);
            int V2 = nodeIDs.Value(2);
            int V3 = nodeIDs.Value(3);

            //! --------------------
            //! create one triangle
            //! --------------------
            Poly_Triangle trig;
            trig.Set(map.FindIndex(V1),map.FindIndex(V2),map.FindIndex(V3));
            trigs.SetValue(i,trig);
            i++;
        }
            break;

        case 4:
        {
            int V1 = nodeIDs.Value(1);
            int V2 = nodeIDs.Value(2);
            int V3 = nodeIDs.Value(3);
            int V4 = nodeIDs.Value(4);

            //! ---------------------
            //! create two triangles
            //! ---------------------
            Poly_Triangle trig, trig1;
            trig.Set(map.FindIndex(V1),map.FindIndex(V2),map.FindIndex(V3));
            trig1.Set(map.FindIndex(V1),map.FindIndex(V3),map.FindIndex(V4));
            trigs.SetValue(i,trig);
            i++;
            trigs.SetValue(i,trig1);
            i++;
        }
            break;
        }
    }

    //! -------------------------
    //! create the triangulation
    //! -------------------------
    aTriangulation = new Poly_Triangulation(nodes,trigs);
    return true;
}

//! ---------------------------
//! function: storeMeshToShape
//! details:
//! ---------------------------
bool MeshTools::storeMeshToShape(const TopoDS_Shape &shape, occHandle(Poly_Triangulation) &theTriangulation)
{
    if(shape.IsNull()) return false;
    for(TopExp_Explorer anExp(shape,TopAbs_FACE);anExp.More();anExp.Next())
    {
        TopoDS_Face curFace = TopoDS::Face(anExp.Current());
        BRepMesh_ShapeTool::AddInFace(curFace, theTriangulation);
    }
    return true;
}

//! -------------------------------------------------------------------
//! function: OCCDSToNegtenSurfaceMesh
//! details:  convert the EMESH surface mesh in the Netgen mesh format
//! -------------------------------------------------------------------
Ng_Mesh* MeshTools::OCCDSToNegtenSurfaceMesh(const TopoDS_Shape &theShape,
                                             const NCollection_Handle<NCollection_List<Standard_Integer>> &badFacesIndices)
{
    cout<<"MeshTools::OCCDSToNegtenSurfaceMesh->____function called____"<<endl;

    //! -------------------
    //! init a Netgen mesh
    //! -------------------
    Ng_Mesh *NgMesh = Ng_NewMesh();

    //! -------------------------------------------------
    //! retrieve the overall surface mesh from the shape
    //! -------------------------------------------------
    Standard_Boolean isComputeNormals = Standard_True;
    Standard_Boolean isGroup = Standard_False;
    occHandle(OMFDS_Mesh) theMesh;
    bool isDone = OMFTools::TriangulationToMesh(theShape,theMesh,isComputeNormals,isGroup);
    if(!isDone)
    {
        Ng_DeleteMesh(NgMesh);
        return NULL;
    }

    Handle (OMFVS_DataSource) theMeshVS_DataSource = new OMFVS_DataSource(theMesh, Standard_True);

    //! --------------------------------
    //! the same: change option isGroup
    //! --------------------------------
    occHandle(OMFDS_Mesh) meshWithSubMeshes;
    isGroup = Standard_True;
    isDone = OMFTools::TriangulationToMesh(theShape,meshWithSubMeshes,isComputeNormals,isGroup);
    if(!isDone)
    {
        Ng_DeleteMesh(NgMesh);
        return NULL;
    }
    occHandle(OMFVS_DataSource) OMFVS_FaceDS = new OMFVS_DataSource(meshWithSubMeshes,Standard_False);

    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(theShape,TopAbs_FACE,faceMap);
    int NbTopologyFaces = faceMap.Extent();

    //! --------------
    //! add the nodes
    //! --------------
    TColStd_PackedMapOfInteger nodeMap = theMeshVS_DataSource->GetAllNodes();
    TColStd_MapIteratorOfPackedMapOfInteger anIter;
    TColStd_IndexedMapOfInteger indexedNodeMap;
    double P[3];
    for(anIter.Initialize(nodeMap);anIter.More();anIter.Next())
    {
        int nodeID = anIter.Key();
        indexedNodeMap.Add(nodeID);

        int NbNodes;
        double aCoordsBuf[3];
        TColStd_Array1OfReal Coords(*aCoordsBuf,1,3);
        MeshVS_EntityType theType;
        OMFVS_FaceDS->GetGeom(nodeID,Standard_False,Coords,NbNodes,theType);
        P[0] = Coords.Value(1);
        P[1] = Coords.Value(2);
        P[2] = Coords.Value(3);
        Ng_AddPoint(NgMesh,P);
    }

    //! ----------------------------
    //! in this position this works
    //! ----------------------------
    Ng_ClearFaceDescriptors(NgMesh);
    int offset = 0;
    for(TopExp_Explorer anExp(theShape,TopAbs_FACE);anExp.More();anExp.Next())
    {
        int faceNr = faceMap.FindIndex(TopoDS::Face(anExp.Current()));

        //! ---------------------
        //! the mesh of the face
        //! ---------------------
        if(badFacesIndices->Contains(faceNr))
        {
            cout<<"MeshTools::OCCDSToNegtenSurfaceMesh()->____jumping over face nr. "<<faceNr<<" which cannot be meshed____"<<endl;
            ccout(QString("MeshTools::OCCDSToNegtenSurfaceMesh()->____jumping over face nr. %1 which cannot be meshed____").arg(faceNr));

            Ng_SetupFaceDescriptor(NgMesh,faceNr,1,0,NbTopologyFaces);
            offset++;
            //!cout<<"offset: "<<offset<<endl;
            continue;
        }
        OMFVS_FaceDS->SetCurrentMesh(faceNr-offset);

        //! ---------------------------------------------------
        //! Ng_SetupFaceDescriptos in this position is correct
        //! ---------------------------------------------------
        Ng_SetupFaceDescriptor(NgMesh,faceNr,1,0,NbTopologyFaces);

        //! --------------------------
        //! scan the surface elements
        //! --------------------------
        int N[3];
        TColStd_PackedMapOfInteger eMap = OMFVS_FaceDS->GetAllElements();
        for(anIter.Initialize(eMap);anIter.More();anIter.Next())
        {
            int eID = anIter.Key();
            Standard_Integer aNodeIDBuf[3];
            TColStd_Array1OfInteger NodeIDs(*aNodeIDBuf,1,3);
            int NbNodes;
            OMFVS_FaceDS->GetNodesByElement(eID,NodeIDs,NbNodes);
            N[0]=indexedNodeMap.FindIndex(NodeIDs.Value(1));
            N[1]=indexedNodeMap.FindIndex(NodeIDs.Value(2));
            N[2]=indexedNodeMap.FindIndex(NodeIDs.Value(3));
            Ng_AddSurfaceElementWithIndex(NgMesh,NG_TRIG,N,faceNr);
        }
    }

    cout<<"MeshTools::OCCDSToNegtenSurfaceMesh->____Express mesh to Netgen mesh conversion done____"<<endl;

    //cout<<"MestTools::OCCDSToNegtenSurfaceMesh->____saving Meshreconstruction_fromOCC.vol____"<<endl;
    //Ng_SaveMesh(NgMesh,"D:\\MeshReconstruction.vol");

    return NgMesh;
}

//! -----------------------------------------------------------------------
//! function: buildColoredMesh
//! details:  create a MeshVS_Mesh object with colored shaded presentation
//! -----------------------------------------------------------------------
bool MeshTools::buildColoredMesh(const occHandle(MeshVS_DataSource) &theMeshVS_DataSource,
                                 const std::map<int,double> &res,
                                 occHandle(MeshVS_Mesh) &aColoredMesh,
                                 double min,
                                 double max,
                                 int numberOfLevels,
                                 bool showEdges,
                                 bool autoscale)
{
    if(theMeshVS_DataSource.IsNull()) return false;

    cout<<"MeshTools::buildColoredMesh->____function called____"<<endl;

    //! -----------------------------------------
    //! build the MeshVS_Mesh interactive object
    //! -----------------------------------------
    aColoredMesh = new MeshVS_Mesh();
    aColoredMesh->SetDataSource(theMeshVS_DataSource);

    //! ----------------
    //! colored builder
    //! ----------------
    occHandle(MeshVS_NodalColorPrsBuilder) nodalColorBuilder =
            new MeshVS_NodalColorPrsBuilder(aColoredMesh, MeshVS_DMF_NodalColorDataPrs | MeshVS_DMF_OCCMask);
    nodalColorBuilder->UseTexture(true);

    //! ----------------------
    //! prepare the color map
    //! ----------------------
    Aspect_SequenceOfColor aColorMap;

    //! -------------------------
    //! this activates autoscale
    //! -------------------------
    int indexOfMin, indexOfMax;
    if(autoscale)
    {
        max = -1e10;
        min = -max;

        for(std::map<int,double>::const_iterator it = res.cbegin(); it!=res.cend(); it++)
        {
            double curVal = it->second;
            int index = it->first;
            if(curVal>=max) { max = curVal; indexOfMax = index; }
            if(curVal<=min) { min = curVal; indexOfMin = index; }
        }
    }

    double range = max-min;
    if(range <=1e-6)
    {
        cout<<"--------------------------------------------------------------"<<endl;
        cout<<" almost constant value found when generating the colored mesh "<<endl;
        cout<<" expanding the scale range to [-10.0, 10.0]                   "<<endl;
        cout<<"--------------------------------------------------------------"<<endl;
        max = max+10.0;
        min = min-10.0;
        range = max - min;
    }
    double delta = range/numberOfLevels;

    //! ---------------------------------------------------------------
    //! Ex: numberOfLevels = 3
    //! => delta = Delta/3
    //! => i = {0, 1, 2, 3}
    //! => val = int(1*(delta/Delta)*360)
    //! => val = int(2*(delta/Delta)*360)
    //! => val = int(3*(delta/Delta)*360) = int(3*(Delta/3/Delta)*360)
    //! ---------------------------------------------------------------
    for(int i=0; i<=numberOfLevels; i++)
    {
        int val = int((i*delta)*(360.0/range));
        int hue = hueFromValue(val,0,360);
        Quantity_Color aColor(hue,1.0,1.0,Quantity_TOC_HLS);
        aColorMap.Append(aColor);
    }

    //! ------------------------------------------------------------------
    //! assign color scale map  values (0.0 <-> 1.0) to nodes
    //! iterate through the nodes and add an appropriate value to the map
    //! scan the result of type (nodeID, scalarValue)
    //! ------------------------------------------------------------------
    std::map<int, double>::const_iterator itNodes;

    TColStd_DataMapOfIntegerReal aScaleMap;

    for(itNodes = res.cbegin(); itNodes!= res.cend(); ++itNodes)
    {
        int nodeID = itNodes->first;
        double aValue;
        if(range != 0.0) aValue = (itNodes->second-min)/range;
        else aValue = 0.0;
        aScaleMap.Bind(nodeID, aValue);
    }

    //! -----------------------------------------------------
    //! pass color map and color scale values to the builder
    //! -----------------------------------------------------
    nodalColorBuilder->SetColorMap(aColorMap);
    nodalColorBuilder->SetInvalidColor (Quantity_NOC_VIOLET);
    nodalColorBuilder->SetTextureCoords(aScaleMap);

    aColoredMesh->AddBuilder(nodalColorBuilder,Standard_False);

    //! ---------------------------------------
    //! configure the drawer: other properties
    //! ---------------------------------------
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DMF_Shading, true);
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, false);
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, showEdges);
    aColoredMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);

    cout<<"MeshTools::buildColoredMesh->____exiting function____"<<endl;

    return true;
}

//! -----------------------------
//! rimuovere il prima possibile
//! -----------------------------
extern int hueFromValue(int,int,int);

//! ------------------------
//! function: buildIsoStrip
//! details:
//! ------------------------
bool MeshTools::buildIsoStrip(const occHandle(MeshVS_DataSource) &theMeshDS,        //! input mesh data source
                              const std::map<int,double> &res,                      //! nodal results
                              double min,                                           //! used for isostrips generation
                              double max,                                           //! used for isostrips generation
                              int NbLevels,                                         //! user for isostrip generation
                              occHandle(MeshVS_Mesh) &aColoredMesh,                 //! result
                              bool showEdges)
{
    Q_UNUSED(showEdges)
    if(theMeshDS.IsNull()) return false;

    //! ---------------------
    //! prepare the isostrip
    //! ---------------------
    cout<<"@ --------------------------"<<endl;
    cout<<"@ - preparing isostrips "<<endl;

    std::vector<isoStrip> vecIsoStrip;
    vecIsoStrip.push_back(isoStrip(-1e10,min));
    double delta = (max-min)/NbLevels;
    for(int n = 0; n<NbLevels; n++)
    {
        double ys = min + n*delta;
        double ye = ys + delta;
        isoStrip anIsoStrip(ys,ye);
        vecIsoStrip.push_back(anIsoStrip);
        cout<<"@ - "<<ys<<"\t"<<ye<<endl;
    }
    vecIsoStrip.push_back(isoStrip(max,1e10));
    cout<<"@ --------------------------"<<endl;

    //! -------------------------------------------------
    //! run the isostrip builder on the mesh data source
    //! -------------------------------------------------
    isoStripBuilder anIsoStripBuilder;
    anIsoStripBuilder.setMeshDataSource(theMeshDS);
    anIsoStripBuilder.setValues(res);
    anIsoStripBuilder.setIsoStrips(vecIsoStrip);

    std::vector<meshElementByCoords> allElements;
    //bool isDone = anIsoStripBuilder.perform(allElements);
    bool isDone = anIsoStripBuilder.perform1(allElements);
    if(isDone == false)
    {
        cout<<"MeshTools::buildIsoStrip()->____cannot build isostrips____"<<endl;
        return false;
    }
    occHandle(Ng_MeshVS_DataSourceFace) finalMesh = new Ng_MeshVS_DataSourceFace(allElements,true,true);
    if(finalMesh.IsNull())
    {
        cout<<"MeshTools::buildIsoStrip()->____cannot build isostrips____"<<endl;
        return false;
    }
    cout<<"@ --------------------------"<<endl;
    cout<<"@ - overall strip mesh"<<endl;
    cout<<"@ - elements: "<<finalMesh->GetAllElements().Extent()<<endl;
    cout<<"@ - nodes: "<<finalMesh->GetAllNodes().Extent()<<endl;
    cout<<"@ --------------------------"<<endl;

    //! -----------------------------------------
    //! build the MeshVS_Mesh interactive object
    //! -----------------------------------------
    aColoredMesh = new MeshVS_Mesh();
    aColoredMesh->SetDataSource(finalMesh);
    occHandle(MeshVS_ElementalColorPrsBuilder) aPrsBuilder =
            new MeshVS_ElementalColorPrsBuilder(aColoredMesh, MeshVS_DMF_ElementalColorDataPrs | MeshVS_DMF_OCCMask);

    int n = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(finalMesh->GetAllElements()); it.More(); it.Next(), n++)
    {
        int isoStripNb = allElements[n].ID;
        int hue = hueFromValue(isoStripNb,1,(int)vecIsoStrip.size()-2);
        if(isoStripNb<1 || isoStripNb > vecIsoStrip.size()-2) hue = 310;
        Quantity_Color aColor(hue,1.0,1.0,Quantity_TOC_HLS);
        aPrsBuilder->SetColor1(it.Key(),aColor);
    }
    aColoredMesh->AddBuilder(aPrsBuilder);
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, false);
    aColoredMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, false);

    cout<<"MeshTools::buildIsoStrip()->____exiting function____"<<endl;
    return true;
}

//! ---------------------------
//! function: buildIsoSurfaces
//! details:
//! ---------------------------
bool MeshTools::buildIsoSurfaces(const occHandle(MeshVS_DataSource) &theMeshDS,        //! input mesh data source
                                 const std::map<int,double> &res,                          //! nodal results
                                 double min,                                           //! used for isostrips generation
                                 double max,                                           //! used for isostrips generation
                                 int NbLevels,                                         //! user for isostrip generation
                                 occHandle(MeshVS_Mesh) &aColoredMesh,                 //! result
                                 bool showEdges)
{
    if(theMeshDS.IsNull()) return false;

    //! ---------------------
    //! prepare the isostrip
    //! ---------------------
    cout<<"@ --------------------------"<<endl;
    cout<<"@ - preparing isostrips "<<endl;

    std::vector<isoStrip> vecIsoStrip;
    vecIsoStrip.push_back(isoStrip(-1e10,min));
    double delta = (max-min)/NbLevels;
    for(int n = 0; n<NbLevels; n++)
    {
        double ys = min + n*delta;
        double ye = ys + delta;
        isoStrip anIsoStrip(ys,ye);
        vecIsoStrip.push_back(anIsoStrip);
        cout<<"@ - "<<ys<<"\t"<<ye<<endl;
    }
    vecIsoStrip.push_back(isoStrip(max,1e10));
    cout<<"@ --------------------------"<<endl;

    //! -------------------------------------------------
    //! run the isostrip builder on the mesh data source
    //! -------------------------------------------------
    isoStripBuilder anIsoStripBuilder;
    anIsoStripBuilder.setMeshDataSource(theMeshDS);
    anIsoStripBuilder.setValues(res);
    anIsoStripBuilder.setIsoStrips(vecIsoStrip);

    int NbIsoSurfaces = NbLevels;
    std::map<int,int> mapElementLevel;
    std::vector<meshElementByCoords> allElements;
    bool isDone = anIsoStripBuilder.performIsoSurface(NbIsoSurfaces,allElements,mapElementLevel);
    if(isDone == false) return false;

    occHandle(Ng_MeshVS_DataSourceFace) finalMesh = new Ng_MeshVS_DataSourceFace(allElements,true,true);

    cout<<"@ --------------------------"<<endl;
    cout<<"@ - overall strip mesh"<<endl;
    cout<<"@ - elements: "<<finalMesh->GetAllElements().Extent()<<endl;
    cout<<"@ - nodes: "<<finalMesh->GetAllNodes().Extent()<<endl;
    cout<<"@ --------------------------"<<endl;

    //! -----------------------------------------
    //! build the MeshVS_Mesh interactive object
    //! -----------------------------------------
    aColoredMesh = new MeshVS_Mesh();
    aColoredMesh->SetDataSource(finalMesh);
    occHandle(MeshVS_ElementalColorPrsBuilder) aPrsBuilder =
            new MeshVS_ElementalColorPrsBuilder(aColoredMesh, MeshVS_DMF_ElementalColorDataPrs | MeshVS_DMF_OCCMask);

    int n = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(finalMesh->GetAllElements()); it.More(); it.Next(), n++)
    {
        std::map<int,int>::iterator itmap = mapElementLevel.find(it.Key());
        if(itmap==mapElementLevel.end()) continue;
        int isoStripNb = itmap->second;
        int hue = hueFromValue(isoStripNb,1,(int)vecIsoStrip.size()-2);
        if(isoStripNb<1 || isoStripNb > vecIsoStrip.size()-2) hue = 310;
        Quantity_Color aColor(hue,1.0,1.0,Quantity_TOC_HLS);
        aPrsBuilder->SetColor1(it.Key(),aColor);
    }
    aColoredMesh->AddBuilder(aPrsBuilder);
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, false);
    aColoredMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, showEdges);
    return true;
}

//! ------------------------
//! function: buildIsoStrip
//! details:
//! ------------------------
bool MeshTools::buildIsoStrip(const occHandle(MeshVS_DataSource) &theMeshDS,        //! input mesh data source
                              const std::map<int,double> &res,                          //! nodal results
                              const std::map<int,gp_Vec> &displacementMap,              //! displacement field
                              double scale,
                              double min,                                           //! used for isostrips generation
                              double max,                                           //! used for isostrips generation
                              int NbLevels,                                         //! user for isostrip generation
                              occHandle(MeshVS_Mesh) &aColoredMesh,                 //! result
                              bool showEdges)                                       //! option
{
    if(theMeshDS.IsNull()) return false;

    //! --------------------------------------------
    //! adjust the scale in case of negative values
    //! --------------------------------------------
    if(scale<0.0) scale = 1.0;

    //! ---------------------
    //! prepare the isostrip
    //! ---------------------
    cout<<"@ --------------------------"<<endl;
    cout<<"@ - preparing isostrips "<<endl;

    std::vector<isoStrip> vecIsoStrip;
    vecIsoStrip.push_back(isoStrip(-1e10,min));
    double delta = (max-min)/NbLevels;
    for(int n = 0; n<NbLevels; n++)
    {
        double ys = min + n*delta;
        double ye = ys + delta;
        isoStrip anIsoStrip(ys,ye);
        vecIsoStrip.push_back(anIsoStrip);
        cout<<"@ - "<<ys<<"\t"<<ye<<endl;
    }
    vecIsoStrip.push_back(isoStrip(max,1e10));
    cout<<"@ --------------------------"<<endl;

    //! --------------------------------------------------------------
    //! the deformed mesh data source
    //! construct with an empty mesh data base of a specialized class
    //! --------------------------------------------------------------
    occHandle(MeshVS_DeformedDataSource) theDeformedDS = new MeshVS_DeformedDataSource(occHandle(Ng_MeshVS_DataSourceFace)(),1.0);

    cout<<"MeshTools::buildIsoStrip()->____input mesh. Elements: "<<theMeshDS->GetAllElements().Extent()<<" Nodes: "<<theMeshDS->GetAllNodes().Extent()<<"____"<<endl;

    theDeformedDS->SetNonDeformedDataSource(theMeshDS);
    for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshDS->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        const gp_Vec &d = displacementMap.at(globalNodeID);
        theDeformedDS->SetVector(globalNodeID,d);
    }

    //! ----------------
    //! apply the scale
    //! ----------------
    theDeformedDS->SetMagnify(scale);

    //! ----------------------------------------------------------
    //! run the isostrip builder on the deformed mesh data source
    //! ----------------------------------------------------------
    isoStripBuilder anIsoStripBuilder;
    anIsoStripBuilder.setMeshDataSource(theDeformedDS);
    anIsoStripBuilder.setValues(res);
    anIsoStripBuilder.setIsoStrips(vecIsoStrip);

    std::vector<meshElementByCoords> allElements;
    //bool isDone = anIsoStripBuilder.perform(allElements);
    bool isDone = anIsoStripBuilder.perform1(allElements);

    Q_UNUSED (isDone)

    occHandle(Ng_MeshVS_DataSourceFace) finalMesh = new Ng_MeshVS_DataSourceFace(allElements,true,true);

    cout<<"@ --------------------------"<<endl;
    cout<<"@ - overall strip mesh"<<endl;
    cout<<"@ - elements: "<<finalMesh->GetAllElements().Extent()<<endl;
    cout<<"@ - nodes: "<<finalMesh->GetAllNodes().Extent()<<endl;
    cout<<"@ --------------------------"<<endl;

    //! -----------------------------------------
    //! build the MeshVS_Mesh interactive object
    //! -----------------------------------------
    aColoredMesh = new MeshVS_Mesh();
    aColoredMesh->SetDataSource(finalMesh);

    occHandle(MeshVS_ElementalColorPrsBuilder) aPrsBuilder = new MeshVS_ElementalColorPrsBuilder(aColoredMesh, MeshVS_DMF_ElementalColorDataPrs | MeshVS_DMF_OCCMask);

    int n = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(finalMesh->GetAllElements()); it.More(); it.Next(), n++)
    {
        int isoStripNb = allElements[n].ID;
        int hue = hueFromValue(isoStripNb,1,(int)vecIsoStrip.size()-2);
        if(isoStripNb<1 || isoStripNb > vecIsoStrip.size()-2) hue = 310;
        Quantity_Color aColor(hue,1.0,1.0,Quantity_TOC_HLS);
        aPrsBuilder->SetColor1(it.Key(),aColor);
    }
    aColoredMesh->AddBuilder(aPrsBuilder);
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, false);
    aColoredMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
    aColoredMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, false);

    return true;
}

//! -----------------
//! function: Normal
//! details:
//! -----------------
static void Normal(const TopoDS_Face& aFace, Poly_Connect& pc, TColgp_Array1OfDir& Nor)
{
    const occHandle(Poly_Triangulation)& T = pc.Triangulation();
    BRepAdaptor_Surface S;
    Standard_Boolean hasUV = T->HasUVNodes();
    Standard_Integer i;
    TopLoc_Location l;
    occHandle(Geom_Surface) GS = BRep_Tool::Surface(aFace, l);

    if (hasUV && !GS.IsNull())
    {
        Standard_Boolean OK = Standard_True;
        gp_Vec D1U,D1V;
        gp_Vec D2U,D2V,D2UV;
        gp_Pnt P;
        Standard_Real U, V;
        CSLib_DerivativeStatus Status;
        CSLib_NormalStatus NStat;
        S.Initialize(aFace, Standard_False);
        const TColgp_Array1OfPnt2d& UVNodes = T->UVNodes();
        if (S.GetType() != GeomAbs_Plane)
        {
            for (i = UVNodes.Lower(); i <= UVNodes.Upper(); i++)
            {
                U = UVNodes(i).X();
                V = UVNodes(i).Y();
                S.D1(U,V,P,D1U,D1V);
                CSLib::Normal(D1U,D1V,Precision::Angular(),Status,Nor(i));
                if (Status != CSLib_Done)
                {
                    S.D2(U,V,P,D1U,D1V,D2U,D2V,D2UV);
                    CSLib::Normal(D1U,D1V,D2U,D2V,D2UV,Precision::Angular(),OK,NStat,Nor(i));
                }
                if (aFace.Orientation() == TopAbs_REVERSED) (Nor(i)).Reverse();
            }
        }
        else
        {
            gp_Dir NPlane;
            U = UVNodes(UVNodes.Lower()).X();
            V = UVNodes(UVNodes.Lower()).Y();
            S.D1(U,V,P,D1U,D1V);
            CSLib::Normal(D1U,D1V,Precision::Angular(),Status,NPlane);
            if (Status != CSLib_Done)
            {
                S.D2(U,V,P,D1U,D1V,D2U,D2V,D2UV);
                CSLib::Normal(D1U,D1V,D2U,D2V,D2UV,Precision::Angular(),OK,NStat,NPlane);
            }
            if (aFace.Orientation() == TopAbs_REVERSED) NPlane.Reverse();
            Nor.Init(NPlane);
        }
    }
    else
    {
        const TColgp_Array1OfPnt& Nodes = T->Nodes();
        Standard_Integer n[3];
        const Poly_Array1OfTriangle& triangles = T->Triangles();

        for (i = Nodes.Lower(); i <= Nodes.Upper(); i++)
        {
            gp_XYZ eqPlan(0, 0, 0);
            for (pc.Initialize(i);  pc.More(); pc.Next())
            {
                triangles(pc.Value()).Get(n[0], n[1], n[2]);
                gp_XYZ v1(Nodes(n[1]).Coord()-Nodes(n[0]).Coord());
                gp_XYZ v2(Nodes(n[2]).Coord()-Nodes(n[1]).Coord());
                eqPlan += (v1^v2).Normalized();
            }
            Nor(i) = gp_Dir(eqPlan);
            if (aFace.Orientation() == TopAbs_REVERSED) (Nor(i)).Reverse();
        }
    }
}

//! ---------------------------------------------------------------------------
//! function: toSTLMesh
//! details:  take a shape and extract the STL mesh of the shaded presentation
//!           it uses disk => change ... to do ...
//! ---------------------------------------------------------------------------
bool MeshTools::toSTLMesh1(const TopoDS_Shape &shape,
                           const QString &surfaceMeshFilePath,
                           occHandle(Ng_MeshVS_DataSource2D) &surfaceMeshDS,
                           std::vector<occHandle(Ng_MeshVS_DataSourceFace)> &vecFaceMeshDS,
                           QProgressIndicator *aProgressIndicator,
                           int done)
{
    //! ----------------------------------------------
    //! This generates the overall 2D mesh datasource
    //! An .stl (extended) format is used, with tags
    //! ----------------------------------------------
    TopTools_IndexedMapOfShape M;
    TopExp::MapShapes(shape,TopAbs_FACE,M);
    int NbFaces = M.Extent();

    //! ----------------------
    //! init the progress bar
    //! ----------------------
    if(aProgressIndicator!=Q_NULLPTR)
    {
        //! --------------------------------
        //! init the secondary progress bar
        //! --------------------------------
        QProgressEvent *pe = new QProgressEvent();
        pe->setVal(done);
        pe->setMessage("Generating the face mesh data sources");
        pe->setAction1(QProgressEventAction::QProgressEvent_Init);
        pe->setRange1(0, NbFaces);
        QApplication::postEvent(aProgressIndicator,pe);
    }

    //! ------------------------------------------------
    //! init a vector of (NbFaces+1) empty StlMesh_Mesh
    //! ------------------------------------------------
    std::vector<occHandle(StlMesh_Mesh)> vecFaceStlMesh;
    for(int faceNr=0; faceNr<=NbFaces; faceNr++)
    {
        occHandle(StlMesh_Mesh) aStlMesh = new StlMesh_Mesh();
        vecFaceStlMesh.push_back(aStlMesh);
    }

    //! ---------------------------------------------------------------------------------------
    //! return:
    //! - the overall .stl mesh
    //! - the vector of the .stl meshes of the faces
    //! - for each face the local to global numbering maps, for the nodes and for the elements
    //! ---------------------------------------------------------------------------------------
    QList<QMap<int,int>> maps_localToGlobal_nodeIDs;
    QList<QMap<int,int>> maps_localToGlobal_elementIDs;
    occHandle(StlMesh_Mesh) overallSTL = ExtendedRWStl::ReadExtendedSTLAscii(surfaceMeshFilePath,
                                                                             vecFaceStlMesh,
                                                                             maps_localToGlobal_nodeIDs,
                                                                             maps_localToGlobal_elementIDs);
    //! -------------
    //! sanity check
    //! -------------
    if(overallSTL.IsNull()) return false;
    if(overallSTL->NbTriangles()<1 || overallSTL->NbVertices()<3) return false;

    //! -----------------------------------
    //! build the surface mesh data source
    //! -----------------------------------
    surfaceMeshDS = new Ng_MeshVS_DataSource2D(overallSTL);

    //! ---------------------------------------------
    //! handle also an error in the mesh constructor
    //! ---------------------------------------------
    if(surfaceMeshDS.IsNull()) return false;

    //! -----------------------------------------------------
    //! generate face mesh data sources using the stl meshes
    //! -----------------------------------------------------
    for(int faceNr = 0; faceNr<=NbFaces; faceNr++)
    {
        occHandle(Ng_MeshVS_DataSourceFace) faceMeshDS;
        occHandle(StlMesh_Mesh) curStlMesh = vecFaceStlMesh.at(faceNr);
        if(curStlMesh.IsNull())
        {
            vecFaceMeshDS[faceNr]=faceMeshDS;
            continue;
        }
        if(curStlMesh->NbTriangles()<=0)
        {
            vecFaceMeshDS[faceNr]=faceMeshDS;
            continue;
        }
        faceMeshDS = new Ng_MeshVS_DataSourceFace(curStlMesh,maps_localToGlobal_nodeIDs,maps_localToGlobal_elementIDs,faceNr);
        //cout<<"MeshTools::toSTLMesh1()->____face nr: "<<faceNr<<" nodes: "<<faceMeshDS->GetAllNodes().Extent()<<"____"<<endl;
        vecFaceMeshDS[faceNr]=faceMeshDS;
    }
    return true;
}

//! -------------------------------------------------------------------
//! function: arrayOfFaceDataSourcesToExtendedStlFile
//! details:  convert an array of face datasources to an extended .stl
//!           file and write it on disk
//! -------------------------------------------------------------------
bool MeshTools::faceDataSourcesToExtendedStlFile(const std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &mapOfFaceMeshDS,
                                                        const std::string &extendedStlFileName)
{
    cout<<"MeshTools::faceDataSourcesToExtendedStlFile()->____function called____"<<endl;

    ofstream s;
    s.open(extendedStlFileName);
    if(!s.is_open()) return false;

    s<<"solid"<<endl;

    for(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>::const_iterator it = mapOfFaceMeshDS.cbegin(); it!=mapOfFaceMeshDS.cend(); it++)
    {
        int faceNr = it->first;
        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceDS = it->second;

        cout<<"____arrayOfFaceDataSourcesToExtendedStlFile()->____working on faceNr: "<<faceNr<<"____"<<endl;

        if(curFaceDS.IsNull()) continue;

        char line[512];
        int NbNodes;
        int k=0;
        double bufn[9];
        double x1,x2,x3,y1,y2,y3,z1,z2,z3;

        TColStd_Array1OfReal coords(*bufn,1,9);
        MeshVS_EntityType type;

        TColStd_PackedMapOfInteger eMap = curFaceDS->GetAllElements();
        TColStd_MapIteratorOfPackedMapOfInteger eIter;

        for(eIter.Initialize(eMap);eIter.More();eIter.Next(),k++)
        {
            int elemID = eIter.Key();
            bool isOK = curFaceDS->GetGeom(elemID,true,coords,NbNodes,type);
            if(isOK)
            {
                x1 = coords(1); y1 = coords(2); z1 = coords(3);
                x2 = coords(4); y2 = coords(5); z2 = coords(6);
                x3 = coords(7); y3 = coords(8); z3 = coords(9);

                double s21x = x2-x1; double s21y = y2-y1; double s21z = z2-z1;
                double s31x = x3-x1; double s31y = y3-y1; double s31z = z3-z1;
                //! ---------------------
                //! i       j       k
                //! s21x    s21y    s21z
                //! s31x    s31y    s31z
                //! ---------------------
                double n1 = s21y*s31z-s21z*s31y;
                double n2 = s21z*s31x-s21x*s31z;
                double n3 = s21x*s31y-s21y*s31x;

                sprintf(line,"facet normal %.6e %.6e %.6e\n"
                             "outer loop\n"
                             "vertex %.6e %.6e %.6e\n"
                             "vertex %.6e %.6e %.6e\n"
                             "vertex %.6e %.6e %.6e\n"
                             "endloop\n"
                             "%d\n"
                             "endfacet\n",
                        n1,n2,n3,
                        x1, y1, z1,
                        x2, y2, z2,
                        x3, y3, z3,
                        faceNr);

                s<<line;
            }
            else continue;
        }
    }
    s<<"end solid"<<endl;
    s.close();
    return true;
}

//! ---------------------------
//! function: sortPointsOnEdge
//! details:
//! ---------------------------
void MeshTools::sortPointsOnEdge(const QList<mesh::meshPoint> &pointsOfTheEdge,
                                 const TopoDS_Edge &anEdge,
                                 QList<mesh::meshPoint> &sortedPointsOfTheEdge)
{
    QMap<double, mesh::meshPoint> pointsMap;

    double edgeTol = BRep_Tool::Tolerance(anEdge);
    double U_start, U_end, U;
    const occHandle(Geom_Curve) &aCurve = BRep_Tool::Curve(anEdge,U_start,U_end);
    //!cout<<"MeshTools::sortPointsOnEdge()->____U start: "<<U_start<<" U end: "<<U_end<<"____"<<endl;

    for(QList<mesh::meshPoint>::const_iterator it = pointsOfTheEdge.cbegin(); it!=pointsOfTheEdge.cend(); ++it)
    {
        mesh::meshPoint aPoint = *it;
        gp_Pnt curPoint(aPoint.x,aPoint.y,aPoint.z);

        //! -------------------------------------------------------------------------------------------------
        //! this algorithm seems to fails in some cases, maybe when the point is not exactly on the 3D curve
        //! -------------------------------------------------------------------------------------------------
        //if(GeomLib_Tool::Parameter(aCurve, curPoint, gp::Resolution(), U))
        if(GeomLib_Tool::Parameter(aCurve, curPoint, edgeTol, U))
        {
            pointsMap.insert(U,aPoint);
        }
        else
        {
            cout<<"____point "<<aPoint.ID<<"("<<aPoint.x<<", "<<aPoint.y<<", "<<aPoint.z<<") on edge not found. Trying projection: ";
            //GeomAPI_ProjectPointOnCurve projector(curPoint,aCurve,U_start,U_end);
            GeomAPI_ProjectPointOnCurve projector(curPoint,aCurve);
            if(projector.NbPoints()!=0)
            {
                U = projector.LowerDistanceParameter();
                pointsMap.insert(U,aPoint);
                cout<<" success____"<<endl;
            }
            else
            {
                cout<<" not found____"<<endl;
            }
        }
    }

    int n=0;
    for(QMap<double,mesh::meshPoint>::iterator it = pointsMap.begin(); it!=pointsMap.end(); ++it)
    {
        const mesh::meshPoint &aPoint = it.value();
        sortedPointsOfTheEdge<<aPoint;
        cout<<"_S__i: "<<++n<<" key: "<<aPoint.ID<<" x= "<<aPoint.x<<" y= "<<aPoint.y<<" z= "<<aPoint.z<<"____"<<endl;
    }
}

//! -----------------------------------------
//! function: buildPLC
//! details:  write a .node and a .poly file
//! -----------------------------------------
bool MeshTools::buildPLC(const NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &arrayOfFaceDS,
                         const QString &nodeFilePath,
                         const QString &polyFilePath,
                         QProgressIndicator *progressIndicator)
{
    cout<<"MeshTools::buildPLC->____function called____"<<endl;
    if(progressIndicator == Q_NULLPTR)
    {
        cout<<"MeshTools::buildPLC->____warning: the progress indicator is NULL. No progress will be shown____"<<endl;
    }

    QList<mesh::meshPoint> pointList_;
    QList<QVector<double>> pointList;
    int NbFacets = 0;
    int globalNodeNumber = 0;

    int startFaceNr = arrayOfFaceDS.Lower();
    int NbFaces = arrayOfFaceDS.Upper();

    //! -------------------
    //! send an init event
    //! -------------------
    if(progressIndicator!=Q_NULLPTR)
    {
        QProgressEvent *progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"Building the PLC",
                                                           QProgressEvent_Init,startFaceNr,NbFaces,1,"MeshTools building PLC on disk");
        QApplication::postEvent(progressIndicator,progressEvent);
        QApplication::processEvents();
    }

    //! ----------------------------------------
    //! iterate over the face mesh data sources
    //! ----------------------------------------
    for(int faceDSNr = startFaceNr; faceDSNr<=arrayOfFaceDS.Upper(); faceDSNr++)
    {
        //cout<<"MeshTools::buildPLC()->____node file: working on face data source nr: "<<faceDSNr<<"____"<<endl;

        //! -------------------
        //! check interruption
        //! -------------------
        if(Global::status().code==0)
        {
            cout<<"MeshTools::buildPLC()->____process interrupted by the user____"<<endl;
            return false;
        }

        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMeshDS = arrayOfFaceDS.Value(faceDSNr);
        if(curFaceMeshDS.IsNull()) continue;

        if(curFaceMeshDS->GetAllElements().Extent()>0)
        {
            int localElementID = 0;
            MeshVS_EntityType type;
            int NbNodes;
            double buf[24];
            TColStd_Array1OfReal coords(*buf,1,24);
            for(TColStd_MapIteratorOfPackedMapOfInteger eIt(curFaceMeshDS->GetAllElements()); eIt.More(); eIt.Next())
            {
                localElementID++;
                int globalElementID = eIt.Key();

                curFaceMeshDS->GetGeom(globalElementID,true,coords,NbNodes,type);
                for(int k=0; k<NbNodes; k++)
                {
                    int s = 3*k;
                    QVector<double> P {coords(s+1),coords(s+2),coords(s+3)};
                    //mesh::meshPoint aP {coords(s+1),coords(s+2),coords(s+3)}; new
                    //if(!pointList_.contains(aP))
                    //{
                    //    globalNodeNumber++;
                    //    pointList_<<aP;
                    //}
                    if(!pointList.contains(P))
                    {
                        globalNodeNumber++;
                        pointList<<P;
                    }
                }
            }
            NbFacets += curFaceMeshDS->GetAllElements().Extent();
        }

        //! ----------------------
        //! post a progress event
        //! ----------------------
        if(progressIndicator!=Q_NULLPTR)
        {
            QProgressEvent *progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"Building the PLC",
                                                               QProgressEvent_Update,0,9999,faceDSNr,"MeshTools building PLC on disk");
            QApplication::postEvent(progressIndicator,progressEvent);
            QApplication::processEvents();
        }
    }

    //! ------------------------------------------------
    //! reset the progress bar after generating the PLC
    //! ------------------------------------------------
    if(progressIndicator!=Q_NULLPTR)
    {
        QProgressEvent *progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"Building the PLC - done",
                                                           QProgressEvent_Update,0,9999,0,"MeshTools building PLC on disk");
        QApplication::postEvent(progressIndicator,progressEvent);
        QApplication::processEvents();
    }

    //! ---------------------
    //! write the .node file
    //! ---------------------
    FILE *nodeFile = fopen(nodeFilePath.toStdString().c_str(),"w");
    if(nodeFile==NULL)
    {
        cout<<"MeshTools::buildPLC()->____node section: error in generating the PLC____"<<endl;
        return false;
    }

    //! ----------------------------------------------------------------------------------------
    //! First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    //! ----------------------------------------------------------------------------------------
    fprintf(nodeFile,"%d\t3\t0\t1\n",pointList.length());

    for(int i=0; i<pointList.length(); i++)
    //for(int i=0; i<pointList_.length(); i++)  new
    {
        const QVector<double> &P = pointList.at(i);
        fprintf(nodeFile,"%d\t%.9e\t%.9e\t%.9e\t%d\n",i+1,P.at(0),P.at(1),P.at(2),1);

        //const mesh::meshPoint &aP_ = pointList_[i];
        //fprintf(nodeFile,"%d\t%.12e\t%.12e\t%.12e\t%d\n",i+1,aP_.x,aP_.y,aP_.z,1);
    }
    fclose(nodeFile);

    //! ---------------------
    //! write the .poly file
    //! ---------------------
    FILE *polyFile = fopen(polyFilePath.toStdString().c_str(),"w");
    if(polyFile==NULL)
    {
        cout<<"MeshTools::buildPLC()->____facet section: error in generating the PLC____"<<endl;
        return false;
    }

    //! ----------------------------------------------------------------------------------------
    //! First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    //! ----------------------------------------------------------------------------------------
    int NbPoints = 0;   //! this means that the nodes are listed in the separate file .node
    int dimension = 3;
    int NbAttributes = 0;
    int boundaryMarkerFlag = 1;

    fprintf(polyFile,"%d\t%d\t%d\t%d\n",NbPoints,dimension,NbAttributes,boundaryMarkerFlag);

    //! ----------------------------------------------------
    //! One line: <# of facets> <boundary markers (0 or 1)>
    //! ----------------------------------------------------
    fprintf(polyFile,"%d\t%d\n",NbFacets,1);
    for(int faceDSNr = arrayOfFaceDS.Lower(); faceDSNr<=arrayOfFaceDS.Upper(); faceDSNr++)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMeshDS = arrayOfFaceDS.Value(faceDSNr);
        if(curFaceMeshDS.IsNull())
        {
            cout<<"MeshTools::buildPLC()->____poly file: the face data source "<<faceDSNr<<" is null____"<<endl;
            continue;
        }

        int localElementID = 0;
        int NbNodes;
        MeshVS_EntityType type;
        double buf[24];
        TColStd_Array1OfReal coords(*buf,1,24);

        for(TColStd_MapIteratorOfPackedMapOfInteger eIt(curFaceMeshDS->GetAllElements());eIt.More();eIt.Next())
        {
            //! ---------------------------------------------------------
            //! One line: <# of polygons> [# of holes] [boundary marker]
            //! ---------------------------------------------------------
            fprintf(polyFile,"%d\t%d\t%d\n",1,0,faceDSNr);

            //! ----------------------------------------------------------------------------------------
            //! Following lines list # of polygons: <# of corners> <corner 1> <corner 2> ... <corner #>
            //! ----------------------------------------------------------------------------------------
            fprintf(polyFile,"%d\t",3);

            localElementID++;
            int globalElementID = eIt.Key();
            curFaceMeshDS->GetGeom(globalElementID,true,coords,NbNodes,type);

            for(int k=0; k<NbNodes-1; k++)
            {
                int s = 3*k;
                QVector<double> P { coords(s+1), coords(s+2), coords(s+3) };
                fprintf(polyFile,"%d\t",pointList.indexOf(P)+1);
            }

            int n = 3*(NbNodes-1);
            QVector<double> P { coords(n+1), coords(n+2), coords(n+3) };
            fprintf(polyFile,"%d\n",pointList.indexOf(P)+1);
        }
    }

    //! ---------
    //! no holes
    //! ---------
    fprintf(nodeFile,"#holes section\n");
    fprintf(polyFile,"%d\n",0);

    //! -----------
    //! no regions
    //! -----------
    fprintf(nodeFile,"#regions section\n");
    fprintf(polyFile,"%d\n",0);
    fclose(polyFile);
    return true;
}

//! -------------------------------------------------------------
//! function: surfaceMeshDSToNetgenMesh
//! details:  Netgen mesh pointer from a surface mesh datasource
//! -------------------------------------------------------------
Ng_Mesh* MeshTools::surfaceMeshDSToNetgenMesh(const occHandle(MeshVS_DataSource) &aSurfaceMeshDS)
{
    Ng_Init();
    Ng_Mesh *NgMesh = Ng_NewMesh();

    //! --------------
    //! add the nodes
    //! --------------
    TColStd_PackedMapOfInteger nodeMap = aSurfaceMeshDS->GetAllNodes();
    TColStd_MapIteratorOfPackedMapOfInteger anIter;
    TColStd_IndexedMapOfInteger indexedNodeMap;
    double P[3];
    for(anIter.Initialize(nodeMap);anIter.More();anIter.Next())
    {
        int nodeID = anIter.Key();
        indexedNodeMap.Add(nodeID);

        Standard_Integer NbNodes;
        Standard_Real aCoordsBuf[3];
        TColStd_Array1OfReal Coords(*aCoordsBuf,1,3);
        MeshVS_EntityType theType;
        aSurfaceMeshDS->GetGeom(nodeID,Standard_False,Coords,NbNodes,theType);
        P[0] = Coords.Value(1);
        P[1] = Coords.Value(2);
        P[2] = Coords.Value(3);
        Ng_AddPoint(NgMesh,P);
    }

    Ng_ClearFaceDescriptors(NgMesh);

    //! ---------------------------------------------------
    //! Ng_SetupFaceDescriptos in this position is correct
    //! ---------------------------------------------------
    Ng_SetupFaceDescriptor(NgMesh,1,0,0,0);

    //! -----------------
    //! add the elements
    //! -----------------
    int N[3];
    TColStd_PackedMapOfInteger eMap = aSurfaceMeshDS->GetAllElements();
    for(anIter.Initialize(eMap);anIter.More();anIter.Next())
    {
        int eID = anIter.Key();
        Standard_Integer aNodeIDBuf[3];
        TColStd_Array1OfInteger NodeIDs(*aNodeIDBuf,1,3);
        int NbNodes;
        aSurfaceMeshDS->GetNodesByElement(eID,NodeIDs,NbNodes);
        N[0]=indexedNodeMap.FindIndex(NodeIDs.Value(1));
        N[1]=indexedNodeMap.FindIndex(NodeIDs.Value(2));
        N[2]=indexedNodeMap.FindIndex(NodeIDs.Value(3));
        Ng_AddSurfaceElement(NgMesh,NG_TRIG,N);
    }

    //! -----------    
    //! diagnostic
    //! -----------
    Ng_Meshing_Parameters mp;
    Ng_GenerateVolumeMesh(NgMesh,&mp);
    Ng_SaveMesh(NgMesh,"D:/testSurfaceMesh.vol");
    Ng_DeleteMesh(NgMesh);
    //! ---------------
    //! end diagnostic
    //! ---------------

    return NgMesh;
}

//! -------------------------------------
//! function: filterVolumeElementsByType
//! details:
//! -------------------------------------
void MeshTools::filterVolumeElementsByType(const occHandle(MeshVS_DataSource) &inputMesh,
                                  ElemType theType,
                                  occHandle(MeshVS_DataSource) &outputMesh)
{
    if(inputMesh.IsNull()) return;
    if(inputMesh->GetAllElements().Extent()<1) return;
    if(inputMesh->GetAllNodes().Extent()<1) return;

    std::vector<meshElementByCoords> listOfElements;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(inputMesh->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        int buf[20], NbNodes;
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        inputMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        ElemType anElemType;
        switch(NbNodes)
        {
        case 4: anElemType = TET; break;
        case 5: anElemType = PYRAM; break;
        case 6: anElemType = PRISM; break;
        case 8: anElemType = HEXA; break;
        }
        if(anElemType!=theType) continue;

        meshElementByCoords ameshElementByCoords;
        ameshElementByCoords.type = anElemType;

        //QList<mesh::meshPoint> pointList;
        for(int n=1; n<=NbNodes; n++)
        {
            int globalNodeID = nodeIDs(n);
            //cout<<"____globalNodeID: "<<globalNodeID<<"____"<<endl;
            double bufd[3];
            TColStd_Array1OfReal coords(*bufd,1,3);
            int NbNodes1;
            MeshVS_EntityType eType1;
            inputMesh->GetGeom(globalNodeID,false,coords,NbNodes1,eType1);
            mesh::meshPoint aP(coords(1),coords(2),coords(3));
            aP.ID = n;
            ameshElementByCoords.pointList<<aP;
        }
        //cout<<"____point list size: "<<pointList.size()<<"____"<<endl;
        //listOfElements<<ameshElementByCoords;
        listOfElements.push_back(ameshElementByCoords);
    }

    //! ------------------
    //! mesh construction
    //! ------------------
    outputMesh = new Ng_MeshVS_DataSource3D(listOfElements);
}

//! -----------------------------
//! function: toListOf3DElements
//! details:
//! -----------------------------
void MeshTools::toListOf3DElements(const occHandle(Ng_MeshVS_DataSource3D) &inputMesh, std::vector<meshElementByCoords> &elements)
{
    if(inputMesh.IsNull()) return;
    if(inputMesh->GetAllElements().Extent()<1) return;
    if(inputMesh->GetAllNodes().Extent()<1) return;

    for(TColStd_MapIteratorOfPackedMapOfInteger it(inputMesh->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        int buf[20], NbNodes;
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        inputMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        ElemType anElemType;
        switch(NbNodes)
        {
        case 4: anElemType = TET; break;
        case 5: anElemType = PYRAM; break;
        case 6: anElemType = PRISM; break;
        case 8: anElemType = HEXA; break;
        }

        meshElementByCoords ameshElementByCoords;
        ameshElementByCoords.type = anElemType;

        for(int n=1; n<=NbNodes; n++)
        {
            int globalNodeID = nodeIDs(n);
            double bufd[3];
            TColStd_Array1OfReal coords(*bufd,1,3);
            int NbNodes1;
            MeshVS_EntityType eType1;
            inputMesh->GetGeom(globalNodeID,false,coords,NbNodes1,eType1);
            mesh::meshPoint aP(coords(1),coords(2),coords(3));
            aP.ID = n;
            ameshElementByCoords.pointList<<aP;
        }
        elements.push_back(ameshElementByCoords);
        //elements<<ameshElementByCoords;
    }
}


//! ----------------------
//! function: groupMeshes
//! details:
//! ----------------------
std::map<GeometryTag, std::vector<occHandle(MeshVS_Mesh)>> MeshTools::groupMeshes(const std::map<GeometryTag, occHandle(MeshVS_Mesh)> &mapOfMeshes)
{
    cout<<"MeshTools::groupMeshes()->____grouping meshes____"<<endl;

    std::map<GeometryTag,std::vector<occHandle(MeshVS_Mesh)>> body2MeshVector;
    for(std::map<GeometryTag,occHandle(MeshVS_Mesh)>::const_iterator it = mapOfMeshes.cbegin(); it!=mapOfMeshes.cend(); it++)
    {
        const GeometryTag &aTag = it->first;
        //const occHandle(MeshVS_Mesh) &aMeshDS = it->second;

        //! ----------------------------------------
        //! a tag representing a main shape (SOLID)
        //! ----------------------------------------
        int bodyIndex = aTag.parentShapeNr;
        GeometryTag bodyTag(bodyIndex,bodyIndex,true,TopAbs_SOLID);
        std::map<GeometryTag,std::vector<occHandle(MeshVS_Mesh)>>::iterator it_ = body2MeshVector.find(bodyTag);
        if(it_==body2MeshVector.end())
        {
            std::vector<occHandle(MeshVS_Mesh)> vecMeshes { it->second };
            body2MeshVector.insert(std::make_pair(bodyTag,vecMeshes));
        }
        else it_->second.push_back(it->second);
    }
    return body2MeshVector;
}

//! --------------------
//! function: addMeshes
//! details:
//! --------------------
occHandle(MeshVS_DataSource) MeshTools::mergeMesh(const occHandle(MeshVS_DataSource) &mesh1, const occHandle(MeshVS_DataSource) &mesh2)
{
    //cout<<"MeshTools::mergeMesh()->____function called____"<<endl;

    occHandle(MeshVS_DataSource) retMeshDS;

    if(mesh1.IsNull()) return mesh2;
    if(mesh2.IsNull()) return mesh1;
    if(mesh2.IsNull() && mesh1.IsNull()) return retMeshDS;

    int meshDim = 0;    //! 0,1,2,3 D

    std::vector<meshElementByCoords> vecElements;
    std::vector<occHandle(MeshVS_DataSource)> vecMeshes {mesh1, mesh2};
    for(int i=0; i<=1; i++)
    {
        const occHandle(MeshVS_DataSource) &curMeshDS = vecMeshes[i];
        for(TColStd_MapIteratorOfPackedMapOfInteger it = curMeshDS->GetAllElements(); it.More(); it.Next())
        {
            int globalElementID = it.Key();

            meshElementByCoords aMeshElement;
            aMeshElement.ID = globalElementID;

            int NbNodes, buf[20];
            TColStd_Array1OfInteger nodeIDs(*buf,1,20);
            curMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

            MeshVS_EntityType aType;
            curMeshDS->GetGeomType(globalElementID,true,aType);

            ElemType anElementType;
            switch(aType)
            {
            case MeshVS_ET_Volume:
            {
                meshDim = 3;
                switch(NbNodes)
                {
                case 4: anElementType = TET; break;
                case 5: anElementType = PYRAM ; break;
                case 6: anElementType = PRISM; break;
                case 8: anElementType = HEXA; break;
                case 10: anElementType = TET10; break;
                case 13: anElementType = PYRAM13 ; break;
                case 15: anElementType = PRISM15; break;
                case 20: anElementType = HEXA20; break;
                }
            }
                break;

            case MeshVS_ET_Face:
            {
                meshDim = 2;
                switch(NbNodes)
                {
                case 3: anElementType = TRIG; break;
                case 4: anElementType = QUAD; break;
                case 6: anElementType = TRIG6; break;
                case 8: anElementType = QUAD8; break;
                }
            }
                break;

            case MeshVS_ET_Link:
            {
                meshDim = 1;
                // to do
            }
                break;
            }

            aMeshElement.type = anElementType;

            for(int i=1; i<=NbNodes; i++)
            {
                int globalNodeID = nodeIDs(i);
                double dbuf[3];
                TColStd_Array1OfReal coords(*dbuf,1,3);
                int NbNodes1;
                curMeshDS->GetGeom(globalNodeID,false,coords,NbNodes1,aType);
                mesh::meshPoint aP(coords(1),coords(2),coords(3),globalNodeID);
                aMeshElement.pointList<<aP;
            }

            //! pile up the elements
            vecElements.push_back(aMeshElement);
        }
    }

    switch(meshDim)
    {
    case 0: break;
    case 1: break;
    case 2: retMeshDS = new Ng_MeshVS_DataSourceFace(vecElements,false,false); break;
    case 3: retMeshDS = new Ng_MeshVS_DataSource3D(vecElements,false,false); break;
    }

    //cout<<"MeshTools::mergeMesh()->____number of nodes: "<<retMeshDS->GetAllNodes().Extent()<<"____"<<endl;
    //cout<<"MeshTools::mergeMesh()->____number of elements: "<<retMeshDS->GetAllElements().Extent()<<"____"<<endl;
    //cout<<"MeshTools::mergeMesh()->____function exiting____"<<endl;

    return retMeshDS;
}


//! --------------------------------
//! function: computeAngleDefectMap
//! details:
//! --------------------------------
void MeshTools::computeAngleDefectMap(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, std::map<int,double> &mapOfAngleDefect)
{
    cout<<"MeshTools::computeAngleDefectMap()->____function called____"<<endl;

    const double PI = 3.1415926538;

    if(aMeshDS.IsNull()) return;
    if(aMeshDS->myNodeToElements.isEmpty()) aMeshDS->computeNodeToElementsConnectivity();   // nodes and elements local numbering
    if(aMeshDS->myNodeNormals.isEmpty()) aMeshDS->computeNormalAtNodes();
    if(aMeshDS->myBoundaryPoints.isEmpty()) aMeshDS->computeFreeMeshSegments();
    cout<<"____number of boundary points: "<<aMeshDS->myBoundaryPoints.length()<<"____"<<endl;

    //! --------------------
    //! scan the mesh nodes
    //! --------------------
    TColStd_MapIteratorOfPackedMapOfInteger it(aMeshDS->GetAllNodes());
    int NbMeshNodes = aMeshDS->GetAllNodes().Extent();
    for(int localNodeID = 1; localNodeID<=NbMeshNodes; localNodeID++, it.Next())
    {
        int globalNodeID = it.Key();

        cout<<" ***************************************"<<endl;
        cout<<"  computing angle deficit @ nodeID: "<<globalNodeID<<endl;

        const QList<int> &attachedElements_local = aMeshDS->myNodeToElements.value(localNodeID);
        int NbAttachedElements = attachedElements_local.length();

        cout<<"  attached elements: "<<NbAttachedElements<<endl;
        cout<<" ***************************************"<<endl;

        std::set<int> surroundingNodes_global;
        for(int n = 0; n<NbAttachedElements; n++)
        {
            int localElementID_attached = attachedElements_local[n];
            int globalElementID_attached = aMeshDS->myElementsMap.FindKey(localElementID_attached);
            int buf[20], NbElementNodes;
            TColStd_Array1OfInteger nodeIDs(*buf,1,20);
            aMeshDS->GetNodesByElement(globalElementID_attached,nodeIDs,NbElementNodes);
            for(int i=1; i<=NbElementNodes; i++)
                surroundingNodes_global.insert(nodeIDs(i));
        }
        for(std::set<int>::iterator it = surroundingNodes_global.begin(); it !=surroundingNodes_global.end(); it++)
            cout<<"*____node ID: "<<*it<<"____"<<endl;

        //! -----------------------------------------------
        //! discard the "globalNodeID" - center of the fan
        //! -----------------------------------------------
        std::set<int>::iterator it__ = surroundingNodes_global.find(globalNodeID);
        if(it__==surroundingNodes_global.end()) exit(20);       // questo non deve mai capitare
        surroundingNodes_global.erase(it__);


        for(std::set<int>::iterator it = surroundingNodes_global.begin(); it !=surroundingNodes_global.end(); it++)
            cout<<"*____node ID: "<<*it<<"____"<<endl;


        //! ---------------------------
        //! center of the triangle fan
        //! ---------------------------
        int NbNodes_;
        double bufd[3];
        TColStd_Array1OfReal coords(*bufd,1,3);
        MeshVS_EntityType aType;
        aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes_,aType);     // center of the triangle fan
        double xP = coords(1);
        double yP = coords(2);
        double zP = coords(3);

        //! ------------------------------
        //! convert the set into a vector
        //! ------------------------------
        std::vector<int> vecOfSurroundingNodes;
        for(std::set<int>::iterator it_ = surroundingNodes_global.begin(); it_ != surroundingNodes_global.end(); it_++)
        {
            cout<<"____surrounding nodes: "<<*it_<<"____"<<endl;
            vecOfSurroundingNodes.push_back(*it_);
        }

        //! ---------------------------------------------------
        //! first point of the vector (first surrounding node)
        //! ---------------------------------------------------
        aMeshDS->GetGeom(vecOfSurroundingNodes[0],false,coords,NbNodes_,aType);
        double xP0 = coords(1);
        double yP0 = coords(2);
        double zP0 = coords(3);

        double vecstartx = xP0-xP;
        double vecstarty = yP0-yP;
        double vecstartz = zP0-zP;

        //! -----------------------------------------------------------
        //! vector for the orientation: the normal at the current node
        //! -----------------------------------------------------------
        QList<double> nodeNormal;
        aMeshDS->myNodeNormals.value(globalNodeID,nodeNormal);
        double nx = nodeNormal[0];
        double ny = nodeNormal[1];
        double nz = nodeNormal[2];

        std::map<double,int> angleMap;
        angleMap.insert(std::make_pair(0.0,vecOfSurroundingNodes[0]));

        double lP0P = sqrt(pow((xP0-xP),2)+pow((yP0-yP),2)+pow((zP0-zP),2));
        size_t Nb = vecOfSurroundingNodes.size();
        for(int i = 0; i<Nb; i++)
        {
            aMeshDS->GetGeom(vecOfSurroundingNodes[i],false,coords,NbNodes_,aType);
            double xPi = coords(1);
            double yPi = coords(2);
            double zPi = coords(3);

            double lP1P = sqrt(pow((xPi-xP),2)+pow((yPi-yP),2)+pow((zPi-zP),2));
            double dot = (xP0-xP)*(xPi-xP)+(yP0-yP)*(yPi-yP)+(zP0-zP)*(zPi-zP);
            dot /= lP0P*lP1P;

            if(dot<-1) dot = -1;
            if(dot>1) dot = 1;
            double angle = std::acos(dot);

            //! ----------------------------------
            //! i           j           k
            //! vecstartx   vecstarty   vecstartz
            //! vecx        vecy        vecz
            //! ----------------------------------
            double vecx = xPi-xP;
            double vecy = yPi-yP;
            double vecz = zPi-zP;

            double crossx = vecstarty*vecz-vecstartz*vecy;
            double crossy = vecstartz*vecx-vecstartx*vecz;
            double crossz = vecstartx*vecy-vecstarty*vecx;

            double flag = crossx*nx+crossy*ny+crossz*nz;
            if(flag < 0) angle = 2*PI - angle;
            angleMap.insert(std::make_pair(angle,vecOfSurroundingNodes[i]));
        }

        //std::vector<int> reordererArray { globalNodeID };
        std::vector<int> reordererArray;
        for(std::map<double,int>::iterator it = angleMap.begin(); it!=angleMap.end(); it++)
        {
            int globalNodeID_ = it->second;
            double angle = it->first;
            reordererArray.push_back(globalNodeID_);
            cout<<"____inserting nodes: "<<globalNodeID_<<" angle: "<<angle*180/PI<<" into the support map____"<<endl;
        }

        size_t N = reordererArray.size();

        for(int i=0; i<N; i++) cout<<"____reordered array: "<<reordererArray[i]<<"____"<<endl;

        double angleDefect = 0;
        for(int i=0; i<N; i++)
        {
            int i1 = i%N;
            aMeshDS->GetGeom(vecOfSurroundingNodes[i1],false,coords,NbNodes_,aType);
            double xi = coords(1);
            double yi = coords(2);
            double zi = coords(3);

            int i2 = (i+1)%N;
            aMeshDS->GetGeom(vecOfSurroundingNodes[i2],false,coords,NbNodes_,aType);
            double xii = coords(1);
            double yii = coords(2);
            double zii = coords(3);

            double dot = (xii-xP)*(xi-xP)+(yii-yP)*(yi-yP)+(zii-zP)*(zi-zP);
            double li = sqrt(pow((xi-xP),2)+pow((yi-yP),2)+pow((zi-zP),2));
            double lii = sqrt(pow((xii-xP),2)+pow((yii-yP),2)+pow((zii-zP),2));
            dot /= li*lii;
            if(dot<-1) dot = -1;
            if(dot>1) dot = 1;

            angleDefect += std::acos(dot);
        }
        angleDefect = 2*PI-angleDefect;
        if(aMeshDS->myBoundaryPoints.contains(globalNodeID)) angleDefect = PI-angleDefect;
        cout<<"____(nodeID, angleDeficit) = ("<<globalNodeID<<", "<<angleDefect*180/PI<<")____"<<endl;
        mapOfAngleDefect.insert(std::make_pair(globalNodeID,angleDefect));
    }
}


//! ------------------
//! function: saveSTL
//! details:
//! ------------------
bool MeshTools::saveSTL(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, const std::string &fileName)
{
    if(aMeshDS.IsNull()) return false;
    if(fileName=="") return false;

    FILE *aFile = fopen(fileName.c_str(),"w");
    if(aFile == NULL) return false;

    fprintf(aFile,"solid surfacemesh\n");
    for(TColStd_MapIteratorOfPackedMapOfInteger it = aMeshDS->GetAllElements(); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        MeshVS_EntityType aType;
        double buf[24];
        TColStd_Array1OfReal coords(*buf,1,24);
        int NbNodes;
        double nx, ny,nz;
        aMeshDS->GetNormal(globalElementID,10,nx,ny,nz);
        if(aMeshDS->GetGeom(globalElementID,true,coords,NbNodes,aType))
        {
            fprintf(aFile,
                    "facet normal %.6e %.6e %.6e\n"
                    "outer loop\n"
                    "vertex %.6e %.6e %.6e\n"
                    "vertex %.6e %.6e %.6e\n"
                    "vertex %.6e %.6e %.6e\n"
                    "endloop\n"
                    "endfacet\n",
                    nx,ny,nz,
                    coords(1),coords(2),coords(3),
                    coords(4),coords(5),coords(6),
                    coords(7),coords(8),coords(9)
                    );
        }
        else
        {
            cerr<<"____the node cannot be found: the STL is not valid____"<<endl;
            fclose(aFile);
            return false;
        }
    }
    fprintf(aFile,"endsolid surfacemesh\n");
    fclose(aFile);
    return true;
}
