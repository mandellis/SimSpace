//! ----------------
//! custom includes
//! ----------------
#include "meshselectnlayers.h"
#include <mesh.h>
#include <meshelement2d.h>
#include <meshelementbycoords.h>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>

//! ----
//! C++
//! ----
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
meshSelectNLayers::meshSelectNLayers(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS, int N):
    myVolumeMeshDS(aVolumeMeshDS),
    myN(N)
{
    ;
}

//! -------------------
//! function: setSteps
//! details:
//! -------------------
void meshSelectNLayers::setSteps(int NbSteps)
{
    myNbSteps = NbSteps;
}

//! ------------------------
//! function: setVolumeMesh
//! details:
//! ------------------------
void meshSelectNLayers::setVolumeMesh(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMesh)
{
    myVolumeMeshDS = aVolumeMesh;
}

//! -------------------------------------
//! function: setStartingSurfaceElements
//! details:
//! -------------------------------------
void meshSelectNLayers::setStartingSurfaceElements(const opencascade::handle<Ng_MeshVS_DataSourceFace> &facesDS)
{
    myStartFace = facesDS;
}

//! ----------------------
//! function: getElements
//! details:
//! ----------------------
bool meshSelectNLayers::getElements(occHandle(Ng_MeshVS_DataSource3D) &theElements)
{
    //! -------------
    //! sanity check
    //! -------------
    if(myVolumeMeshDS.IsNull()) return false;
    if(myStartFace.IsNull()) return false;

    //! ---------------------------------------------------------
    //! build the volume element to volume elements connectivity
    //! ---------------------------------------------------------
    std::map<int,std::vector<int>> elementsAttachedToElement;
    myVolumeMeshDS->buildElementToElementConnectivity(elementsAttachedToElement);
    if(elementsAttachedToElement.size()==0)
    {
        cerr<<"meshSelectNLayers::getElements()->____something very strange has happened when building element to element connectivity____"<<endl;
        return false;
    }

    //! ----------------------------------------------
    //! the volume elements building the n-step layer
    //! ----------------------------------------------
    std::vector<int> vecGlobalElementIDs;

    //! ---------------------------------------------------------------
    //! first step: iterate over the surface elements of "myStartFace"
    //! ---------------------------------------------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myStartFace->GetAllElements()); it.More(); it.Next())
    {
        int globalSurfaceElementID = it.Key();
        int bufn[8];
        TColStd_Array1OfInteger nodeIDs(*bufn,1,8);
        int NbNodes;
        myStartFace->GetNodesByElement(globalSurfaceElementID,nodeIDs,NbNodes);

        meshElement2D aMeshElement2D;
        for(int i=1; i<=NbNodes; i++) aMeshElement2D.nodeIDs<<nodeIDs(i);

        //! --------------------
        //! this is unessential
        //! --------------------
        switch(NbNodes)
        {
        case 3: aMeshElement2D.type = TRIG; break;
        case 4: aMeshElement2D.type = QUAD; break;
        }

        //! -------------------------------------------------------------------
        //! search the meshElement2D into the face to element connectivity map
        //! at this first iteration only one volume element corresponds to a
        //! face mesh element (because they have been taken on the surface by
        //! definition).
        //! -------------------------------------------------------------------
        QList<int> volumeElementsAttachedToSurfaceElement = myVolumeMeshDS->myFaceToElements.value(aMeshElement2D);
        for(int i=0; i<volumeElementsAttachedToSurfaceElement.length(); i++)
        {
            int globalElementID = volumeElementsAttachedToSurfaceElement.at(i);
            //! ------------------------------------------------------
            //! check if the volume element has been already recorded
            //! ------------------------------------------------------
            if(std::find(vecGlobalElementIDs.begin(),vecGlobalElementIDs.end(),globalElementID)==vecGlobalElementIDs.end())
                vecGlobalElementIDs.push_back(globalElementID);
        }
    }

    //! --------------------
    //! the remaining steps
    //! --------------------
    for(int n = 0; n<myNbSteps; n++)
    {
        std::vector<int> elementsToScan = vecGlobalElementIDs;

        //! -----------------------------
        //! add elements while iterating
        //! -----------------------------
        for(int k=0; k<elementsToScan.size(); k++)
        {
            int globalElementID = elementsToScan.at(k);

            //! -------------------------------------------------------------------------------
            //! the current element is not contained into the element to attached elements map
            //! -------------------------------------------------------------------------------
            std::map<int,std::vector<int>>::iterator it1 = elementsAttachedToElement.find(globalElementID);
            if(it1==elementsAttachedToElement.end()) continue;

            //! ---------------------------------------------------------------------------
            //! the current element is contained into the element to attached elements map
            //! ---------------------------------------------------------------------------
            std::pair<int,std::vector<int>> aPair = *it1;
            std::vector<int> attachedVolumeElements = aPair.second;
            for(int h=0; h<attachedVolumeElements.size(); h++)
            {
                int attachedGlobalElementID = attachedVolumeElements.at(h);
                if(std::find(vecGlobalElementIDs.begin(), vecGlobalElementIDs.end(), attachedGlobalElementID) == vecGlobalElementIDs.end())
                    vecGlobalElementIDs.push_back(attachedGlobalElementID);
            }
        }
    }

    if(vecGlobalElementIDs.size()==0) return false;

    //! ---------------------------------------------------------------
    //! build the resulting mesh - now only for visualization purposes
    //! ---------------------------------------------------------------
    cout<<"meshSelectNLayers::getElements()->____calling volume mesh constructor____"<<endl;
    QList<meshElementByCoords> listVolumeElements;
    for(int n=0; n<vecGlobalElementIDs.size(); n++)
    {
        int globalElementID = vecGlobalElementIDs[n];
        int buf[20], NbNodes;
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        bool isDone = myVolumeMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        if(!isDone) continue;

        meshElementByCoords aVolumeElement;
        aVolumeElement.ID = globalElementID;
        switch(NbNodes)
        {
        case 4: aVolumeElement.type = TET; break;
        case 5: aVolumeElement.type = PYRAM; break;
        case 6: aVolumeElement.type = PRISM; break;
        case 8: aVolumeElement.type = HEXA; break;
        }

        for(int n=1; n<=NbNodes; n++)
        {
            int globalNodeID = nodeIDs(n);
            int localNodeID = myVolumeMeshDS->myNodesMap.FindIndex(globalNodeID);
            const std::vector<double> &P = myVolumeMeshDS->getNodeCoordinates(localNodeID);
            mesh::meshPoint aP(P[0],P[1],P[2],globalNodeID);
            aVolumeElement.pointList<<aP;
        }
        listVolumeElements<<aVolumeElement;
    }

    //! ---------------
    //! the final mesh
    //! ---------------
    //theElements = new Ng_MeshVS_DataSource3D(listVolumeElements);
    theElements = new Ng_MeshVS_DataSource3D(listVolumeElements,false,false);
    return true;
}
