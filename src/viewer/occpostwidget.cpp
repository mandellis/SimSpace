//! ----------------
//! custom includes
//! ----------------
#include <occpostwidget.h>
#include "graphicstools.h"
#include "meshtools.h"
#include "simulationdatabase.h"
#include "simulationmanager.h"
#include "resultstoolbar.h"
#include <meshslicer.h>

//! ----
//! OCC
//! ----
#include <NCollection_TListIterator.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DisplayModeFlags.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include <MeshVS_MeshPrsBuilder.hxx>

//! ----------------------
//! function: constructor
//! detail:
//! ----------------------
occPostWidget::occPostWidget(QWidget *parent):occPreGLWidget(parent)
{
    cout<<"occPostWidget::occPostWidget()->____CONSTRUCTOR CALLED I____"<<endl;

    //! -------------------------------------
    //! the occ context for the results view
    //! -------------------------------------
    //occPostContext = new AIS_InteractiveContext(occViewer);
}

//! ----------------------
//! function: constructor
//! detail:
//! ----------------------
occPostWidget::occPostWidget(meshDataBase *mDB, QWidget *parent):occPreGLWidget(parent),
    myMeshDataBase(mDB)
{
    cout<<"occPostWidget::occPostWidget()->____CONSTRUCTOR CALLED II____"<<endl;

    //! -------------------------------------
    //! the occ context for the results view
    //! -------------------------------------
    //occPostContext = new AIS_InteractiveContext(occViewer);

    //! ---------------------------------
    //! the result presentation settings
    //! ---------------------------------
    myResultPresentation = Global::status().myResultPresentation;
}

//! ---------------
//! function: init
//! details:
//! ---------------
void occPostWidget::init()
{
    occPreGLWidget::init();
}

//! ---------------------
//! function: paintEvent
//! details:
//! ---------------------
void occPostWidget::paintEvent(QPaintEvent *e)
{
    Q_UNUSED(e);
    occPreGLWidget::paintEvent(e);
    if(occPostContext.IsNull()) occPostContext = new AIS_InteractiveContext(occViewer);
}

//! -------------------------
//! function: createColorBox
//! detail:
//! -------------------------
void occPostWidget::createColorBox(double min, double max, int Nintervals)
{
    cout<<"occPostWidget::createColorBox->____function called____"<<endl;

    //! ---------------------------------------------------------
    //! https://tracker.dev.opencascade.org/view.php?id=25785
    //!
    //! GB: setColor sets the color of the text and of the frame
    //! the code of the class has been copied from OCC source
    //! since it was not possible to change the width of the box
    //! ---------------------------------------------------------
    graphicsTools::createColorBox(min,max,Nintervals,myColorBox);
    occHandle(Graphic3d_TransformPers) trs = new Graphic3d_TransformPers(Graphic3d_TMF_2d);
    trs->SetCorner2d(Aspect_TOTP_LEFT_UPPER);

    Graphic3d_Vec2i offset(20,425);

    trs->SetOffset2d(offset);
    myColorBox->SetZLayer(Graphic3d_ZLayerId_TopOSD);
    myColorBox->SetTransformPersistence(trs);
    myColorBox->SetToUpdate();
}

//! --------------------------
//! function: displayColorBox
//! detail:
//! --------------------------
void occPostWidget::displayColorBox(bool isVisible)
{
    if(!myColorBox.IsNull())
    {
        if(isVisible)
        {
            occPostContext->Display(myColorBox,1);
            occViewer->Update();
        }
        else
        {
            occPostContext->Remove(myColorBox,true);
        }
    }
}

//! ------------------------
//! function: displayResult
//! details:
//! ------------------------
void occPostWidget::displayResult(sharedPostObject &aPostObject)
{
    cout<<"occPostWidget::displayResult()->____function called____"<<endl;

    static resultPresentation aResultPresentationOld;

    //! ----------------
    //! read the status
    //! ----------------
    resultPresentation aResultPresentation = myResultPresentation;

    //! ----------------------------------------
    //! delete all the objects from the display
    //! ----------------------------------------
    occMeshContext->RemoveAll(false);
    occPostContext->RemoveAll(false);

    //! ---------------------------
    //! display the colored meshes
    //! ---------------------------
    const std::map<GeometryTag,occHandle(MeshVS_Mesh)> &coloredMeshes = aPostObject->getColoredMeshes();
    for(std::map<GeometryTag,occHandle(MeshVS_Mesh)>::const_iterator it = coloredMeshes.cbegin(); it != coloredMeshes.cend(); ++it)
    {
        const occHandle(MeshVS_Mesh) &aColoredMesh = it->second;
        occPostContext->Display(aColoredMesh,false);
    }

    //! ----------------------------------------------
    //! handle the combined view of results and model
    //! ----------------------------------------------
    switch(aResultPresentation.theCombinedView)
    {
    case resultPresentation::combinedView_resultOnly:
    {
        ;
    }
        break;

    case resultPresentation::combinedView_meshVisible:
    {
        const std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)> &mapOfMeshDS = aPostObject->getMeshDataSourcesForView();
        for(std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)>::const_iterator it = mapOfMeshDS.cbegin(); it!=mapOfMeshDS.cend(); it++)
        {
            occHandle(MeshVS_Mesh) aMeshObject = new MeshVS_Mesh();
            aMeshObject->SetDataSource(it->second);

            aMeshObject->SetDisplayMode(MeshVS_DMF_WireFrame);
            aMeshObject->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges,false);
            aMeshObject->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
            aMeshObject->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes,false);
            occHandle(MeshVS_MeshPrsBuilder) aB = new MeshVS_MeshPrsBuilder(aMeshObject);
            aMeshObject->AddBuilder(aB,false);

            occMeshContext->Display(aMeshObject,false);
        }
    }
        break;

    case resultPresentation::combinedView_undeformedWireFrame:
    {
        const std::vector<GeometryTag> &vecLocs = aPostObject->getLocations();
        for(std::vector<GeometryTag>::const_iterator it = vecLocs.cbegin(); it!=vecLocs.cend(); it++)
        {
            int bodyIndex = it->parentShapeNr;
            const TopoDS_Shape &aShape = myDS2->bodyMap.value(bodyIndex);
            const occHandle(AIS_Shape) &anAISShape = new AIS_Shape(aShape);
            occPostContext->SetColor(anAISShape,Quantity_NOC_BLACK,false);
            occPostContext->Display(anAISShape,AIS_WireFrame,-1,false,false);
        }
    }
        break;

    case resultPresentation::combinedView_undeformedModel:
    {
        const std::vector<GeometryTag> &vecLocs = aPostObject->getLocations();
        for(std::vector<GeometryTag>::const_iterator it = vecLocs.cbegin(); it!=vecLocs.cend(); it++)
        {
            int bodyIndex = it->parentShapeNr;
            const TopoDS_Shape &aShape = myDS2->bodyMap.value(bodyIndex);
            const occHandle(AIS_Shape) &anAISShape = new AIS_Shape(aShape);
            occPostContext->SetColor(anAISShape,Quantity_NOC_AQUAMARINE1,false);
            occPostContext->SetTransparency(anAISShape,0.9,true);
            occPostContext->Display(anAISShape,AIS_Shaded,-1,false,false);
        }
    }
       break;
    }

    //! ----------------------
    //! display the color box
    //! ----------------------
    const occHandle(AIS_ColorScaleExtended) &colorBox = aPostObject->getColorBox();
    occPostContext->Display(colorBox,false);

    //! ------------------
    //! update the viewer
    //! ------------------
    occPostContext->UpdateCurrentViewer();
    occMeshContext->UpdateCurrentViewer();
    occContext->UpdateCurrentViewer();

    aResultPresentationOld = myResultPresentation;
}

//! ---------------------
//! function: clipResult
//! details:
//! ---------------------
void occPostWidget::clipResult()
{
    cout<<"occPostWidget::clipResult()->____function called____"<<endl;

    //! --------------
    //! a mesh slicer
    //! --------------
    meshSlicer aSlicer;

    //! ------------------------------
    //! iterate over the post objects
    //! ------------------------------
    AIS_ListOfInteractive listOfIO;
    occPostContext->ObjectsInside(listOfIO);
    for(AIS_ListIteratorOfListOfInteractive it(listOfIO); it.More(); it.Next())
    {
        const occHandle(MeshVS_Mesh) &aMesh = occHandle(MeshVS_Mesh)::DownCast(it.Value());
        if(aMesh.IsNull()) continue;

        //! --------------------------
        //! reset the hidden elements
        //! --------------------------
        aMesh->SetHiddenNodes(new TColStd_HPackedMapOfInteger());

        const occHandle(MeshVS_DataSource) &aMeshDS = aMesh->GetDataSource();
        if(aMeshDS.IsNull()) continue;
        aSlicer.setMeshDataSource(aMeshDS);

        TColStd_PackedMapOfInteger hiddenElementIDs;
        for(QMap<int,occHandle(Graphic3d_ClipPlane)>::const_iterator itplane = myMapOfClipPlanes.cbegin(); itplane != myMapOfClipPlanes.cend(); itplane++)
        {
            const occHandle(Graphic3d_ClipPlane) &aClipPlane = itplane.value();
            if(aClipPlane->IsOn()==false) continue;
            Graphic3d_ClipPlane::Equation eq = aClipPlane->GetEquation();
            double a = eq.GetData()[0];
            double b = eq.GetData()[1];
            double c = eq.GetData()[2];
            double d = eq.GetData()[3];

            occHandle(TColStd_HPackedMapOfInteger) HHiddenElementIDs;
            bool isDone = aSlicer.perform(a,b,c,d,HHiddenElementIDs);
            if(isDone == false) return;
            hiddenElementIDs.Unite(HHiddenElementIDs->Map());
        }
        occHandle(TColStd_HPackedMapOfInteger) mapOfHiddenElements = new TColStd_HPackedMapOfInteger;
        mapOfHiddenElements->ChangeMap() = hiddenElementIDs;
        aMesh->SetHiddenElems(mapOfHiddenElements);
        occPostContext->RecomputePrsOnly(aMesh,false,false);
    }
    occPostContext->UpdateCurrentViewer();
}

//! -------------------------------------
//! function: updateViewerStatus
//! details:  update the status variable
//! -------------------------------------
void occPostWidget::updateViewerStatus()
{
    myResultPresentation = Global::status().myResultPresentation;

    //switch(myResultPresentation.theCombinedView)
    //{
    //case resultPresentation::combinedView_resultOnly: cout<<"occPostWidget::updateViewerStatus()->____results only____"<<endl; break;
    //case resultPresentation::combinedView_meshVisible: cout<<"occPostWidget::updateViewerStatus()->____results with mesh____"<<endl; break;
    //case resultPresentation::combinedView_undeformedModel: cout<<"occPostWidget::updateViewerStatus()->____results and undeformed model____"<<endl; break;
    //case resultPresentation::combinedView_undeformedWireFrame: cout<<"occPostWidget::updateViewerStatus()->____results and undeformed wireframe____"<<endl; break;
    //}
    //cout<<"occPostWidget::updateViewerStatus()->____scale: "<<myResultPresentation.theScale<<"____"<<endl;

    //! --------------------------------------------------------------
    //! this method calls the simulation manager, which in turn calls
    //! <this>::displayResult()
    //! --------------------------------------------------------------
    emit resultsPresentationChanged();
}

//! -------------------------
//! function: hideAllResults
//! details:
//! -------------------------
void occPostWidget::hideAllResults()
{
    //cout<<"occPostWidget::hideAllResults()->____function called____"<<endl;
    occPostContext->RemoveAll(true);
}

//! -----------------------
//! function: updateResult
//! details:
//! -----------------------
void occPostWidget::updateResult(postObject &aPostObject)
{
    aPostObject.init(static_cast<meshDataBase*>(myDS2));
}

//! ----------------------------------
//! function: setWorkingMode_Solution
//! details:
//! ----------------------------------
void occPostWidget::setWorkingMode_Solution()
{
    if(myCurWorkingMode!=curWorkingMode_onSolution)
    {
        cout<<"occPreGLWidget::setWorkingMode_Solution()->____ON SOLUTION____"<<endl;

        //! ---------------------------
        //! change the internal status
        //! ---------------------------
        myCurWorkingMode = curWorkingMode_onSolution;

        //! -------------------------------------------------------
        //! display the bodies in wireframe mode (this allow to
        //! select the center of rotation) and hide all the meshes
        //! -------------------------------------------------------
        this->setWireframeView();
        this->hideAllBodies();
        this->hideAllMeshes();

        //! ------------------------------
        //! disable the selection buttons
        //! ------------------------------
        emit requestDisableSelectionButtons();

        //! --------------------------------------
        //! deactivate the current selection mode
        //! --------------------------------------
        this->clearGeometrySelection();

        AIS_ListOfInteractive listOfIO;
        for(AIS_ListIteratorOfListOfInteractive it(listOfIO); it.More(); it.Next())
        {
            this->setSelectionMode(CurSelection_Nothing);
        }

        //switch(myCurSelectionMode)
        //{
        //case CurSelection_Solid: occContext->DeactivateStandardMode(TopAbs_SOLID); break;
        //case CurSelection_Face: occContext->DeactivateStandardMode(TopAbs_FACE); break;
        //case CurSelection_Edge: occContext->DeactivateStandardMode(TopAbs_EDGE); break;
        //case CurSelection_Vertex: occContext->DeactivateStandardMode(TopAbs_VERTEX); break;
        //default: break;
        //}

        //! --------------------------
        //! show the results tool bar
        //! --------------------------
        ResultsToolBar *theResultsToolBar = static_cast<ResultsToolBar*>(tools::getWidgetByName("resultsToolBar"));
        theResultsToolBar->setVisible(true);

        theResultsToolBar->blockSignals(true);
        //theResultsToolBar->setStatus(myResultPresentation);
        theResultsToolBar->blockSignals(false);
    }
}

//! --------------------------
//! function: refreshMeshView
//! details:
//! --------------------------
void occPostWidget::refreshMeshView(bool onlyExterior)
{
    occPreGLWidget::refreshMeshView(onlyExterior);
    //this->clipResult();
    isMeshViewVolume = (onlyExterior == true? false:true);
}

//! ----------------
//! function: reset
//! details:
//! ----------------
void occPostWidget::reset()
{
    occPreGLWidget::reset();
}
