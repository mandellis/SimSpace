//! ----------------
//! custom includes
//! ----------------
#include <occpostwidget.h>
#include "graphicstools.h"
#include "meshtools.h"
#include "simulationdatabase.h"
#include "simulationmanager.h"
#include "resultstoolbar.h"

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
    myResultPresentation.theCombinedView = resultPresentation::combinedView_resultOnly;
    myResultPresentation.isDeformedView = false;
    myResultPresentation.theScale = 1.0;
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
        cout<<"occPostWidget::displayColorBox(isVisible)->____function called____"<<endl;
        if(isVisible)
        {
            //occContext->Display(myColorBox,1);
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
//void occPostWidget::displayResult(postObject &aPostObject)
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
    occMeshContext->RemoveAll(true);
    occPostContext->RemoveAll(true);
    occContext->RemoveAll(true);

    //! -----------------
    //! update the scale
    //! -----------------
    if(aResultPresentation.theScale != aResultPresentationOld.theScale)
    {
        cout<<"@ --------------------------------"<<endl;
        cout<<"@ - SCALE CHANGED"<<endl;
        cout<<"@ --------------------------------"<<endl;
        aPostObject->setScale(aResultPresentation.theScale);
        aPostObject->updateScaledView();
    }

    //! ---------------------------
    //! display the colored meshes
    //! ---------------------------
    const QMap<GeometryTag,occHandle(MeshVS_Mesh)> &coloredMeshes = aPostObject->getColoredMeshes();
    for(QMap<GeometryTag,occHandle(MeshVS_Mesh)>::const_iterator it = coloredMeshes.cbegin(); it != coloredMeshes.cend(); ++it)
    {
        const occHandle(MeshVS_Mesh) &aColoredMesh = it.value();
        occPostContext->Display(aColoredMesh,false);
    }

    //! ----------------------------------------------
    //! handle the combined view of results and model
    //! ----------------------------------------------
    switch(aResultPresentation.theCombinedView)
    {
    case resultPresentation::combinedView_resultOnly:
    {
        //! ----------------------------------------------
        //! remove the wireframe or transparent body view
        //! ----------------------------------------------
        //occContext->RemoveAll(false);

        //! ---------------------
        //! remove the mesh view
        //! ---------------------
        //occMeshContext->RemoveAll(false);
    }
        break;

    case resultPresentation::combinedView_meshVisible:
    {
        //! ----------------------------------------------
        //! remove the wireframe or transparent body view
        //! ----------------------------------------------
        //occContext->RemoveAll(false);

        const QMap<GeometryTag,occHandle(MeshVS_DeformedDataSource)> &mapOfMeshDS = aPostObject->getMeshDataSources();
        for(QMap<GeometryTag,occHandle(MeshVS_DeformedDataSource)>::const_iterator it = mapOfMeshDS.cbegin(); it!=mapOfMeshDS.cend(); it++)
        {
            cout<<"@ ----------------------------------"<<endl;
            cout<<"@ - number of nodes: "<<it.value()->GetAllNodes().Extent()<<endl;
            cout<<"@ - number of elements: "<<it.value()->GetAllElements().Extent()<<endl;
            cout<<"@ ----------------------------------"<<endl;

            occHandle(MeshVS_Mesh) aMeshObject = new MeshVS_Mesh();
            aMeshObject->SetDataSource(it.value());
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
        AIS_ListOfInteractive objInside;

        //! ---------------------------------------------------------------------
        //! remove the transparent body view: "0" is the signature for AIS_Shape
        //! ---------------------------------------------------------------------
        occContext->ObjectsInside(objInside,AIS_KOI_Shape,0);
        for(AIS_ListIteratorOfListOfInteractive it(objInside); it.More(); it.Next()) occContext->Remove(it.Value(),false);

        //! ---------------------
        //! remove the mesh view
        //! ---------------------
        //occMeshContext->RemoveAll(false);

        const QVector<GeometryTag> &vecLocs = aPostObject->getLocations();
        for(QVector<GeometryTag>::const_iterator it = vecLocs.cbegin(); it!=vecLocs.cend(); it++)
        {
            int bodyIndex = it->parentShapeNr;
            const TopoDS_Shape &aShape = myDS2->bodyMap.value(bodyIndex);
            const occHandle(AIS_Shape) &anAISShape = new AIS_Shape(aShape);
            occContext->SetColor(anAISShape,Quantity_NOC_BLACK,false);
            occContext->SetTransparency(anAISShape,0.0,false);
            occContext->Display(anAISShape,AIS_WireFrame,-1,false,false);
        }
    }
        break;

    case resultPresentation::combinedView_undeformedModel:
    {
        //! -------------------------------------------------------------------
        //! remove the wireframe body view: "0" is the signature for AIS_Shape
        //! -------------------------------------------------------------------
        //occContext->RemoveAll(false);

        //! ---------------------
        //! remove the mesh view
        //! ---------------------
        //occMeshContext->RemoveAll(false);

        const QVector<GeometryTag> &vecLocs = aPostObject->getLocations();
        for(QVector<GeometryTag>::const_iterator it = vecLocs.cbegin(); it!=vecLocs.cend(); it++)
        {
            int bodyIndex = it->parentShapeNr;
            const TopoDS_Shape &aShape = myDS2->bodyMap.value(bodyIndex);
            const occHandle(AIS_Shape) &anAISShape = new AIS_Shape(aShape);
            occContext->SetColor(anAISShape,Quantity_NOC_AQUAMARINE1,false);
            occContext->SetTransparency(anAISShape,0.9,true);
            occContext->Display(anAISShape,AIS_Shaded,-1,false,false);
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

//! --------------------------------
//! function: setResultPresentation
//! details:
//! --------------------------------
void occPostWidget::setResultPresentation(const resultPresentation &aResPresentation)
{
    myResultPresentation = aResPresentation;

    switch(myResultPresentation.theCombinedView)
    {
    case resultPresentation::combinedView_resultOnly: cout<<"occPostWidget::setResultPresentation()->____results only____"<<endl; break;
    case resultPresentation::combinedView_meshVisible: cout<<"occPostWidget::setResultPresentation()->____results with mesh____"<<endl; break;
    case resultPresentation::combinedView_undeformedModel: cout<<"occPostWidget::setResultPresentation()->____results and undeformed model____"<<endl; break;
    case resultPresentation::combinedView_undeformedWireFrame: cout<<"occPostWidget::setResultPresentation()->____results and undeformed wireframe____"<<endl; break;
    }

    cout<<"occPostWidget::setResultPresentation()->____scale: "<<aResPresentation.theScale<<"____"<<endl;

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
        this->emptyTheSelection();

        switch(myCurSelectionMode)
        {
        case CurSelection_Solid: occContext->DeactivateStandardMode(TopAbs_SOLID); break;
        case CurSelection_Face: occContext->DeactivateStandardMode(TopAbs_FACE); break;
        case CurSelection_Edge: occContext->DeactivateStandardMode(TopAbs_EDGE); break;
        case CurSelection_Vertex: occContext->DeactivateStandardMode(TopAbs_VERTEX); break;
        default: break;
        }

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
    //! ------------------------------------------------------
    //! replace all the post object meshes with volume meshes
    //! (when they are defined)
    //! ------------------------------------------------------
    occPreGLWidget::refreshMeshView(onlyExterior);
    isMeshViewVolume = (onlyExterior == true? false:true);
}
