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
    myResultPresentation.isMeshVisible = false;
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
void occPostWidget::displayResult(const postObject &aPostObject)
{
    cout<<"occPostWidget::displayResult()->____function called____"<<endl;

    //! -----------------------------------
    //! display the list of colored meshes
    //! -----------------------------------
    QMap<GeometryTag,occHandle(MeshVS_Mesh)>::const_iterator anIt;
    QMap<GeometryTag,occHandle(MeshVS_Mesh)> coloredMeshes = aPostObject.getColoredMeshes();
    for(anIt = coloredMeshes.cbegin(); anIt != coloredMeshes.cend(); ++anIt)
    {
        const occHandle(MeshVS_Mesh) &aColoredMesh = anIt.value();
        occPostContext->Display(aColoredMesh,false);
    }

    //! ----------------------
    //! display the color box
    //! ----------------------
    const occHandle(AIS_ColorScaleExtended) &colorBox = aPostObject.getColorBox();
    occPostContext->Display(colorBox,false);

    //! ------------------
    //! update the viewer
    //! ------------------
    occPostContext->UpdateCurrentViewer();
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

//! ---------------------------
//! function: hideSingleResult
//! details:
//! ----------------------------
void occPostWidget::hideSingleResult(postObject curPostObject)
{
    //cout<<"occPostWidget::hideSingleResult()->____function called____"<<endl;
    occPostContext->RemoveAll(true);
}

//! -----------------------
//! function: updateResult
//! details:
//! -----------------------
void occPostWidget::updateResult(postObject &aPostObject)
{
    aPostObject.update(static_cast<meshDataBase*>(myDS2));
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
        //this->setWireframeView();
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
