//! ----------------
//! custom includes
//! ----------------
#include "meshtoolbar.h"
#include "qpushbuttonextended.h"
#include "global.h"

//! ---
//! Qt
//! ---
#include <QMenu>
#include <QAction>
#include <QMessageBox>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
MeshToolBar::MeshToolBar(const QString &title, QWidget *parent):QToolBar(title,parent)
{
    //! ----------------
    //! set object name
    //! ----------------
    this->setObjectName("meshToolBar");

    meshViewSelector = new QPushButtonExtended(this);
    meshViewSelector->setIcon(QIcon(":/icons/icon_surface mesh.png"));
    meshViewSelector->setToolTip("Enable disable mesh interior view");

    //! ----------------
    //! the button menu
    //! ----------------
    QMenu *menu = new QMenu(this);

    //! --------------------------
    //! action show interior mesh
    //! --------------------------
    QAction *showInteriorMesh = new QAction("Show volume mesh",this);
    showInteriorMesh->setIcon(QIcon(":/icons/icon_volume mesh.png"));
    connect(showInteriorMesh,SIGNAL(triggered(bool)),this,SLOT(emitShowInteriorMeshRequested()));

    //! ------------------------------
    //! action show only surface mesh
    //! ------------------------------
    QAction *showOnlySurfaceMesh = new QAction("Show surface mesh",this);
    showOnlySurfaceMesh->setIcon(QIcon(":/icons/icon_surface mesh.png"));
    connect(showOnlySurfaceMesh,SIGNAL(triggered(bool)),this,SLOT(emitShowOnlySurfaceMeshRequested()));

    menu->addAction(showOnlySurfaceMesh);
    menu->addAction(showInteriorMesh);

    meshViewSelector->setMenu(menu);
    this->addWidget(meshViewSelector);

    //! --------------------
    //! clear meshes button
    //! --------------------
    actionClearMesh = new QAction(this);
    actionClearMesh->setIcon(QIcon(":/icons/icon_clear data.png"));
    actionClearMesh->setToolTip("Clear the mesh on the selected body");
    this->addAction(actionClearMesh);
    connect(actionClearMesh,SIGNAL(triggered(bool)),this,SLOT(emitRequestClearMesh()));

    //! ----------------------
    //! generate surface mesh
    //! ----------------------
    actionGenerateSurfaceMesh = new QAction(this);
    actionGenerateSurfaceMesh->setIcon(QIcon(":/icons/icon_surface mesh.png"));
    actionGenerateSurfaceMesh->setToolTip("Generate the surface mesh on the selected body");
    this->addAction(actionGenerateSurfaceMesh);
    connect(actionGenerateSurfaceMesh,SIGNAL(triggered(bool)),this,SLOT(emitRequestGenerateSurfaceMesh()));

    //! ---------------------
    //! generate volume mesh
    //! ---------------------
    actionGenerateVolumeMesh = new QAction(this);
    actionGenerateVolumeMesh->setIcon(QIcon(":/icons/icon_volume mesh.png"));
    actionGenerateVolumeMesh->setToolTip("Generate the volume mesh on the selected body");
    this->addAction(actionGenerateVolumeMesh);
    connect(actionGenerateVolumeMesh,SIGNAL(triggered(bool)),this,SLOT(emitRequestGenerateVolumeMesh()));
}

//! ----------------------------------------
//! function: emitShowInteriorMeshRequested
//! details:
//! ----------------------------------------
void MeshToolBar::emitShowInteriorMeshRequested()
{
    meshViewSelector->setIcon(QIcon(":/icons/icon_volume mesh.png"));
    Global::status().myResultPresentation.useExteriorMeshForVolumeResults = false;
    emit showExteriorMeshRequest(false);
    emit requestUpdateViewerStatus();
}

//! -------------------------------------------
//! function: emitShowOnlySurfaceMeshRequested
//! details:
//! -------------------------------------------
void MeshToolBar::emitShowOnlySurfaceMeshRequested()
{
    meshViewSelector->setIcon(QIcon(":/icons/icon_surface mesh.png"));
    Global::status().myResultPresentation.useExteriorMeshForVolumeResults = true;
    emit showExteriorMeshRequest(true);
    emit requestUpdateViewerStatus();
}

//! -------------------------------
//! function: enableMeshViewButton
//! details:
//! -------------------------------
void MeshToolBar::enableMeshViewButton(bool isEnabled)
{
    meshViewSelector->setEnabled(isEnabled);
}
