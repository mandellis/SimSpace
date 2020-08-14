#ifndef MESHTOOLBAR_H
#define MESHTOOLBAR_H

//! ---
//! Qt
//! ---
#include <QToolBar>
#include <QString>
#include <QWidget>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

class QPushButtonExtended;

class MeshToolBar: public QToolBar
{
    Q_OBJECT

public:

    MeshToolBar(const QString &title, QWidget *parent = 0);

private:

    QPushButtonExtended* meshViewSelector;
    QAction *actionClearMesh;
    QAction *actionGenerateSurfaceMesh;
    QAction *actionGenerateVolumeMesh;

public slots:

    void enableMeshViewButton(bool isEnabled);
    void enableClearMesh(bool isEnabled) {actionClearMesh->setEnabled(isEnabled); }
    void enableSurfaceMeshButton(bool isEnabled) { actionGenerateSurfaceMesh->setEnabled(isEnabled); }
    void enableVolumeMeshButton(bool isEnabled) { actionGenerateVolumeMesh->setEnabled(isEnabled); }

private slots:

    void emitShowInteriorMeshRequested();
    void emitShowOnlySurfaceMeshRequested();
    void emitRequestClearMesh() { emit requestClearMesh(); }
    void emitRequestGenerateSurfaceMesh() { emit requestGenerateSurfaceMesh(); }
    void emitRequestGenerateVolumeMesh() { emit requestGenerateVolumeMesh(); }

signals:

    void showExteriorMeshRequest(bool showExterior);
    void requestClearMesh();
    void requestGenerateSurfaceMesh();
    void requestGenerateVolumeMesh();
    void requestUpdateViewerStatus();
};

#endif // MESHTOOLBAR
