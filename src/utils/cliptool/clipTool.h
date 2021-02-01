#ifndef CLIPTOOL_H
#define CLIPTOOL_H

#define CLIPPLANE_NAME_COLUMN 0
#define CLIPPLANE_ID_COLUMN 1
#define CLIPPLANE_STATUS_COLUMN 2
#define CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN 3
#define CLIPPLANE_BASE_PLANE_DATA_COLUMN 4
#define CLIPPLANE_TRANSLATION_COLUMN 5
#define CLIPPLANE_SHIFTED_PLANE_COEFFICIENTS_COLUMN 6

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>
#include <MeshVS_DataSource.hxx>
#include <MeshVS_Mesh.hxx>
#include <Quantity_NameOfColor.hxx>
#include <TColStd_HPackedMapOfInteger.hxx>

//! ---
//! Qt
//! ---
#include <QWidget>
#include <QTableView>
#include <QStandardItemModel>
#include <qextendedstandarditem.h> //! this type is registered
#include <QMouseEvent>
#include <QHBoxLayout>

//! ----------------
//! custom includes
//! ----------------
#include "occpostwidget.h"
//#include "occPreGLwidget.h"
#include <meshdatabase.h>


//! ----
//! C++
//! ----
#include <map>
#include <iostream>
using namespace std;

class clipTool: public QTableView
{
    Q_OBJECT

public:

    explicit clipTool(QWidget *parent = 0);
    ~clipTool(){ cout<<"clipTool::~clipTool()->____DESTRUCTOR CALLED____"<<endl;}

protected:

    //! -------------------
    //! handle right click
    //! -------------------
    virtual void mousePressEvent(QMouseEvent *event)
    {
        if (event->button() == Qt::RightButton)
        {
            myCurPoint = event->pos();
        }
        if (event->button() == Qt::LeftButton)
        {
            QModelIndex index = indexAt(event->pos());
            if (index.column() == 1)
            {
                edit(index);
            }
        }

        // pass on other buttons to base class
        QTableView::mousePressEvent(event);
    }

private slots:

    void showContextMenu(QPoint aPoint);
    void updateClipPlaneOfRow();
    void handleClipPlaneOfRowStatus();

    void updateCSTranslation(int sliderPosition);
    void updateCSDefinition();

public slots:

    void updateCSDataByExternalCSChange(QStandardItem *theCurrentModifiedCS);
    //void updateClippedMeshView(bool onlyExterior);

private:

    void addItemToTable();
    void removeItemFromTable();

private:

    QPoint myCurPoint;
    QStandardItemModel *internalModel;
    QExtendedStandardItem *myCoordinateSystemRoot;
    //occPreGLWidget *myOCCViewer;
    occPostWidget *myOCCViewer;
    int myMaxClipPlanes;
    int myCurNumberOfClipPlanes;

    meshDataBase *myMDB;
    std::map<int,occHandle(TColStd_HPackedMapOfInteger)> myHiddenElements;

public:

    //! -------------------------------
    //! set the coordinate system root
    //! -------------------------------
    void setCoordinateSystemRoot(QStandardItem *itemCSRoot)
    {
        if(itemCSRoot!=NULL) myCoordinateSystemRoot = static_cast<QExtendedStandardItem*>(itemCSRoot);
    }

    //! ---------------
    //! set the viewer
    //! ---------------
    //void setViewer(occPreGLWidget *anOCCViewer)
    void setViewer(occPostWidget *anOCCViewer)
    {
        myOCCViewer = anOCCViewer;

        //! get the limit of clipping planes for the current view
        myMaxClipPlanes = myOCCViewer->getClipPlanesNbLimit();
    }

    //! --------------------------------
    //! get the coordinate systems root
    //! --------------------------------
    QExtendedStandardItem* getCoordinateSystemRoot() { return myCoordinateSystemRoot; }

    //! -------------------------------
    //! get the status of a clip plane
    //! -------------------------------
    bool isCurrentRowActive();

private:

    //! return the plane equation coefficients
    QVector<double> getPlaneCoefficients(QStandardItem *aCSItem);

    //! ------------------------------------------------
    //! working mode
    //! 0 - Mesh; 1 - contacts; 2 - model; 3 - solution
    //! ------------------------------------------------
    int myWorkingMode;

public slots:

    //! set the geometry/mesh database
    void setMeshDataBase(meshDataBase *aMeshDataBase);

    //! set the working mode
    void setWorkingMode(int workingMode);

    //! compute hidden elements
    //void computeHiddenElements();

signals:

    void clipPlaneAdded();
    void clipPlaneRemoved();
    void clipPlaneChanged();
    void clipPlaneEnabled(int,bool);

private:

    //! retrieve the map of clip planes (ID, {coefficients}
    int retrieveActiveClipPlanes(std::map<int, std::vector<double> > &mapOfClipPlanes);

    //! set current clip plane
    void setCurrentClipPlane(int curClipPlaneID);

    //! translate plane
    void translatePlane(double &a, double &b, double &c, double &d, const double &t);
};

#endif // CLIPTOOL_H
