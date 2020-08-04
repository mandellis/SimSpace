#ifndef OCCPOSTWIDGET_H
#define OCCPOSTWIDGET_H

//! ----------------
//! custom includes
//! ----------------
#include <occPreGLwidget.h>
#include "ais_colorscaleextended.h"
#include "listofmesh.h"
#include "mydefines.h"
#include "postobject.h"
#include "frdreader.h"
#include <meshdatabase.h>
#include "resultpresentation.h"

//! ---
//! Qt
//! ---
#include <QWidget>

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>
#include <AIS_ListOfInteractive.hxx>
#include <MeshVS_DataSource.hxx>

class occPostWidget: public occPreGLWidget
{
    Q_OBJECT

public:

    //! constructor
    occPostWidget(QWidget *parent=0);

    //! constructor
    occPostWidget(meshDataBase *mDB, QWidget *parent=0);

    //! destructor
    virtual ~occPostWidget()
    {
        cout<<"occPostWidget::~occPostWidget()->____DESTRUCTOR CALLED____"<<endl;
    }

    //! the mesh data base
    meshDataBase *myMeshDataBase;

    //! set database
    void setDataDase(meshDataBase *mDB) { myMeshDataBase = mDB; }

private:

    //! Creates an additional context for the mesh view
    //occHandle(AIS_InteractiveContext) occPostContext;

    //! the color box
    occHandle(AIS_ColorScaleExtended) myColorBox;

public slots:

    void displayColorBox(bool isVisible);

    void createColorBox(double min, double max, int Nintervals);

    //! it read the private member myResultPresentation before displaying
    void displayResult(const postObject &aPostObject);

    void hideAllResults();

    void hideSingleResult(postObject);

    void updateResult(postObject &aPostObject);

    //! Activate the working mode solution
    void setWorkingMode_Solution();

    //! refresh mesh view
    virtual void refreshMeshView(bool onlyExterior);

    //! set the status variable
    void setResultPresentation(const resultPresentation &aResPresentation) { myResultPresentation = aResPresentation; }

private:

    //! status variable
    resultPresentation myResultPresentation;
};

#endif // OCCPOSTWIDGET_H
