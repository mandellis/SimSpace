#ifndef DOCKABLEVIEWPORT_H
#define DOCKABLEVIEWPORT_H

//! Qt
#include <QWidget>
#include <QDockWidget>
#include <QFocusEvent>

//! custom includes
#include <occPreGLwidget.h>
#include "mydefines.h"
#include "simulationdatabase.h"

//! OCC
#include <AIS_InteractiveContext.hxx>

class dockableViewPort: public QDockWidget
{
    Q_OBJECT

public:

    dockableViewPort(QWidget *parent=0);

public:

    //! interfacing functions
    void setAction3D_Rotation();
    void setAction3D_Pan();
    void setAction3D_WindowZooming();
    void FitAll();
    void setSelectionMode(CurSelectionMode selectionMode);

    //! set the simulation database - interface
    void setSimulationdataBase(simulationDataBase *sDB) { myViewPort->setSimulationdataBase(sDB); }

    //! create the interactive shapes - interface
    void createInteractiveShapes() { myViewPort->createInteractiveShapes(); }

    //! display CAD - interface
    void displayCAD() { myViewPort->displayCAD(); }

private:

    occPreGLWidget *myViewPort;
    QWidget *myWidget;

private:

    void createContent();

protected:

    virtual void focusInEvent(QFocusEvent *focusEvent) override;

public:

    //! get the view-port: access function
    occPreGLWidget* getViewPort() { return myViewPort; }

public slots:

    //! hide body - interface
    void hideBody(const TColStd_ListOfInteger &bodyListNbs) { myViewPort->hideBody(bodyListNbs); }

    //! show body - interface
    void showBody(const TColStd_ListOfInteger &bodyListNbs) { myViewPort->showBody(bodyListNbs); }

    //! get the interactive context - interface
    const occHandle(AIS_InteractiveContext)& getContext() const { return myViewPort->getContext(); }

    //! displayShapeCopy - interface
    void displayShapeCopy(const TopTools_ListOfShape &list1,
                          const TopTools_ListOfShape &list2,
                          Quantity_NameOfColor color1,
                          Quantity_NameOfColor color2,
                          QVariant options=QVariant())
    {
        myViewPort->displayShapeCopy(list1,list2,color1,color2,options);
    }

    //! displayShapeCopy1 - interface
    void displayShapeCopy1(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color)
    {
        myViewPort->displayShapeCopy1(listShapes,color);
    }

    //! displayShapeCopy2 - interface
    void displayShapeCopy2(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color)
    {
        myViewPort->displayShapeCopy2(listShapes,color);
    }

signals:

    void requestChangeDelegateContext(const occHandle(AIS_InteractiveContext) &aCTX);
};

#endif // DOCKABLEVIEWPORT_H
