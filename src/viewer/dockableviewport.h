#ifndef DOCKABLEVIEWPORT_H
#define DOCKABLEVIEWPORT_H

//! ---
//! Qt
//! ---
#include <QWidget>
#include <QDockWidget>
#include <QFocusEvent>

//! ----------------
//! custom includes
//! ----------------
#include <occPreGLwidget.h>
#include "src/main/mydefines.h"
#include "simulationdatabase.h"
#include "occhandle.h"

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>

class dockableViewPort: public QDockWidget
{
    Q_OBJECT

public:

    dockableViewPort(QWidget *parent=0);

public:

    //! wrappers
    void setAction3D_Rotation();
    void setAction3D_Pan();
    void setAction3D_WindowZooming();
    void FitAll();
    void setSelectionMode(CurSelectionMode selectionMode);
    void setSimulationdataBase(simulationDataBase *sDB) { myViewPort->setSimulationdataBase(sDB); }
    void createInteractiveShapes() { myViewPort->createInteractiveShapes(); }
    void displayCAD();

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

    //! wrappers
    void hideBody(const TColStd_ListOfInteger &bodyListNbs);
    void showBody(const TColStd_ListOfInteger &bodyListNbs);
    const occHandle(AIS_InteractiveContext)& getContext() const;
    void displayShapeCopy(const TopTools_ListOfShape &list1,
                          const TopTools_ListOfShape &list2,
                          Quantity_NameOfColor color1,
                          Quantity_NameOfColor color2,
                          QVariant options=QVariant());
    void displayShapeCopy1(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color);
    void displayShapeCopy2(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color);

signals:

    void requestChangeDelegateContext(const occHandle(AIS_InteractiveContext) &aCTX);
};

#endif // DOCKABLEVIEWPORT_H
