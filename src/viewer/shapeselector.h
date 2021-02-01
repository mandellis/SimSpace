#ifndef SHAPESELECTOR_H
#define SHAPESELECTOR_H

//! ---
//! Qt
//! ---
#include <QWidget>

//! ----------------
//! custom includes
//! ----------------
#include "shapeselectorbox.h"
#include "src/registeredMetatypes/listofshape.h"
#include "geometrytag.h"
#include "occhandle.h"

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>
#include <TopAbs_ShapeEnum.hxx>

class QPushButton;
class QHBoxLayout;
class QGridLayout;
class geometryDataBase;

class ShapeSelector : public ShapeSelectorBox
{
    Q_OBJECT

public:

    explicit ShapeSelector(QWidget *parent=0);
    explicit ShapeSelector(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent=0);
    virtual ~ShapeSelector();

private:

    QPushButton *bap;
    QPushButton *bcn;
    QHBoxLayout *h;
    QGridLayout *g;

    ListOfShape myShapes;
    std::vector<GeometryTag> myVecLoc;

    occHandle(AIS_InteractiveContext) myCTX;

    geometryDataBase *gdb;

private:

    void showPushButtons();
    void hidePushButtons();
    void clearContext(bool updateViewer=true);

private slots:

    void createContent();
    void setAccepted();
    void setRejected();

public:

    void setShape(const std::vector<GeometryTag> &vecLoc);
    void setContext(const occHandle(AIS_InteractiveContext) &aCTX);

    ListOfShape getShape() const;
    std::vector<GeometryTag> getVecLoc() const;

signals:

    //void editingFinished();
    void editingSelectionFinished();
};

#endif // SHAPESELECTOR_H
