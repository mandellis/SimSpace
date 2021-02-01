#ifndef DIRECTIONSELECTOR_H
#define DIRECTIONSELECTOR_H

#define DIRECTION_SELECTOR_COLOR Quantity_NOC_GREEN

//! Qt
#include <QWidget>

//! custom includes
#include "shapeselectorbox.h"
#include "src/registeredMetatypes/listofshape.h"
#include "ext/occ_extended/handle_ais_arrowmarker_reg.h"

//! OCC
#include <AIS_InteractiveContext.hxx>

class QPushButton;
class QHBoxLayout;
class QGridLayout;

class DirectionSelector : public ShapeSelectorBox
{
    Q_OBJECT

public:

    explicit DirectionSelector(QWidget *parent=0);
    explicit DirectionSelector(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent=0);
    ~DirectionSelector();

private:

    QPushButton *bap;
    QPushButton *bcn;
    QHBoxLayout *h;
    QGridLayout *g;

    ListOfShape myShape;
    QVector<double> myDirection;
    AIS_ArrowMarker_handle_reg myMarker;
    occHandle(AIS_InteractiveContext) myCTX;

private:

    void showPushButtons();
    void hidePushButtons();

private slots:

    void createContent();
    void setAccepted();
    void setRejected();
    void showMarker();

public:

    void setDirection(QVector<double> vec);
    void setContext(const occHandle(AIS_InteractiveContext) &aCTX);
    QVector<double> getDirection();

signals:

    void editingFinished();
};

#endif // DIRECTIONSELECTOR_H
