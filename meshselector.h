#ifndef MESHSELECTOR_H
#define MESHSELECTOR_H

#include <AIS_InteractiveContext.hxx>
#include "shapeselectorbox.h"
#include "occhandle.h"
#include <QWidget>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QPushButton>

#include <meshelementbycoords.h>

class MeshSelector : public ShapeSelectorBox
{
    Q_OBJECT

public:

    explicit MeshSelector(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent);
    explicit MeshSelector(QWidget *parent = 0);
    virtual ~MeshSelector();

private:

    QPushButton *bap;
    QPushButton *bcn;
    QHBoxLayout *h;
    QGridLayout *g;

    occHandle(AIS_InteractiveContext) myCTX;
    std::vector<meshElementByCoords> myElements;
    std::vector<int> myIDs;

private:

    void showPushButtons();
    void hidePushButtons();
    void clearContext(bool updateViewer=true);

private slots:

    void createContent();
    void setAccepted();
    void setRejected();

public:

    void setElements(const std::vector<meshElementByCoords> &vecElements);
    void setContext(const occHandle(AIS_InteractiveContext) &aCTX);

    std::vector<meshElementByCoords> getElements() const;

signals:

    void editingSelectionFinished();
};

#endif // MESHSELECTOR_H
