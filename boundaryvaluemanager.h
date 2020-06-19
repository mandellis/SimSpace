#ifndef BOUNDARYVALUEMANAGER_H
#define BOUNDARYVALUEMANAGER_H

//! custom includes
#include "shapeselectorbox.h"
#include "myenumvariables.h"
#include "property.h"

//! Qt
#include <QWidget>

class QMenu;
class QPushButton;
class ShapeSelectorBox;
class QAction;

class BoundaryValueManager: public QWidget
{
    Q_OBJECT

public:

    //! constructor I
    explicit BoundaryValueManager(QWidget *parent=0);

private:

    QMenu *myContextMenu;
    QPushButton *but;
    ShapeSelectorBox *lineEdit;

    QAction* actionConstant;
    QAction* actionTabular;
    QAction* actionFunction;
    QAction* actionFree;

    Property::typeOfValue myTypeOfValue;

private slots:

    void showContextMenu();
    void buildContextMenu();

public slots:

    void setValue(Property::typeOfValue theTypeOfValue);
    Property::typeOfValue getValue() const
    {
        return myTypeOfValue;
    }

signals:

    void typeOfValue_constant_selected();
    void typeOfValue_tabular_selected();
    void typeOfValue_function_selected();
    void typeOfValue_free_selected();
    void editingFinished();
};

#endif // BOUNDARYVALUEMANAGER_H
