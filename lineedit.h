#ifndef LINEEDIT_H
#define LINEEDIT_H

//! Qt
#include <QLineEdit>

//! custom includes
#include <myenumvariables.h>
#include "property.h"

class QToolButton;

class LineEdit : public QLineEdit
{
    Q_OBJECT

public:

    LineEdit(bool hasFreeOption, QWidget *parent = 0);
    QString getText() const;

protected:

    void resizeEvent(QResizeEvent *);

private slots:

    void createContent(bool hasFreeOption);
    void showContextMenu();
    void emitEditingFinished();
    void emitEditingFinishedAfterTextChange();
    void changeMyText();
    void unselect() { this->setSelection(0,0); }

public slots:

    void setData(Property::loadDefinition theLoadDefinition);
    Property::loadDefinition getLoadDefinition() const;

private:

    QToolButton *Button;
    QMenu *myContextMenu;
    QAction* actionConstant;
    QAction* actionTabular;
    QAction* actionFunction;
    QAction* actionFree;
    QString myTextData;

    bool hasFreeOption;
    Property::loadDefinition myLoadDefinition;

signals:

    void editingFinished(const QString& theTextValue);
};

#endif // LIENEDIT_H
