#ifndef ITEMSELECTOR_H
#define ITEMSELECTOR_H

#include <QPushButton>
#include <QLineEdit>
#include <QEvent>

//! ------------------------
//! class QLineEditExtended
//! ------------------------
class QLineEditExtended: public QLineEdit
{
    Q_OBJECT

public:

    QLineEditExtended(QWidget *parent):QLineEdit(parent)
    {
        this->installEventFilter(this);
    }

    //! destructor
    virtual ~QLineEditExtended()
    {
        //! remove the event filter
        this->removeEventFilter(this);
    }

    //! event filter
    bool eventFilter(QObject* object, QEvent* event)
    {
        //! the left mouse button has been pressed
        if(object == this && event->type() == QEvent::MouseButtonPress)
        {
            //! bring up your custom edit
            emit activateSelector();
            return false; // lets the event continue to the edit
        }
        return QObject::event(event);
    }

signals:

    //! activate selector
    void activateSelector();
};


//! -------------------
//! class itemSelector
//! -------------------
#include <QStandardItem>
#include <QTreeView>
#include <QList>

class itemSelector: public QLineEditExtended
{
    Q_OBJECT

public:

    itemSelector(QWidget *parent=0);
    itemSelector(QTreeView *aTreeView, QWidget *parent=0);

public:

    QList<void*> getObjects() { return objects; }
    void setObjects(const QList<void*> objectList);
    void setTreeView(QTreeView *aTreeView) { treeView = aTreeView; }
    QTreeView* getTreeView(){ return treeView; }

private:

    QPushButton *bap;
    QPushButton *bcn;
    QList<void*> objects;
    QTreeView *treeView;

private:

    void createContent();
    void createConnections();
    void updateText();

private slots:

    void showPushButtons();
    void hidePushButtons();
    void setAccepted();
    void setRejected();

    void setStrongFocus();

signals:

};

#endif // ITEMSELECTOR_H
