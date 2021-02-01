#ifndef SHAPESELECTORBOX_H
#define SHAPESELECTORBOX_H

#include <QLineEdit>

class ShapeSelectorBox : public QLineEdit
{
    Q_OBJECT

public:

    //! constructor
    ShapeSelectorBox(QWidget *parent = 0);

    //! destructor
    virtual ~ShapeSelectorBox()
    {
        //! remove the event filter
        this->removeEventFilter(this);
    }

    //! event filter
    bool eventFilter(QObject* object, QEvent* event);

signals:

    //! activate selector
    void activateSelector();
};

#endif // SHAPESELECTORBOX_H
