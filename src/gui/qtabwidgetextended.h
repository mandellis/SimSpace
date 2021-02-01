#ifndef QTABWIDGETEXTENDED_H
#define QTABWIDGETEXTENDED_H

#include <QTabWidget>
#include <iostream>

class QTabWidgetExtended : public QTabWidget
{
    Q_OBJECT

public:

    QTabWidgetExtended(QWidget* parent=0):QTabWidget(parent)
    {
        std::cout<<"QTabWidgetExtended()->____CONSTRUCTOR CALLED____"<<std::endl;
    }
    QTabWidgetExtended::~QTabWidgetExtended()
    {
        std::cout<<"QTabWidgetExtended::~QTabWidgetExtended()->____DESTRUCTOR CALLED____"<<std::endl;
    }

protected:

    virtual void resizeEvent(QEvent *e)
    {
        Q_UNUSED(e)
        emit resized();
    }

public slots:

    void setCurrentTab(const QString &tabName)
    {
        foreach (QWidget* child, findChildren<QWidget*>())
        {
            if (child->objectName() == tabName)
            {
                setCurrentIndex(indexOf(child));
            }
        }
    }

signals:

    void resized();
};

#endif // QTABWIDGETEXTENDED_H
