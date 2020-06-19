//! custom includes
#include "tabledelegate.h"
#include "colorselector.h"
#include "celldata.h"
#include "viewoptions.h"

//! Qt
#include <QColor>
#include <QLineEdit>
#include <QComboBox>

//! C++
#include <iostream>

using namespace std;

//! ------------------------------------------------------------------------------
//! function: constructor
//! details:
//! ------------------------------------------------------------------------------
tableDelegate::tableDelegate(QWidget *parent): QStyledItemDelegate(parent)
{
    ;
}

//! ------------------------------------------------------------------------------
//! function: createEditor
//! details:
//! ------------------------------------------------------------------------------
QWidget* tableDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    cout<<"tableDelegate::createEditor()->____function called___"<<endl;
    Q_UNUSED(option)

    //! allow modification only for second column
    if(!index.isValid() || index.column()!=1) return 0;

    QVariant data = index.data(Qt::UserRole);
    if(data.isValid())
    {
        cout<<data.typeName()<<endl;
        if(data.canConvert<viewOptions>())
        {
            QComboBox *cb = new QComboBox(parent);
            cb->addItem("Uniform",0);
            cb->addItem("Horizontal",1);
            cb->addItem("Vertical",2);
            return cb;
        }
        else if(data.canConvert<QColor>())
        {
            cout<<"tableDelegate::createEditor____creating editor for colors____"<<endl;
            ColorSelector *colSelector = new ColorSelector(parent);
            return colSelector;
        }
        else if(data.canConvert<QString>())
        {
            cout<<"tableDelegate::createEditor()->____create editors for strings____"<<endl;
            QLineEdit *le = new QLineEdit(parent);
            return le;
        }
        else return 0;
    }
    return 0;
}

//! ------------------------------------------------------------------------------
//! function: setEditorData
//! details:
//! ------------------------------------------------------------------------------
void tableDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    cout<<"tableDelegate::setEditorData()->____function called____"<<endl;
    QVariant data = index.model()->data(index, Qt::UserRole);
    cout<<"tableDelegate::setEditorData()->____initialize editor for: "<<data.typeName()<<"____"<<endl;

    if(data.canConvert<QColor>())
    {
        ColorSelector *cs = static_cast<ColorSelector*>(editor);
        cs->setColor(data.value<QColor>());
    }
    else if(data.canConvert<QString>())
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        le->setText(data.toString());
    }
    else if(data.canConvert<viewOptions>())
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(data.value<viewOptions>());
    }
}

//! ------------------------------------------------------------------------------
//! function: setModelData
//! details:
//! ------------------------------------------------------------------------------
void tableDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    if(!index.isValid())return;
    QVariant data = index.data(Qt::UserRole);
    if(data.canConvert<QColor>())
    {
        ColorSelector *cs = static_cast<ColorSelector*>(editor);
        data.setValue(cs->color());
    }
    else if(data.canConvert<QString>())
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        data.setValue(le->text());
    }
    else if(data.canConvert<viewOptions>())
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        data.setValue(static_cast<viewOptions>(cb->currentIndex()));
    }

    //! change the data of the model @index position
    model->setData(index, data, Qt::UserRole);
}
