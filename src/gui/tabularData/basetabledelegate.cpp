#include "basetabledelegate.h"

//! ---
//! Qt
//! ---
#include <QLineEdit>
#include <QComboBox>

//! ----------------
//! custom includes
//! ----------------
#include "basetablemodel.h"

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
baseTableDelegate::baseTableDelegate(QWidget *parent):QStyledItemDelegate(parent) {}

//! -----------------------
//! function: createEditor
//! -----------------------
QWidget* baseTableDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    Q_UNUSED(option)
    if(!index.isValid()) return nullptr;

    QVariant data = index.data();
    switch(data.type())
    {
    case QMetaType::Bool:
    {
        QComboBox *cb = new QComboBox(parent);
        cb->addItem("false",false);
        cb->addItem("true",true);
        return cb;
    }
        break;
    case QMetaType::Int:
    {
        QLineEdit *le = new QLineEdit(parent);
        QIntValidator *v = new QIntValidator(le);
        le->setValidator(v);
        return le;
    }
        break;
    case QMetaType::Float:
    {
        QLineEdit *le = new QLineEdit(parent);
        QDoubleValidator *v = new QDoubleValidator(le);
        le->setValidator(v);
        return le;
    }
        break;
    case QMetaType::Double:
    {
        QLineEdit *le = new QLineEdit(parent);
        QDoubleValidator *v = new QDoubleValidator(le);
        le->setValidator(v);
        return le;
    }
        break;
    case QMetaType::QString:
    {
        QLineEdit *le = new QLineEdit(parent);
        return le;
    }
        break;
    default:
        return nullptr;
        break;
    }
}

//! ------------------------
//! function: setEditorData
//! details:
//! ------------------------
void baseTableDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    if(!index.isValid()) return;
    QVariant data = index.data(Qt::EditRole);
    switch(data.type())
    {
    case QMetaType::Bool:
    {
        QComboBox *cb = (QComboBox*)editor;
        if(data.toBool()==true) cb->setCurrentIndex(1);
        else cb->setCurrentIndex(0);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
        break;
    case QMetaType::Int:
    {
        QLineEdit *le = (QLineEdit*)editor;
        le->setText(QString("%1").arg(data.toInt()));
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdit()));
    }
        break;
    case QMetaType::Float:
    {
        QLineEdit *le = (QLineEdit*)editor;
        le->setText(QString("%1").arg(data.toInt()));
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdit()));
    }
        break;
    case QMetaType::Double:
    {
        QLineEdit *le = (QLineEdit*)editor;
        le->setText(QString("%1").arg(data.toInt()));
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdit()));
    }
        break;
    case QMetaType::QString:
    {
        QLineEdit *le = (QLineEdit*)editor;
        le->setText(data.toString());
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdit()));
    }
        break;
    default:
        break;
    }
}

//! -----------------------
//! function: setModelData
//! details:
//! -----------------------
void baseTableDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    if(!index.isValid()) return;
    QVariant data = index.data(Qt::EditRole);
    switch(data.type())
    {
    case QMetaType::Bool:
    {
        QComboBox *cb = (QComboBox*)editor;
        if(cb->currentIndex()==0) data.setValue(false);
        else data.setValue(true);
    }
        break;
    case QMetaType::Int:
    {
        QLineEdit *le = (QLineEdit*)editor;
        data.setValue(le->text().toInt());
        //cout<<"____"<<data.toInt()<<"____"<<endl;
    }
        break;
    case QMetaType::Float:
    {
        QLineEdit *le = (QLineEdit*)editor;
        data.setValue(le->text().toFloat());
    }
        break;
    case QMetaType::Double:
    {
        QLineEdit *le = (QLineEdit*)editor;
        data.setValue(le->text().toDouble());
    }
        break;
    case QMetaType::QString:
    {
        QLineEdit *le = (QLineEdit*)editor;
        data.setValue(le->text());
    }
        break;
    default: break;
    }
    if(!data.isValid()) return;
    if(data.isNull()) return;

    //cout<<"____setting data at index("<<index.row()<<", "<<index.column()<<")____"<<endl;
    (BaseTableModel*)model->setData(index, data, Qt::EditRole);
    //(BaseTableModel*)model->setData(index, data, Qt::DisplayRole);
}

//! ---------------------------------
//! function: commitAndCloseComboBox
//! details:
//! ---------------------------------
void baseTableDelegate::commitAndCloseComboBox()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit dataChanged();
}

//! ---------------------------------
//! function: commitAndCloseLineEdit
//! details:
//! ---------------------------------
void baseTableDelegate::commitAndCloseLineEdit()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit dataChanged();
}
