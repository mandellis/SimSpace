//! ----------------
//! custom includes
//! ----------------
#include "cliptooldelegate.h"
#include "clipTool.h"

//! ---
//! Qt
//! ---
#include<QComboBox>
#include<QLineEdit>
#include<QSlider>

//! ----
//! C++
//! ----
#include<iostream>
using namespace std;

const QString sliderStyleSheet = QString("QSlider::groove:horizontal { "
                                 "border: 1px solid #999999; "
                                 "background: white; "
                                 "} "
                                 "QSlider::handle:horizontal { "
                                 "background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #b4b4b4, stop:1 #8f8f8f); "
                                 "border: 1px solid #5c5c5c; "
                                 "width: 18px; "
                                 "margin: 0px 0px; "
                                 "} ");

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
clipToolDelegate::clipToolDelegate(QWidget *parent):QStyledItemDelegate(parent)
{
    cout<<"clipToolDelegate::clipToolDelegate()->____CONSTRUCTOR CALLED____"<<endl;
}

//! -----------------------
//! function: createEditor
//! details:
//! -----------------------
QWidget* clipToolDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    Q_UNUSED(option)

    if(!index.isValid()) return 0;
    clipTool *theClipTool = static_cast<clipTool*>(this->parent());

    //! -----------------
    //! resize the cells
    //! -----------------
    //theClipTool->resizeRowsToContents();
    //theClipTool->resizeColumnsToContents();

    switch(index.column())
    {
    case CLIPPLANE_NAME_COLUMN:
    {
        QLineEdit *le = new QLineEdit(parent);
        return le;
    }
        break;

    case CLIPPLANE_ID_COLUMN:
    {
        return 0;
    }
        break;

    case CLIPPLANE_STATUS_COLUMN:
    {
        QVariant data;
        QComboBox *comboBox = new QComboBox(parent);
        data.setValue(false); comboBox->addItem("Off",data); //! combo-box index "0"
        data.setValue(true); comboBox->addItem("On",data);   //! combo-box index "1"
        return comboBox;
    }
        break;

    case CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN:
    {
        cout<<"clipToolDelegate::createEditor()->____create editor for column: "<<index.column()<<"____"<<endl;
        //! -----------------------------------------------
        //! return the editor only if the status is active
        //! -----------------------------------------------
        if(theClipTool->isCurrentRowActive())
        {
            //! -----------------------------
            //! clip plane coordinate system
            //! -----------------------------
            QComboBox *comboBox = new QComboBox(parent);

            QVariant data;
            QExtendedStandardItem *itemCSRoot = theClipTool->getCoordinateSystemRoot();
            for(int row=0; row<itemCSRoot->rowCount(); row++)
            {
                QExtendedStandardItem *curCSItem = static_cast<QExtendedStandardItem*>(itemCSRoot->child(row,0));
                void *p = (void*)curCSItem;
                data.setValue(p);
                const QIcon &anIcon = curCSItem->data(Qt::DecorationRole).value<QIcon>();
                comboBox->addItem(anIcon,curCSItem->data(Qt::DisplayRole).toString(),data);
            }
            return comboBox;
        }
        else return 0;
    }
        break;

    case CLIPPLANE_BASE_PLANE_DATA_COLUMN:
    {
        //! ----------------
        //! clip plane data
        //! ----------------
        return 0;
    }
        break;

    case CLIPPLANE_TRANSLATION_COLUMN:
    {
        //! -----------------------------------------------
        //! return the editor only if the status is active
        //! -----------------------------------------------
        if(theClipTool->isCurrentRowActive())
        {
            //! ------------------
            //! plane translation
            //! ------------------
            QSlider *editor = new QSlider(parent);
            editor->setStyleSheet(sliderStyleSheet);
            //editor->setFixedWidth(75);
            editor->setOrientation(Qt::Horizontal);
            editor->setMinimum(-100);
            editor->setMaximum(100);
            return editor;
        }
        else return 0;
    }
        break;

    default:
        return 0;
        break;
    }
}

//! ------------------------
//! function: setEditorData
//! details:
//! ------------------------
void clipToolDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    int col = index.column();
    switch(col)
    {
    case CLIPPLANE_NAME_COLUMN:
    {
        //! ---------------------------------------------
        //! initialize th editor for the clip plane name
        //! ---------------------------------------------
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        le->setText(index.data(Qt::UserRole).toString());
    }
        break;

    case CLIPPLANE_ID_COLUMN:
    {
        //! --------------------------
        //! clip plane ID: do nothing
        //! --------------------------
        ;
    }
        break;

    case CLIPPLANE_STATUS_COLUMN:
    {
        //! --------------------------------------------
        //! initialize the editor for clip plane status
        //! --------------------------------------------
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(index.data(Qt::UserRole).toInt());
        connect(comboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(handleCSStatusChanged()));
    }
        break;

    case CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN:
    {
        //! ------------------------------------------------
        //! initialize the editor for clip plane definition
        //! ------------------------------------------------
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QVariant data = index.data(Qt::UserRole);
        int comboBoxIndex = comboBox->findData(data,Qt::UserRole);
        if(comboBoxIndex!=-1) comboBox->setCurrentIndex(comboBoxIndex);
        else comboBox->setCurrentIndex(0);
        connect(comboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(handleCSChanged()));
    }
        break;

    case CLIPPLANE_BASE_PLANE_DATA_COLUMN:
    {
        //! ----------------------------
        //! clip plane data: do nothing
        //! ----------------------------
    }
        break;

    case CLIPPLANE_TRANSLATION_COLUMN:
    {
        //! ------------------
        //! plane translation
        //! ------------------
        QSlider *slider = static_cast<QSlider*>(editor);
        slider->setValue(index.data(Qt::UserRole).toInt());
        //connect(slider,SIGNAL(sliderMoved(int)),this,SLOT(handleCSTranslation(int)));
        connect(slider,SIGNAL(sliderReleased()),this,SLOT(handleCSTranslationHandleReleased()));
        connect(slider,SIGNAL(valueChanged(int)),this,SLOT(handleCSTranslation(int)));
    }
        break;
    }
}

//! -----------------------
//! function: setModelData
//! details:
//! -----------------------
void clipToolDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    QVariant data;
    switch(index.column())
    {
    case CLIPPLANE_NAME_COLUMN:
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        data.setValue(le->text());
        model->setData(index,data,Qt::UserRole);
    }
        break;

    case CLIPPLANE_ID_COLUMN:
    {
        ;
    }
        break;

    case CLIPPLANE_STATUS_COLUMN:
    {
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        data = comboBox->currentData(Qt::UserRole);
        model->setData(index,data,Qt::UserRole);
        if(data.toBool()) data.setValue(QString("On"));
        else data.setValue(QString("Off"));
        model->setData(index,data,Qt::DisplayRole);
    }
        break;

    case CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN:
    {
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QVariant data = comboBox->currentData(Qt::UserRole);
        model->setData(index,data,Qt::UserRole);
        QExtendedStandardItem *itemInVoidPointer = static_cast<QExtendedStandardItem*>(data.value<void*>());
        SimulationNodeClass *nodeInItem = itemInVoidPointer->data(Qt::UserRole).value<SimulationNodeClass*>();
        data.setValue(nodeInItem->getName());
        model->setData(index,data,Qt::DisplayRole);
    }
        break;

    case CLIPPLANE_BASE_PLANE_DATA_COLUMN:
    {
        ;
    }
        break;

    case CLIPPLANE_TRANSLATION_COLUMN:
    {
        QSlider *slider = static_cast<QSlider*>(editor);
        int val = slider->value();
        QVariant data;
        data.setValue(val);
        model->setData(index,data,Qt::UserRole);
        data.setValue(QString("%1").arg(val));
        model->setData(index,data,Qt::DisplayRole);
    }
        break;
    }
}

//! --------------------------
//! function: handleCSChanged
//! details:
//! --------------------------
void clipToolDelegate::handleCSChanged()
{
    cout<<"clipToolDelegate::handleCSChanged()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit currentCSChanged();
}

//! --------------------------------
//! function: handleCSStatusChanged
//! details:
//! --------------------------------
void clipToolDelegate::handleCSStatusChanged()
{
    cout<<"clipToolDelegate::handleCSStatusChanged()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit currentCSStatusChanged();
}

//! ------------------------------
//! function: handleCSTranslation
//! details:
//! ------------------------------
void clipToolDelegate::handleCSTranslation(int curZ)
{
    emit currentCSTranslationApplied(curZ);
}

//! --------------------------------------------
//! function: handleCSTranslationHandleReleased
//! details:
//! --------------------------------------------
void clipToolDelegate::handleCSTranslationHandleReleased()
{
    QSlider *editor = qobject_cast<QSlider*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}
