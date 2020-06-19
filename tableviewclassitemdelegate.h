#ifndef TABLEVIEWCLASSITEMDELEGATE_H
#define TABLEVIEWCLASSITEMDELEGATE_H

//! Qt
#include <QWidget>
#include <QStyledItemDelegate>

class tableViewClassItemDelegate: public QStyledItemDelegate
{
    Q_OBJECT

public:

    //! constructor
    explicit tableViewClassItemDelegate(QWidget *parent=0);

    //! Create Editor when constructing the delegate
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

    //! Set the Editor
    void setEditorData(QWidget *editor, const QModelIndex &index) const override;

    //! When data are modified,  this model reflect the change
   void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;

private slots:

   //! commit and close line edit
   void commitAndCloseLineEdit();

   //! commit and close comboBox
   void commitAndCloseComboBox();

   //! commit and close comboBox
   void commitAndCloseBoltStatusDefineBy();

signals:

   void boltStatusDefineByChanged();


};

#endif // TABLEVIEWCLASSITEMDELEGATE_H
