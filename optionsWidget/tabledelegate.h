#ifndef TABLEDELEGATE_H
#define TABLEDELEGATE_H

#include <QStyledItemDelegate>

class tableDelegate: public QStyledItemDelegate
{
public:

    //! constructor
    explicit tableDelegate(QWidget *parent=0);

    //! Create editor
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

    //! Set the editor
    void setEditorData(QWidget *editor, const QModelIndex &index) const override;

    //! Set model data
   void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;
};

#endif // TABLEDELEGATE_H
