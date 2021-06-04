#ifndef BASETABLEDELEGATE_H
#define BASETABLEDELEGATE_H

#include "basetabledelegate.h"

//! ---
//! Qt
//! ---
#include <QWidget>
#include <QStyledItemDelegate>

class baseTableDelegate: public QStyledItemDelegate
{
    Q_OBJECT

public:

    //! constructor
    baseTableDelegate(QWidget *parent=0);

    //! destructor
    ~baseTableDelegate(){}

    //! create editor
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

    //! set editor data
    void setEditorData(QWidget *editor, const QModelIndex &index) const override;

    //! set model data
    void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;

signals:

    //! data changed
    void dataChanged();

private slots:

    //! commint and close line edit editor
    void commitAndCloseLineEdit();

    //! commit and close combobox editor
    void commitAndCloseComboBox();

};

#endif // BASETABLEDELEGATE_H
