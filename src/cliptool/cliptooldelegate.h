#ifndef CLIPTOOLDELEGATE_H
#define CLIPTOOLDELEGATE_H

#include <QStyledItemDelegate>

class clipToolDelegate: public QStyledItemDelegate
{
    Q_OBJECT

public:

    //! constructor
    explicit clipToolDelegate(QWidget *parent=0);

    //! Create editor
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

    //! Set the editor
    void setEditorData(QWidget *editor, const QModelIndex &index) const override;

    //! Set model data
   void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;

private slots:

   void handleCSChanged();
   void handleCSStatusChanged();
   void handleCSTranslation(int sliderValue);
   void handleCSTranslationHandleReleased();

signals:

   void currentCSChanged();
   void currentCSStatusChanged();
   void currentCSTranslationApplied(int);
};

#endif // CLIPTOOLDELEGATE_H
