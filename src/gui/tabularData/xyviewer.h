#ifndef XYVIEWER_H
#define XYVIEWER_H

#define _MATH_DEFINES_DEFINED

//! ---
//! Qt
//! ---
#include <QWidget>
#include <QAbstractItemView>
#include <QVector>

//! ----------------
//! custom includes
//! ----------------
#include <QCustomPlot/qcp/qcustomplot.h>
#include "basetablemodel.h"

//! --------
//! STL C++
//! --------
#include <map>

class QComboBox;
class QMenu;

class XYViewer: public QAbstractItemView
{
    Q_OBJECT

public:

    XYViewer(QWidget *parent = nullptr);
    virtual void setModel(BaseTableModel *model);

protected:

    QSize viewportSizeHint() const override;
    virtual void paintEvent(QPaintEvent *event) override;

    //! ---------------------------------
    //! re-implemented virtual functions
    //! ---------------------------------
    virtual int horizontalOffset() const override;
    virtual int verticalOffset() const override;
    virtual QModelIndex indexAt(const QPoint &point) const override;
    virtual bool isIndexHidden(const QModelIndex &index) const override;
    virtual QModelIndex moveCursor(CursorAction cursorAction, Qt::KeyboardModifiers modifiers) override;
    virtual void scrollTo(const QModelIndex &index, ScrollHint hint = EnsureVisible) override;
    virtual void setSelection(const QRect &rect, QItemSelectionModel::SelectionFlags flags) override;
    virtual QRect visualRect(const QModelIndex &index) const override;
    virtual QRegion visualRegionForSelection(const QItemSelection &selection) const override;

private:

    QCustomPlot *m_plot;
    //std::map<int,std::pair<int,int>> m_graphs;          //! key => graph number value=> h series index, v series index
    std::map<int,QVector<double>> m_series;             //! key => column value => series values

    //! table data model
    BaseTableModel *m_model;

    //! table of graph controls
    BaseTableModel *m_controls;

    //! controls
    QTableView *m_tableOfControls;

    //! menu
    QMenu *m_contextMenu;

    //! indexed map of graphs (from "0")
    std::map<QCPGraph*,int> m_mapOfGraphs;

    //! current selected plot
    QCPAbstractPlottable* m_selectedPlottable;

public:

    //! add plot
    bool addGraph(int xDataCol, int yDataCol);

private:

    //! build custom menu
    void buildCustomMenu();

    //! refresh all series
    void refreshAllSeries();

    //! remove additional axes
    void removeAdditionalAxes();


protected slots:

    //! data changed
    virtual void dataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight, const QVector<int> &roles = QVector<int>()) override;

private slots:

    //! show context menu
    void showContextMenu(const QPoint &);

    //! hide all plot
    void hideAll();

    //! hide selected plot
    void hideSelected();

    //! add plot
    void addPlot();

    //! selection changed by user
    void selectionChangedByUser();

    //! hide controls
    void hideControls();

    //! show controls
    void showControls();

    //! replot all
    void replotAll();
};


//! ----------------
//! class: m_ticker
//! ----------------
class m_ticker: public QCPAxisTicker
{
private:

    double m_min;
    double m_max;

public:

    m_ticker(double min, double max): m_min(min), m_max(max) {}

protected:

    virtual void generate(const QCPRange &range, const QLocale &locale,
                          QChar formatChar, int precision,
                          QVector<double> &ticks,
                          QVector<double> *subTicks,
                          QVector<QString> *tickLabels) override
    {
        Q_UNUSED(locale)
        QCPAxisTicker::generate(range,locale,formatChar,precision,ticks,subTicks,tickLabels);

        tickLabels->clear();
        int NbTicks = ticks.size();
        double increment = (m_max-m_min)/(NbTicks-1);
        for(int n=0; n<NbTicks; n++)
        {
            double tickValue = m_min+n*increment;
            const QString &tickLabel = QString("%1").arg(tickValue);
            tickLabels->push_back(tickLabel);
        }
    }
};

#endif // XYVIEWER_H
