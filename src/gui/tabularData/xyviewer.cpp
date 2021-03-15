//! ----------------
//! custom includes
//! ----------------
#include "xyviewer.h"
#include "basetabledelegate.h"

//! ---
//! Qt
//! ---
#include <QVariant>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMenu>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! --------
//! STL C++
//! --------
#include <set>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
XYViewer::XYViewer(QWidget *parent):
    QAbstractItemView(parent)
{    
    //! table data
    m_model = nullptr;

    //! series of normalized data
    m_series.clear();

    //! context menu
    m_contextMenu = new QMenu("XY viewer",this);

    //! plotting widget
    m_plot = new QCustomPlot(this);
    m_plot->axisRect()->setupFullAxesBox();
    m_plot->setMinimumSize(500,100);

    //! context menu
    m_plot->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(m_plot,SIGNAL(customContextMenuRequested(const QPoint &)),this,SLOT(showContextMenu(const QPoint&)));

    //! enable multiselection of plots
    //m_plot->setInteractions(QCP::iSelectPlottables);
    //m_plot->setMultiSelectModifier(Qt::ControlModifier);

    //! connections for handling selection in the viewer
    //connect(m_plot,SIGNAL(selectionChangedByUser()),this,SLOT(selectionChangedByUser()));

    //! layout for automatic resize
    QHBoxLayout *h = new QHBoxLayout(this);
    h->setContentsMargins(0,0,0,0);

    //! add plot widget
    h->addWidget(m_plot);

    //! ------------------------------------
    //! table view of controls and delegate
    //! ------------------------------------
    m_tableOfControls = new QTableView(this);                                       // viewer
    m_tableOfControls->setAlternatingRowColors(true);                               // cosmetic
    baseTableDelegate *tableOfControlsDelegate = new baseTableDelegate(this);       // delegate
    m_tableOfControls->setItemDelegate(tableOfControlsDelegate);
    connect(tableOfControlsDelegate,SIGNAL(dataChanged()),this,SLOT(replotAll()));
    m_tableOfControls->hide();                                                      // initially hide the table of controls

    //! --------------------------------
    //! table model containing controls
    //! --------------------------------
    m_controls = new BaseTableModel(this);                                          // model
    m_tableOfControls->setModel(m_controls);                                        // set model
    //connect(m_controls,SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),this,SLOT(replotAll()));
    //connect(m_controls,SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),this,SLOT(dataChanged(QModelIndex,QModelIndex,QVector<int>)));
    h->addWidget(m_tableOfControls);

    //! set the header for columns
    std::vector<QString> labels { "Number", "Visible", "X col", "Y col", "X visible", "Y visible" };
    for(int n=0; n<labels.size(); n++) m_controls->setColumnHeader(n,labels[n]);
}

//! ---------------------------
//! function: buildContextMenu
//! details:
//! ---------------------------
void XYViewer::buildCustomMenu()
{
    m_contextMenu->clear();

    //! -------------------
    //! show/hide controls
    //! -------------------
    if(m_tableOfControls->isVisible())
    {
        QAction *actionHideControls = new QAction("Hide controls", m_contextMenu);
        connect(actionHideControls,SIGNAL(triggered(bool)),this,SLOT(hideControls()));
        m_contextMenu->addAction(actionHideControls);

        m_contextMenu->addSeparator();

        //! ---------
        //! add plot
        //! ---------
        QAction *actionAddPlot = new QAction("Add plot", m_contextMenu);
        m_contextMenu->addAction(actionAddPlot);
        connect(actionAddPlot,SIGNAL(triggered(bool)),this,SLOT(addPlot()));
    }
    else
    {
        QAction *actionShowControls = new QAction("Show controls", m_contextMenu);
        m_contextMenu->addAction(actionShowControls);
        connect(actionShowControls,SIGNAL(triggered(bool)),this,SLOT(showControls()));
    }
}

//! --------------------------
//! function: showContextMenu
//! details:
//! --------------------------
void XYViewer::showContextMenu(const QPoint &pos)
{
   cout<<"XYViewer::showContextMenu()->____function called____"<<endl;

   this->buildCustomMenu();
   m_contextMenu->exec(this->mapToGlobal(pos));
}

//! --------------------------------
//! function: removeAdditionalAxes
//! details:  at bottom and af left
//! --------------------------------
void XYViewer::removeAdditionalAxes()
{
    if(m_plot->axisRect()->axis(QCPAxis::atBottom,0)==nullptr) exit(9999);

    int Nh = m_plot->axisRect()->axisCount(QCPAxis::atBottom);
    while(m_plot->axisRect()->axisCount(QCPAxis::atBottom)!=1)
    {
        QCPAxis* hAxis = m_plot->axisRect()->axis(QCPAxis::atBottom,--Nh);
        m_plot->axisRect()->removeAxis(hAxis);
    }
    int Nv = m_plot->axisRect()->axisCount(QCPAxis::atLeft);
    while(m_plot->axisRect()->axisCount(QCPAxis::atLeft)!=1)
    {
        QCPAxis* vAxis = m_plot->axisRect()->axis(QCPAxis::atLeft,--Nv);
        m_plot->axisRect()->removeAxis(vAxis);
    }
}

//! --------------------
//! function: replotAll
//! details:
//! --------------------
void XYViewer::replotAll()
{
    //cout<<"XYViewer::replotAll()->____function called____"<<endl;

    if(!m_model) return;
    m_plot->clearGraphs();                                  //! remove all graphs
    this->removeAdditionalAxes();                           //! remove additional axes

    //cout<<"(1)____number of axes at bottom: "<<m_plot->axisRect()->axisCount(QCPAxis::atBottom)<<"____"<<endl;
    //cout<<"(1)____number of axes at left: "<<m_plot->axisRect()->axisCount(QCPAxis::atLeft)<<"____"<<endl;

    m_plot->xAxis->setLabel("");
    m_plot->yAxis->setLabel("");

    //! -------------------------------------------
    //! refresh all series (reread the model data)
    //! -------------------------------------------
    this->refreshAllSeries();

    //! -----------------------------------
    //! normalize data and prepare tickers
    //! -----------------------------------
    std::map<int,QSharedPointer<m_ticker>> tickers;
    for(int n=0; n<m_series.size(); n++)
    {
        QVector<double> hSeries = m_series.at(n);
        double min = 1e20, max = -1e20;
        for(int s=0; s<hSeries.size(); s++) { double val = hSeries[s]; if(val<=min) min=val; if(val>=max) max = val; }
        double delta = max-min;
        if(delta==0) delta = fabs(max)*0.10;
        for(int s=0; s<hSeries.size(); s++)
        {
            double x_s = (hSeries.at(s)-min)/delta;
            m_series.at(n).replace(s,x_s);
        }
        QSharedPointer<m_ticker> aTicker (new m_ticker(min,max));
        tickers.insert(std::make_pair(n,aTicker));
    }

    //! -------------------------------------
    //! read settings from table of controls
    //! -------------------------------------
    for(int n=0; n<m_controls->rowCount(); n++)     // cycle over the number of plots
    {
        //! X/Y axis column
        int Xcol = m_controls->dataByRC(n,1).toInt();
        int Ycol = m_controls->dataByRC(n,2).toInt();

        //! add a graph
        QCPGraph *aGraph = m_plot->addGraph();

        //! add data to the graph
        QVector<double> xSeries = m_series.at(Xcol);
        QVector<double> ySeries = m_series.at(Ycol);
        aGraph->setData(xSeries,ySeries);

        //! axis labels
        QString HLabel = m_model->headerData(Xcol,Qt::Horizontal,Qt::DisplayRole).toString();
        QString VLabel = m_model->headerData(Ycol,Qt::Horizontal,Qt::DisplayRole).toString();

        //! add horizontal and vertical axes, and set range and labels
        //! set visibility according to table of controls
        if(n!=0)
        {
            QCPAxis *hAxis = m_plot->axisRect()->addAxis(QCPAxis::atBottom);
            hAxis->setRange(0,1);
            hAxis->setLabel(HLabel);
            hAxis->setTicker(tickers.at(Xcol));

            bool isHAxisVisible = m_controls->dataByRC(n,3).toBool();
            if(isHAxisVisible) hAxis->setVisible(true);
            else hAxis->setVisible(false);

            QCPAxis *vAxis = m_plot->axisRect()->addAxis(QCPAxis::atLeft);
            vAxis->setRange(0,1);
            vAxis->setLabel(VLabel);
            vAxis->setTicker(tickers.at(Ycol));

            bool isVAxisVisible = m_controls->dataByRC(n,4).toBool();
            if(isVAxisVisible) vAxis->setVisible(true);
            else vAxis->setVisible(false);
        }
        else
        {
            m_plot->xAxis->setRange(0,1);
            m_plot->xAxis->setLabel(HLabel);
            m_plot->xAxis->setTicker(tickers.at(Xcol));
            bool isHAxisVisible = m_controls->dataByRC(n,3).toBool();
            if(!isHAxisVisible) m_plot->xAxis->setVisible(false);
            else m_plot->xAxis->setVisible(true);

            m_plot->yAxis->setRange(0,1);
            m_plot->yAxis->setLabel(VLabel);
            m_plot->yAxis->setTicker(tickers.at(Ycol));
            bool isVAxisVisible = m_controls->dataByRC(n,4).toBool();
            if(!isVAxisVisible) m_plot->yAxis->setVisible(false);
            else m_plot->yAxis->setVisible(true);
        }
    }
    //cout<<"(2)____number of axes at bottom: "<<m_plot->axisRect()->axisCount(QCPAxis::atBottom)<<"____"<<endl;
    //cout<<"(2)____number of axes at left: "<<m_plot->axisRect()->axisCount(QCPAxis::atLeft)<<"____"<<endl;

    //! -----------------------------------------------
    //! number of visible horizontal and vertical axes
    //! eliminate duplicated axes
    //! -----------------------------------------------
    int NbHAxes = m_plot->axisRect()->axisCount(QCPAxis::atBottom);
    int NbVAxes = m_plot->axisRect()->axisCount(QCPAxis::atTop);
    int NbHAxesVisible = 0;
    int NbVAxesVisible = 0;
    for(int n=0; n<NbHAxes; n++) if(m_plot->axisRect()->axis(QCPAxis::atBottom,n)->visible()) NbHAxesVisible++;
    for(int n=0; n<NbVAxes; n++) if(m_plot->axisRect()->axis(QCPAxis::atLeft,n)->visible()) NbVAxesVisible++;

    std::set<QString> labels;
    for(int n=0; n<NbHAxesVisible; n++)
    {
        QCPAxis *axis = m_plot->axisRect()->axis(QCPAxis::atBottom,n);
        if(labels.count(axis->label())==0)
        {
            labels.insert(axis->label());
            continue;
        }
        m_plot->axisRect()->axis(QCPAxis::atBottom,n)->setVisible(false);
    }
    labels.clear();
    for(int n=0; n<NbVAxesVisible; n++)
    {
        QCPAxis *axis = m_plot->axisRect()->axis(QCPAxis::atLeft,n);
        if(labels.count(axis->label())==0)
        {
            labels.insert(axis->label());
            continue;
        }
        m_plot->axisRect()->axis(QCPAxis::atLeft,n)->setVisible(false);
    }

    //! -------------------------------------------------
    //! show/hide plots according to visibility in table
    //! -------------------------------------------------
    for(int n=0; n<m_controls->rowCount(); n++)     // cycle over the number of plots
    {
        bool isPlotVisible = m_controls->dataByRC(n,0).toBool();
        m_plot->graph(n)->setVisible(isPlotVisible);
    }

    //! -------
    //! replot
    //! -------
    m_plot->replot();
}

//! ----------------------
//! function: dataChanged
//! details:
//! ----------------------
void XYViewer::dataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight, const QVector<int> &roles)
{
    //cout<<"XYViewer::dataChanged()->____function called____"<<endl;
    Q_UNUSED(topLeft)
    Q_UNUSED(bottomRight)
    if(roles.contains(Qt::EditRole))
    {
        //cout<<"XYViewer::dataChanged()->____edit role____"<<endl;
        this->replotAll();
    }
    QAbstractItemView::dataChanged(topLeft,bottomRight,roles);
}

//! ---------------------------
//! function: refreshAllSeries
//! details:
//! ---------------------------
void XYViewer::refreshAllSeries()
{
    //cout<<"XYViewer::refreshAllSeries()->____function called____"<<endl;

    //! ----------------
    //! check data size
    //! ----------------
    int NbColumns = m_model->columnCount();
    int NbRows = m_model->rowCount();
    if(NbColumns<1 || NbRows==0) return;

    //cout<<"XYViewer::refreshAllSeries()->____(rows, cols) = ("<<NbRows<<", "<<NbColumns<<")____"<<endl;

    //! -------------------
    //! refresh the series
    //! -------------------
    m_series.clear();

    //! ------------------------------------------
    //! read the data from the table model
    //! and create a series for each table column
    //! ------------------------------------------
    for(int col = 0; col<NbColumns; col++)
    {
        QVector<double> aSeries;
        for(int row = 0; row<NbRows; row++)
        {
            //cout<<"____reading col: "<<col<<" row: "<<row;
            QVariant data = m_model->dataByRC(row,col);
            double value;

            //! here you can "translate" into "double" other
            //! types of data, different from "double"

            if(data.isNull()) value = 0;
            if(!data.isValid()) value = 0;
            if(data.canConvert<double>()) value = data.toDouble();
            else value = 0;
            aSeries.push_back(value);
            //cout<<" => done"<<" => "<<value<<endl;
        }
        m_series.insert(std::make_pair(col,aSeries));
        //cout<<"____done____"<<endl;
    }
}

//! -------------------
//! function: setModel
//! details:
//! -------------------
void XYViewer::setModel(BaseTableModel *model)
{
    if(model==nullptr) return;

    //! private member init
    m_model = model;

    //! ------------------------------------------------------------------------
    //! isVisible, x axis column, y axis column, x axis visible, y axis visible
    //! ------------------------------------------------------------------------
    for(int n=1; n<m_model->columnCount(); n++)
    {
        std::vector<QVariant> aRowOfControls { true, 0, n, true, true };
        m_controls->appendRow(aRowOfControls);
    }

    //! ---------------------------------------------
    //! set the header labels for control table rows
    //! ---------------------------------------------
    int NbPlots = m_model->columnCount()-1;
    for(int n=0; n<NbPlots; n++)
    {
        QVariant data;
        data.setValue(QString("Graph # %1").arg(n));
        m_controls->setHeaderData(n,Qt::Vertical,data,Qt::DisplayRole);
    }

    //! ------------------------------------------------
    //! set the header labels for control table columns
    //! ------------------------------------------------
    QVariant data;
    std::vector<QString> labels { "Visible", "H axis", "V axis", "HA visible", "VA visible"};
    for(int n=0; n<labels.size(); n++)
    {
        data.setValue(labels.at(n));
        m_controls->setHeaderData(n,Qt::Horizontal,data,Qt::DisplayRole);
    }

    //! -----------
    //! replot all
    //! -----------
    this->replotAll();

    QAbstractItemView::setModel(model);
}

//! ---------------------------------------------------
//! function: horizontalOffset
//! details:  return the horizontal offset of the view
//! ---------------------------------------------------
int XYViewer::horizontalOffset() const
{
    return 0;
}

//! --------------------------------------------------
//! function: verticalOffset
//! details:  Returns the vertical offset of the view
//! --------------------------------------------------
int XYViewer::verticalOffset() const
{
    return 0;
}

//! --------------------------------------------------------------------------------
//! function: indexAt
//! details:  Returns the model index of the item at the viewport coordinates point
//! --------------------------------------------------------------------------------
QModelIndex XYViewer::indexAt(const QPoint &point) const
{
    Q_UNUSED(point)
    return QModelIndex();
}

//! ----------------------------------------------------------------------
//! function: isIndexHidden
//! details:  Returns true if the item referred to by the given index
//!           is hidden in the view, otherwise returns false.
//!           Hiding is a view specific feature. For example in TableView
//!           a column can be marked as hidden or a row in the TreeView.
//! ----------------------------------------------------------------------
bool XYViewer::isIndexHidden(const QModelIndex &index) const
{
    Q_UNUSED(index)
    return false;
}

//! -----------------------------------------------------------------------------------------
//! function: moveCursor
//! details:  Returns a QModelIndex object pointing to the next object in the view,
//!           based on the given cursorAction and keyboard modifiers specified by modifiers.
//! -----------------------------------------------------------------------------------------
QModelIndex XYViewer::moveCursor(CursorAction cursorAction, Qt::KeyboardModifiers modifiers)
{
    Q_UNUSED(cursorAction)
    Q_UNUSED(modifiers)
    return QModelIndex();
}

//! -------------------------------------------------------------------------------------
//! function: scrollTo
//! details:  Scrolls the view if necessary to ensure that the item at index is visible.
//!           The view will try to position the item according to the given hint.
//! -------------------------------------------------------------------------------------
void XYViewer::scrollTo(const QModelIndex &index, ScrollHint hint)
{
    Q_UNUSED(index)
    Q_UNUSED(hint)
    return;
}

//! ----------------------------------------------------------------------------------------------
//! function: setSelection
//! details:  Applies the selection flags to the items in or touched by the rectangle, rect.
//!           When implementing your own itemview setSelection should call
//!           selectionModel()->select(selection, flags) where selection is either an empty
//!           QModelIndex or a QItemSelection that contains all items that are contained in rect.
//! ----------------------------------------------------------------------------------------------
void XYViewer::setSelection(const QRect &rect, QItemSelectionModel::SelectionFlags flags)
{
    Q_UNUSED(rect)
    Q_UNUSED(flags)
    return;
}

//! --------------------------------------------------------------------------------
//! function: visualRect
//! details:  Returns the rectangle on the viewport occupied by the item at index.
//!           If your item is displayed in several areas then visualRect should
//!           return the primary area that contains index and not the complete area
//!           that index might encompasses, touch or cause drawing.
//! --------------------------------------------------------------------------------
QRect XYViewer::visualRect(const QModelIndex &index) const
{
    Q_UNUSED(index)
    return QRect();
}

//! ------------------------------------------------------------------------------------
//! function: visualRegionForSelection
//! details:  Returns the region from the viewport of the items in the given selection.
//! ------------------------------------------------------------------------------------
QRegion XYViewer::visualRegionForSelection(const QItemSelection &selection) const
{
    Q_UNUSED(selection)
    return QRegion();
}

//! ---------------------------
//! function: viewportSizeHint
//! details:
//! ---------------------------
QSize XYViewer::viewportSizeHint() const
{
    return m_plot->size();
}

//! ----------------------
//! function: paintEvent
//! details:
//! ----------------------
void XYViewer::paintEvent(QPaintEvent *event)
{
    QAbstractItemView::paintEvent(event);
}

//! ------------------
//! function: hideAll
//! details:
//! ------------------
void XYViewer::hideAll()
{
    m_plot->clearGraphs();
    this->removeAdditionalAxes();
    m_plot->xAxis->setRange(0,1);
    m_plot->yAxis->setRange(0,1);
    m_plot->xAxis->setLabel("");
    m_plot->yAxis->setLabel("");
    m_plot->replot();
}

//! -----------------------
//! function: hideSelected
//! details:
//! -----------------------
void XYViewer::hideSelected()
{
    QCPGraph *graph = (QCPGraph*)m_selectedPlottable;
    graph->setVisible(false);
    int index = m_mapOfGraphs.find(graph)->second;
    m_plot->axisRect()->axis(QCPAxis::atLeft,index)->setVisible(false);
    m_plot->axisRect()->axis(QCPAxis::atBottom,index)->setVisible(false);
    m_plot->replot();
    m_selectedPlottable=nullptr;
}

//! ------------------
//! function: addPlot
//! details:
//! ------------------
void XYViewer::addPlot()
{
    //cout<<"XYViewer::addPlot()->____function called____"<<endl;
    int NbRows = m_controls->rowCount();
    std::vector<QVariant> aRowOfControls { true, 0, m_controls->rowCount()-1, true, true };
    m_controls->appendRow(aRowOfControls);
    m_controls->setRowHeader(m_controls->rowCount()-1,QString("Graph # %1").arg(NbRows));
    this->replotAll();
}

//! -------------------
//! function: addGraph
//! details:
//! -------------------
bool XYViewer::addGraph(int xDataCol, int yDataCol)
{
    if(xDataCol>m_model->columnCount()-1 || yDataCol>m_model->columnCount()-1) return false;
    this->replotAll();
    return true;
}

//! ---------------------------------
//! function: selectionChangedByUser
//! details:
//! ---------------------------------
void XYViewer::selectionChangedByUser()
{
    m_selectedPlottable = m_plot->selectedPlottables().at(0);
    if(m_selectedPlottable==nullptr) cout<<"____NULL____"<<endl;
}

//! -----------------------
//! function: hideControls
//! details:
//! -----------------------
void XYViewer::hideControls()
{
    m_tableOfControls->hide();
}

//! -----------------------
//! function: showControls
//! details:
//! -----------------------
void XYViewer::showControls()
{
    m_tableOfControls->show();
}
