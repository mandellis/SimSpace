#ifndef QHISTOGRAM_H
#define QHISTOGRAM_H

//! ----------------
//! custom includes
//! ----------------
#include <qhistogramdata.h>

//! ----
//! C++
//! ----
#include <iostream>

//! ---
//! Qt
//! ---
#include <QWidget>
#include <QDockWidget>
#include <QCloseEvent>

class QCustomPlot;
class QCPBars;

class QHistogram: public QWidget
{
    Q_OBJECT

public:

    //! constructor
    QHistogram(QWidget *parent=0);
    QHistogram(const histogramData &hdata, QWidget *parent=0);

    //! set data
    void setData(const histogramData &hdata, bool updateViewer=true);

public slots:

    //! update plot
    void updatePlot();

private:

    //! internal data
    histogramData myHistogramData;

    //! plot
    QCustomPlot *myPlot;
    QCPBars *myBars;

    //! create content
    void createContent();

    //! clear
    void clear();

    //! make fake data
    histogramData fakeData();

protected:

};

#endif // QHISTOGRAM_H
