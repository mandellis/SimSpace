#ifndef RESULTSTOOLBAR_H
#define RESULTSTOOLBAR_H

//! Qt
#include <QToolBar>
#include <QWidget>

//! C++
#include <iostream>
using namespace std;

//! ----------------
//! custom includes
//! ----------------
#include "qpushbuttonextended.h"
#include "resultpresentation.h"

class QMenu;
class QAction;

class ResultsToolBar: public QToolBar
{
    Q_OBJECT

public:

    ResultsToolBar(const QString &title, QWidget *parent=0);

private:

    //! --------
    //! actions
    //! --------
    QAction *actionNoWireframe;
    QAction *actionShowUndeformedWireframe;
    QAction *actionShowUndeformedModel;
    QAction *actionShowElements;

    //! -----------------
    //! button with menu
    //! -----------------
    QPushButtonExtended *combinedViewSelectorButton;

    //! --------------------
    //! menu for the button
    //! --------------------
    QMenu *combinedViewSelectorMenu;

    //! status variable
    //resultPresentation myResultPresentation;

    //! ------------
    //! save status
    //! ------------
    bool saveStatus(const std::string &fileName);

    //! -----------
    //! set status
    //! -----------
    void setStatus(const resultPresentation &aResultPresentation);

private slots:

    void emitRequestNoWireframe();
    void emitRequestShowUndeformedModel();
    void emitRequestShowUndeformedWireframe();
    void emitRequestShowMeshElements();
    void ResultsToolBar::emitRequestUpdatePostObjectScale(double scale);

    void updateIcon(QAction *action);

signals:

    void requestUpdateViewerStatus();
    void requestUpdatePostObjectScale(double scale);

public:


};

#endif // RESULTSTOOLBAR_H
