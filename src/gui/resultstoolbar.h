#ifndef RESULTSTOOLBAR_H
#define RESULTSTOOLBAR_H

//! ---
//! Qt
//! ---
#include <QToolBar>
#include <QWidget>

//! ----
//! C++
//! ----
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

    //! --------------------------------
    //! type of presentation
    //! menu for the button and actions
    //! --------------------------------
    QPushButtonExtended *typeOfPresentationButton;
    QMenu *typeOfPresentationMenu;

    QAction *actionUseIsoStrips;
    QAction *actionUseIsoSurfaces;
    QAction *actionUseSmoothNodal;
    QAction *actionUseIsoLines;

    //! --------------------------------
    //! combined views
    //! menu for the button and actions
    //! --------------------------------
    QPushButtonExtended *combinedViewSelectorButton;
    QMenu *combinedViewSelectorMenu;
    QAction *actionNoWireframe;
    QAction *actionShowUndeformedWireframe;
    QAction *actionShowUndeformedModel;
    QAction *actionShowElements;

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
    void emitRequestUpdatePostObjectScale(double scale);

    void emitRequstUseIsoStrips();
    void emitRequestUseIsoSurface();
    void emitRequestUseSmoothNodal();
    void emitRequestUseIsoLines();

    void updateIcon(QAction *action);

signals:

    void requestUpdateViewerStatus();
    void requestUpdatePostObjectScale(double scale);

public:


};

#endif // RESULTSTOOLBAR_H
