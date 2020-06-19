#ifndef TIMESTEPBUILDER_H
#define TIMESTEPBUILDER_H

//! ---
//! Qt
//! ---
#include <QDialog>
#include <QVector>

//! ----------------
//! custom includes
//! ----------------
#include "qcustomplotextended.h"

class QAction;
class QMenuBar;
class QGridLayout;
class QVBoxLayout;
class QToolBar;

class TimeStepBuilder: public QDialog
{
    Q_OBJECT

public:

    //TimeStepBuilder(QWidget *parent = 0);
    TimeStepBuilder(const QString &currentDirectory, const QString &outputDir, QWidget *parent = 0);
    ~TimeStepBuilder();

    //void setWorkingDirectory(const QString &workingDirectory) { myWorkingDir = workingDirectory; }
    //QString getWorkingDirectory() const  { return myWorkingDir; }

private:

    QCustomPlotExtended *myPlot;

    QStatusBar *myStatusBar;
    QMenuBar* myMenuBar;
    QMenu *myContextMenu;
    QMenu *menuFile;

    QToolBar *myToolBar;

    QAction *actionOpen;
    QAction *actionClear;
    QAction *actionAccept;
    QAction *actionExit;

    QAction *actionZoom;
    QAction *actionPan;
    QAction *actionFit;

    QAction *actionSaveImage;
    QAction *actionFitRange;

    QString timeHistoryFile;

    QVector<double> timeList;
    QVector<double> valueList;

    QString firstColumnName, secondColumnName;
    QString myInputFile;
    QString myOutputDir;

private:

    void createContent();
    void createActions();
    void createPlot();
    void writeTimeStepFile(const QString &fileName);

private slots:

    void showContextMenu(QPoint pos);
    void updateStatusBar(double time);
    bool openFile();
    void saveImage();
    void fitRange();
    void acceptTimeDivision();
    void pan();
    void zoom();
    void fit();
    void clearPanel();
};

#endif // TIMESTEPBUILDER_H
