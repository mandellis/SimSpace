#ifndef TIMESTEPBUILDER_H
#define TIMESTEPBUILDER_H

//! Qt
#include <QDialog>
#include <QVector>

class QAction;
class QMenuBar;
class QMenu;
class QGridLayout;
class QVBoxLayout;
class QToolBar;
class QCustomPlotExtended;
class QStatusBar;

class TimeStepBuilder: public QDialog
{
    Q_OBJECT

public:

    TimeStepBuilder(QWidget *parent = 0, Qt::WindowFlags f =Qt::Dialog);
    TimeStepBuilder(const QString &currentDirectory, QWidget *parent = 0, Qt::WindowFlags f =Qt::Dialog);
    ~TimeStepBuilder();

    void setWorkingDirectory(const QString &workingDirectory) { myWorkingDir = workingDirectory; }
    QString getWorkingDirectory() const  { return myWorkingDir; }

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
    QString myWorkingDir;

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

protected:

};

#endif // TIMESTEPBUILDER_H
