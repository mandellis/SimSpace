#ifndef SYSTEMCONSOLE_H
#define SYSTEMCONSOLE_H

//! Qt
#include <QWidget>
#include <QPlainTextEdit>
#include <QMenuBar>

//! C++
#include <iostream>
#include <fstream>
using namespace std;

class QVBoxLayout;

class systemConsole : public QWidget
{
    Q_OBJECT

public:

    //! constructor
    explicit systemConsole(bool isCreatedAsRunning = true,
                           bool isMenuBarVisible = true,
                           bool isContinuosLog = false,
                           const QString &myLogFile = QString(),
                           QWidget *parent = 0);

    //! destructor
    virtual ~systemConsole();

    //! operator <<
    void operator <<(const QString &aText)
    {
        myTextEdit->appendPlainText(aText);
    }

    //! set log file name
    virtual void setLogFile(const QString& logFilePathAbsolute)
    {
        myLogFile = logFilePathAbsolute;
    }

private:

    QVBoxLayout *vl;
    QMenuBar *menuBar;
    QFont myFont;

    QAction *actionStart;
    QAction *actionStop;
    QAction *actionSave;
    QAction *actionSaveAs;
    QAction *actionClear;
    QAction *actionExit;


    bool isCreatedAsRunning;
    bool myIsMenuBarVisible;

protected:

    QString myLogFile;
    QPlainTextEdit *myTextEdit;
    bool isRunning;
    bool myContinuoslyLogging;

    ofstream *myOutFileStream;

private:

    void createContent();

protected:

    virtual bool eventFilter(QObject* object, QEvent* event);

public:

    bool IsRunning() {return isRunning;}
    void setContinuoslyLogging(bool isOn, const QString &logFile);

public slots:

    void start();
    void stop();
    void clear();
    void save();
    void saveAs();

    void setMenuHidden()
    {
        menuBar->hide();
        myIsMenuBarVisible=false;
    }
    void setMenuVisible()
    {
        menuBar->show();
        myIsMenuBarVisible=true;
    }

public slots:

    void showContextMenu(QPoint p);
    void setText(const QString &text);

signals:

    void systemConsoleClosed();
    void systemConsoleStopped();
    void systemConsoleStarted();
};

#endif // SYSTEMCONSOLE_H
