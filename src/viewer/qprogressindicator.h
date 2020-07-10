#ifndef QPROGRESSINDICATOR_H
#define QPROGRESSINDICATOR_H

//! -------------------------------------------
//! redefinition of opencascade Handle() macro
//! -------------------------------------------
#include "occhandle.h"

//! ---
//! Qt
//! ---
#include <QWidget>
#include <QString>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

class QProgressBar;
class QPushButton;
class QLabel;
class QThread;

class QProgressIndicator: public QWidget
{
    Q_OBJECT

public:

    //! constructor
    QProgressIndicator(bool additionalBar=false, QWidget *parent=0);

    //! destructor
    QProgressIndicator::~QProgressIndicator();

protected:

    //! ------------------------------
    //! create the content with style
    //! 0 => vertical layout
    //! 1 => horizontal layout
    //! ------------------------------
    virtual void createContent(bool additionalBar=true);

    //! event filter
    virtual bool eventFilter(QObject *object, QEvent *event);

protected:

    QPushButton *myStopButton;
    QProgressBar *myProgressBar;
    QProgressBar *myProgressBar1;
    QLabel *myLabel;

private:

    QObject *myCurCallerObject;
    QString myCurrentRunningTask;

public:

    void setPrimaryBarVisible(bool isVisible);
    void setSecondaryBarVisible(bool isVisible);
    void setAbortButtonVisible(bool isVisible);

    int getCurrentProgress();
    int getCurrentProgress1();
    QString getText();
    void setText(const QString &text);
    void setTask(const QString &task);

public slots:

    void enableStop();
    void disableStop();
    void handleStopPressed();

signals:

    void abortNetgenSTLPressed();
    void requestStopCCX();
};

#endif // QPROGRESSINDICATOR_H
