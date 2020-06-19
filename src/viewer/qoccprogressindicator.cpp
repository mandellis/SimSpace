//! ----------------
//! custom includes
//! ----------------
#include "qoccprogressindicator.h"

//! ---
//! Qt
//! ---
#include <QPushButton>
#include <QApplication>
#include <QProgressDialog>

IMPLEMENT_STANDARD_HANDLE(QOccProgressIndicator,Message_ProgressIndicator)
IMPLEMENT_STANDARD_RTTIEXT(QOccProgressIndicator,Message_ProgressIndicator)

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
QOccProgressIndicator::QOccProgressIndicator(int min, int max, QWidget *parent, Qt::WindowFlags theFlags)
{
    myProgress = new QProgressDialog (parent, theFlags);
    myProgress->setWindowModality (Qt::ApplicationModal);
    myProgress->setFixedWidth(300);
    myProgress->setCancelButtonText("Stop");

    //! --------------------------
    //! style of the progress bar
    //! --------------------------
    QString progressBarStyleSheet = "QProgressBar"
                                    "{"
                                    "border: 1px solid grey;"
                                    "border-radius: 5px;"
                                    "text-align: center;"
                                    "}";

    myProgress->setStyleSheet(progressBarStyleSheet);

    //! set title
    myProgress->setWindowTitle("Progress indicator");

    //! set the range for Message_ProgressIndicator
    this->SetScale(min,max,1);

    //! synch with Qt
    myProgress->setRange(min,max);

    //! the dialog will appear only after 500 ms
    myProgress->setMinimumDuration(500);
}

//! ---------------
//! function: show
//! details:
//! ---------------
Standard_Boolean QOccProgressIndicator::Show(const Standard_Boolean theForce)
{
    Q_UNUSED(theForce)

    //! update the progress:
    //! [1] get the position (double, always within [0,1]) from the Message_ProgressIndicator)
    Standard_Real position = GetPosition();
    int theVal = int(myProgress->minimum()+ position *(myProgress->maximum()-myProgress->minimum()));
    myProgress->setValue(theVal);

    //! [2] set the label
    occHandle(TCollection_HAsciiString) aName = GetScope(1).GetName();
    if(!aName.IsNull()) myProgress->setLabelText (aName->ToCString());

    //! [3] to let redraw and keep GUI responsive
    QApplication::processEvents();

    return Standard_True;
}

//! --------------------
//! function: userBreak
//! details:
//! --------------------
Standard_Boolean QOccProgressIndicator::UserBreak()
{
    Standard_Boolean wasCanceled = myProgress->wasCanceled();
    if(wasCanceled)
    {
        cout<<"QOccProgressIndicator::UserBreak()->____Stop request has been received____"<<endl;
        myProgress->setLabelText("Your stop request has been received");
    }
    return wasCanceled;
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
QOccProgressIndicator::~QOccProgressIndicator()
{
    if (myProgress)
    {
        delete myProgress;
        myProgress = Q_NULLPTR;
    }
}
