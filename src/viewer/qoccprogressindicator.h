#ifndef QOCCPROGRESSINDICATOR_H
#define QOCCPROGRESSINDICATOR_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

//! Qt
#include <QObject>

//! OCC
#include <Message_ProgressIndicator.hxx>
#include <Standard_DefineHandle.hxx>

class QProgressDialog;

class QOccProgressIndicator: public Message_ProgressIndicator
{

public:

    DEFINE_STANDARD_RTTIEXT(QOccProgressIndicator,Message_ProgressIndicator)

    //! constructor
    Standard_EXPORT QOccProgressIndicator(int min=0, int max=100, QWidget *parent = 0, Qt::WindowFlags theFlags = 0);

    //! Deletes the object
    Standard_EXPORT virtual ~QOccProgressIndicator ();

    //! Updates presentation of the object
    Standard_EXPORT virtual bool Show(const Standard_Boolean theForce = Standard_True);

    //! Returns True if the user has signaled to cancel the process
    Standard_EXPORT virtual Standard_Boolean UserBreak();

protected:

    QProgressDialog *myProgress;
};

DEFINE_STANDARD_HANDLE(QOccProgressIndicator,Message_ProgressIndicator)

#endif // QOCCPROGRESSINDICATOR_H
