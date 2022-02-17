#include "steptools.h"
#include <XSControl_WorkSession.hxx>
#include <Interface_InterfaceModel.hxx>
#include <StepRepr_Representation.hxx>

//! -------------------------
//! function: dumpSTEPLabels
//! details:
//! -------------------------
void STEPTools::dumpSTEPLabels(const STEPControl_Reader &reader, TColStd_ListOfAsciiString &listOfNames)
{
    const occHandle(XSControl_WorkSession)& theSession = reader.WS();
    const occHandle(Interface_InterfaceModel)& theModel = theSession->Model();
    Standard_Integer nb = theModel->NbEntities();
    for(Standard_Integer i=1; i<=nb; i++)
    {
        occHandle(StepRepr_Representation) ent = occHandle(StepRepr_Representation)::DownCast(theModel->Value(i));
        if (ent.IsNull()) continue;
        if (ent->Name().IsNull()) continue;
        //! --------------
        //! HAscii string
        //! --------------
        const occHandle(TCollection_HAsciiString) &theHName = ent->Name();

        //! -------------
        //! Ascii string
        //! -------------
        const TCollection_AsciiString &theName = theHName->String();
        cout<<"STEPTools::dumpSTEPLabels::_____name "<<theName.ToCString()<<endl;
        listOfNames.Append(theName);
    }
}
