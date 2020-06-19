#ifndef STLAPIWRITER_H
#define STLAPIWRITER_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>

#include <Standard_Boolean.hxx>
#include <StlAPI_ErrorStatus.hxx>
#include <Standard_CString.hxx>

class StlMesh_Mesh;
class TopoDS_Shape;

//! This class creates and writes
//! STL files from Open CASCADE shapes. An STL file can be
//! written to an existing STL file or to a new one..

class STLAPIWriter
{
public:

    DEFINE_STANDARD_ALLOC

    //! Creates a writer object with
    //! default parameters: ASCIIMode, canonical STL
    Standard_EXPORT STLAPIWriter();

    //! Returns the address to the
    //! flag defining the mode for writing the file. This address
    //! may be used to either read or change the flag.
    //! If the mode returns True (default value) the generated
    //! file is an ASCII file. If the mode returns False, the
    //! generated file is a binary file.
    Standard_EXPORT Standard_Boolean& ASCIIMode();

    //! Converts a given shape to STL format and writes it to file with a given filename.
    //! \return the error state.
    Standard_EXPORT bool Write (const TopoDS_Shape& aShape, const Standard_CString aFileName);
    Standard_EXPORT bool WriteExtended(const TopoDS_Shape& aShape, const Standard_CString aFileName);


protected:

private:

    Standard_Boolean theASCIIMode;
    occHandle(StlMesh_Mesh) theStlMesh;
};

#endif // STLAPIWRITER_H
