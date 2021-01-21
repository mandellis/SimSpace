#ifndef EXTENDEDRWSTL_H
#define EXTENDEDRWSTL_H

// Created on: 1994-10-13
// Created by: Marc LEGAY
// Copyright (c) 1994-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <Standard.hxx>
#include <Standard_Handle.hxx>
#include <Message_ProgressIndicator.hxx>
#include <NCollection_Array1.hxx>
#include <QList>
#include <QString>
#include <QMap>
#include <vector>

class StlMesh_Mesh;
class OSD_Path;

class ExtendedRWStl
{
public:

    /*
    Standard_EXPORT static occHandle(StlMesh_Mesh) ReadExtendedSTLAscii (const QString& thePath,
                                                                         NCollection_Array1<occHandle(StlMesh_Mesh)> &vecFaceStlMesh_Mesh,
                                                                         QList<QMap<int,int>> &maps_localToGlobal_nodeIDs,
                                                                         QList<QMap<int,int>> &maps_localToGlobal_elementIDs);
    */
    Standard_EXPORT static occHandle(StlMesh_Mesh) ReadExtendedSTLAscii (const QString& thePath,
                                                                         std::vector<occHandle(StlMesh_Mesh)> &vecFaceStlMesh_Mesh,
                                                                         QList<QMap<int,int>> &maps_localToGlobal_nodeIDs,
                                                                         QList<QMap<int,int>> &maps_localToGlobal_elementIDs);


    Standard_EXPORT static occHandle(StlMesh_Mesh) ReadAscii(const OSD_Path& aPath,
                                                             QMap<int,int> &triangleTagMap,
                                                             const occHandle(Message_ProgressIndicator)& aProgInd = NULL);


    Standard_EXPORT static void WriteAscii (const occHandle(StlMesh_Mesh)& theMesh,
                                            const QMap<int,int> &triangleTagMap,
                                            const OSD_Path& thePath);

};

#endif // EXTENDEDRWSTL_H
