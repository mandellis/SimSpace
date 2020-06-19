#include "extendedrwstl.h"

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

#include <BRepBuilderAPI_CellFilter.hxx>
#include <BRepBuilderAPI_VertexInspector.hxx>
#include <gp.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Message_ProgressIndicator.hxx>
#include <Message_ProgressSentry.hxx>
#include <OSD.hxx>
#include <OSD_File.hxx>
#include <OSD_Host.hxx>
#include <OSD_OpenFile.hxx>
#include <OSD_Path.hxx>
#include <OSD_Protection.hxx>
#include <Precision.hxx>
#include <RWStl.hxx>
#include <Standard_NoMoreObject.hxx>
#include <Standard_TypeMismatch.hxx>
#include <StlMesh_Mesh.hxx>
#include <StlMesh_MeshExplorer.hxx>
#include <TCollection_AsciiString.hxx>
#include <NCollection_Array1.hxx>
#include <stdio.h>

#include <vector>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <QMap>

//! ------------------------------------------------------------------------------
//! A static method adding nodes to a mesh and keeping coincident (sharing) nodes
//! ------------------------------------------------------------------------------
static Standard_Integer AddVertex(occHandle(StlMesh_Mesh)& mesh,
                                  BRepBuilderAPI_CellFilter& filter,
                                  BRepBuilderAPI_VertexInspector& inspector,
                                  const gp_XYZ& p)
{
    Standard_Integer index;
    inspector.SetCurrent(p);
    gp_XYZ minp = inspector.Shift(p, -Precision::Confusion());
    gp_XYZ maxp = inspector.Shift(p, +Precision::Confusion());
    filter.Inspect(minp, maxp, inspector);
    const TColStd_ListOfInteger& indices = inspector.ResInd();
    if (indices.IsEmpty() == Standard_False)
    {
        index = indices.First(); // it should be only one
        inspector.ClearResList();
    }
    else
    {
        index = mesh->AddVertex(p.X(), p.Y(), p.Z());
        filter.Add(index, p);
        inspector.Add(p);
    }
    return index;
}

//! ----------
//! constants
//! ----------
static const size_t ASCII_LINES_PER_FACET =   8;    //! one more line because of triangle tag number
static const int IND_THRESHOLD = 1000;              //! increment the indicator every 1k triangles

//! ----------------------------------------------------------------------------------------
//! function: ReadExtendedSTLAscii
//! details:  read and extended .stl file from disk, and return
//!           - the StlMesh_Mesh of the surface
//!           - the StlMesh_Mesh of each face
//!           - for each face two maps: 1) local to global nodeID for the face nodes
//!                                     2) local to global elementID for the face elements
//!
//!           The Ng_MeshVS_DataSource constructor will be in the form:
//!           Ng_MeshVS_DataSource(const StlMesh_Mesh &aMesh,
//!                                const QList<int,QMap<int,int>> maps_localToGlobalNodeID,
//!                                int faceNr)
//! ----------------------------------------------------------------------------------------
/*
occHandle(StlMesh_Mesh) ExtendedRWStl::ReadExtendedSTLAscii (const QString& thePath,
                                                             NCollection_Array1<occHandle(StlMesh_Mesh)> &vecFaceStlMesh_Mesh,
                                                             QList<QMap<int,int>> &maps_localToGlobal_nodeIDs,
                                                             QList<QMap<int,int>> &maps_localToGlobal_elementIDs)
{
    //cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->____function called on .stl mesh file "<<thePath.toStdString()<<"____"<<endl;
    long ipos;
    Standard_Integer nbLines = 0;
    Standard_Integer nbTris = 0;
    Standard_Integer iTri;
    Standard_Integer i1,i2,i3;
    occHandle(StlMesh_Mesh) ReadMesh;

    //! --------------
    //! Open the file
    //! --------------
    FILE* file = fopen(thePath.toStdString().c_str(),"r");
    if(file==NULL)
    {
        cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->____error in reading the file____"<<endl;
        return occHandle(StlMesh_Mesh)();
    }
    fseek(file,0L,SEEK_END);
    long filesize = ftell(file);
    rewind(file);

    //! --------------------------
    //! count the number of lines
    //! --------------------------
    for (ipos = 0; ipos < filesize; ++ipos)
    {
        if(getc(file)=='\n') nbLines++;
    }

    //! ----------------------------
    //! compute number of triangles
    //! ----------------------------
    nbTris = (nbLines / ASCII_LINES_PER_FACET);

    //! -------------------------------------
    //! go back to the beginning of the file
    //! -------------------------------------
    rewind(file);

    //! ----------------
    //! skip the header
    //! ----------------
    while (getc(file) != '\n');

    //! -----------------------------
    //! the overall surface stl mesh
    //! -----------------------------
    ReadMesh = new StlMesh_Mesh();
    ReadMesh->AddDomain();

    //! ----------------------------------------------------
    //! fill the arrays of maps
    //! total number of faces = (faces in the geometry + 1)
    //! ----------------------------------------------------
    for(int i=0; i<=vecFaceStlMesh_Mesh.Size()-1;i++)
    {
        vecFaceStlMesh_Mesh.Value(i)->AddDomain();
        //! an empty map
        QMap<int,int> aMap;
        maps_localToGlobal_elementIDs.append(aMap);
        maps_localToGlobal_nodeIDs.append(aMap);
    }

    //! ------------------------------------------------------
    //! Filter unique vertices to share the nodes of the mesh
    //! ------------------------------------------------------
    BRepBuilderAPI_CellFilter uniqueVertices(Precision::Confusion());
    BRepBuilderAPI_VertexInspector inspector(Precision::Confusion());

    std::vector<BRepBuilderAPI_CellFilter> vecCellFilter;
    std::vector<BRepBuilderAPI_VertexInspector> vecVertexInspector;

    //! ----------------------------------------------------------
    //! identical technique for each face, including the fake "0"
    //! ----------------------------------------------------------
    for(int i=0; i<=vecFaceStlMesh_Mesh.Size()-1;i++)
    {
        BRepBuilderAPI_CellFilter anUniqueVertices(Precision::Confusion());
        BRepBuilderAPI_VertexInspector anInspector(Precision::Confusion());
        vecCellFilter.push_back(anUniqueVertices);
        vecVertexInspector.push_back(anInspector);
    }

    //! --------------
    //! main reading
    //! --------------
    for (iTri = 0; iTri < nbTris;)
    {
        char x[256]="", y[256]="", z[256]="";

        //! ----------------------------------
        //! reading the facet normal
        //! error should be properly reported
        //! ----------------------------------
        if (3 != fscanf(file,"%*s %*s %80s %80s %80s\n", x, y, z)) break;
        gp_XYZ aN (Atof(x), Atof(y), Atof(z));

        //! -------------------------------
        //! skip the keywords "outer loop"
        //! -------------------------------
        if (fscanf(file,"%*s %*s") < 0) break;

        //! -----------------------------------
        //! reading vertex
        //! errors should be properly reported
        //! -----------------------------------
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z)) break;
        gp_XYZ aV1 (Atof(x), Atof(y), Atof(z));
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z)) break;
        gp_XYZ aV2 (Atof(x), Atof(y), Atof(z));
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z)) break;
        gp_XYZ aV3 (Atof(x), Atof(y), Atof(z));

        //! ---------------------------------------------------------------
        //! here the facet must be built and put in the mesh datastructure
        //! ---------------------------------------------------------------
        i1 = AddVertex(ReadMesh, uniqueVertices, inspector, aV1);
        i2 = AddVertex(ReadMesh, uniqueVertices, inspector, aV2);
        i3 = AddVertex(ReadMesh, uniqueVertices, inspector, aV3);

        //! global element number: number of the element within the overall surface mesh
        int globalElementNumber = ReadMesh->AddTriangle (i1, i2, i3, aN.X(), aN.Y(), aN.Z());

        //! ----------------------------
        //! skip the keywords "endloop"
        //! ----------------------------
        if (fscanf(file,"%*s") < 0) break;

        //! --------------------------------------------
        //! read the tag
        //! ii1, ii2, ii3 are local (face) node numbers
        //! --------------------------------------------
        unsigned int tag;
        if(fscanf(file,"%u",&tag)==1)
        {
            //cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->____tag (face number): "<<tag<<"____"<<endl;

            int ii1 = AddVertex(vecFaceStlMesh_Mesh(tag), vecCellFilter.at(tag), vecVertexInspector.at(tag), aV1);
            int ii2 = AddVertex(vecFaceStlMesh_Mesh(tag), vecCellFilter.at(tag), vecVertexInspector.at(tag), aV2);
            int ii3 = AddVertex(vecFaceStlMesh_Mesh(tag), vecCellFilter.at(tag), vecVertexInspector.at(tag), aV3);

            cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->local ____("<<ii1<<", "<<ii2<<", "<<ii3<<")____"<<endl;
            cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->global____("<<i1<<", "<<i2<<", "<<i3<<")____"<<endl;

            const occHandle(StlMesh_Mesh) &faceSTLMesh = vecFaceStlMesh_Mesh(tag);
            int localElementNumber = faceSTLMesh->AddTriangle(ii1, ii2, ii3, aN.X(), aN.Y(), aN.Z());

            const QMap<int,int> &localToGlobal_nodeID = maps_localToGlobal_nodeIDs.value(tag);
            if(!localToGlobal_nodeID.contains(ii1)) maps_localToGlobal_nodeIDs.operator [](tag).insert(ii1,i1);
            if(!localToGlobal_nodeID.contains(ii2)) maps_localToGlobal_nodeIDs.operator [](tag).insert(ii2,i2);
            if(!localToGlobal_nodeID.contains(ii3)) maps_localToGlobal_nodeIDs.operator [](tag).insert(ii3,i3);

            maps_localToGlobal_elementIDs.operator [](tag).insert(localElementNumber,globalElementNumber);
        }

        //! -----------------------------
        //! skip the keywords "endfacet"
        //! -----------------------------
        if (fscanf(file,"%*s") < 0) break;
    }

    //! -------------------------------------------------------------
    //! diagnostic: print out the face numbers without triangulation
    //! -------------------------------------------------------------
    //unsigned int NbInvalidFaces = 0;
    //for(int tag=0; tag<vecFaceStlMesh_Mesh.Size(); tag++)
    //{
    //    if(vecFaceStlMesh_Mesh.Value(tag)->NbTriangles()==0)
    //    {
    //        NbInvalidFaces++;
    //        //!cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->____The face Nr: "<<tag<<" ("<<NbInvalidFaces<<"-th) has not a triangulation____"<<endl;
    //    }
    //}

    fclose(file);
    return ReadMesh;
}
*/

//! -------------------------------------
//! function: ReadAscii
//! details:  read an extended .stl file
//! -------------------------------------
occHandle(StlMesh_Mesh) ExtendedRWStl::ReadAscii(const OSD_Path& thePath,
                                                 QMap<int,int> &triangleTagMap,
                                                 const occHandle(Message_ProgressIndicator)& theProgInd)
{
    TCollection_AsciiString filename;
    long ipos;
    Standard_Integer nbLines = 0;
    Standard_Integer nbTris = 0;
    Standard_Integer iTri;
    Standard_Integer i1,i2,i3;
    occHandle(StlMesh_Mesh) ReadMesh;

    thePath.SystemName (filename);

    //! Open the file
    FILE* file = OSD_OpenFile(filename.ToCString(),"r");

    fseek(file,0L,SEEK_END);

    long filesize = ftell(file);

    rewind(file);

    //! count the number of lines
    for (ipos = 0; ipos < filesize; ++ipos)
    {
        if (getc(file) == '\n') nbLines++;
    }

    //! compute number of triangles
    nbTris = (nbLines / ASCII_LINES_PER_FACET);

    //! go back to the beginning of the file
    rewind(file);

    //! skip header
    while (getc(file) != '\n');
    ReadMesh = new StlMesh_Mesh();
    ReadMesh->AddDomain();

    //! Filter unique vertices to share the nodes of the mesh.
    BRepBuilderAPI_CellFilter uniqueVertices(Precision::Confusion());
    BRepBuilderAPI_VertexInspector inspector(Precision::Confusion());

    //! main reading
    Message_ProgressSentry aPS (theProgInd, "Triangles", 0, (nbTris - 1) * 1.0 / IND_THRESHOLD, 1);
    for (iTri = 0; iTri < nbTris && aPS.More();)
    {
        char x[256]="", y[256]="", z[256]="";

        //! reading the facet normal
        if (3 != fscanf(file,"%*s %*s %80s %80s %80s\n", x, y, z))
            break; // error should be properly reported
        gp_XYZ aN (Atof(x), Atof(y), Atof(z));

        //! skip the keywords "outer loop"
        if (fscanf(file,"%*s %*s") < 0)
            break;

        //! reading vertex
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z))
            break; // error should be properly reported
        gp_XYZ aV1 (Atof(x), Atof(y), Atof(z));
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z))
            break; // error should be properly reported
        gp_XYZ aV2 (Atof(x), Atof(y), Atof(z));
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z))
            break; // error should be properly reported
        gp_XYZ aV3 (Atof(x), Atof(y), Atof(z));

        //! here the facet must be built and put in the mesh datastruc§§§§§§§§§ture
        //! ---------------------------------------------------------------
        i1 = AddVertex(ReadMesh, uniqueVertices, inspector, aV1);
        i2 = AddVertex(ReadMesh, uniqueVertices, inspector, aV2);
        i3 = AddVertex(ReadMesh, uniqueVertices, inspector, aV3);
        int triangleRank = ReadMesh->AddTriangle (i1, i2, i3, aN.X(), aN.Y(), aN.Z());

        // skip the keywords "endloop"
        if (fscanf(file,"%*s") < 0)
            break;

        //! ---------------------------------
        //! build the map (triangleRank,tag)
        //! ---------------------------------
        unsigned int tag;
        if(fscanf(file,"%u",&tag)==1)
        {
            triangleTagMap.insert(triangleRank,tag);
            //cout<<"____Map(triangleRank,Tag) = ("<<triangleRank<<", "<<tag<<")____"<<endl;
        }

        // skip the keywords "endfacet"
        if (fscanf(file,"%*s") < 0)
            break;

        // update progress only per 1k triangles
        if (++iTri % IND_THRESHOLD == 0)
            aPS.Next();
    }
    fclose(file);
    return ReadMesh;
}

//! ---------------------------------------
//! function: WriteAscii
//! details:  it works only for one domain
//! ---------------------------------------
void ExtendedRWStl::WriteAscii(const occHandle(StlMesh_Mesh)& theMesh,
                               const QMap<int,int> &triangleTagMap,
                               const OSD_Path& thePath)
{
    OSD_File theFile (thePath);
    theFile.Build(OSD_WriteOnly,OSD_Protection());
    TCollection_AsciiString buf ("solid\n");
    theFile.Write (buf,buf.Length());buf.Clear();

    Standard_Real x1, y1, z1;
    Standard_Real x2, y2, z2;
    Standard_Real x3, y3, z3;
    int count = 0;
    char sval[512];

    StlMesh_MeshExplorer aMexp (theMesh);
    for (aMexp.InitTriangle (); aMexp.MoreTriangle(); aMexp.NextTriangle())
    {
        count++;
        int tag = triangleTagMap.value(count);
        aMexp.TriangleVertices (x1,y1,z1,x2,y2,z2,x3,y3,z3);

        gp_XYZ Vect12 ((x2-x1), (y2-y1), (z2-z1));
        gp_XYZ Vect23 ((x3-x2), (y3-y2), (z3-z2));
        gp_XYZ Vnorm = Vect12 ^ Vect23;
        Standard_Real Vmodul = Vnorm.Modulus ();
        if (Vmodul > gp::Resolution())
        {
            Vnorm.Divide (Vmodul);
        }
        else
        {
            // si Vnorm est quasi-nul, on le charge a 0 explicitement
            Vnorm.SetCoord (0., 0., 0.);
        }
        sprintf (sval,
                 " facet normal %.9e %.9e %.9e\n"
                 "   outer loop\n"
                 "     vertex %.9e %.9e %.9e\n"
                 "     vertex %.9e %.9e %.9e\n"
                 "     vertex %.9e %.9e %.9e\n"
                 "   endloop\n"
                 "   %d\n"
                 " endfacet\n",
                 Vnorm.X(), Vnorm.Y(), Vnorm.Z(),
                 x1, y1, z1,
                 x2, y2, z2,
                 x3, y3, z3,
                 tag);
        buf += sval;
        theFile.Write (buf, buf.Length()); buf.Clear();
    }

    buf += "endsolid\n";
    theFile.Write (buf, buf.Length()); buf.Clear();
    theFile.Close();
}


//! ---------------------------------------------------------------------------------------------
//! function: ReadExtendedSTLAscii
//! details:  read and extended .stl file from disk, and return the StlMesh_Mesh of the surface,
//!           the StlMesh_Mesh of each face, for each face two maps: 1) local to global nodeID
//!           for the face nodes 2) local to global elementID for the face elements
//!           The mesh data source constructor will use these info, being called as:
//!           Ng_MeshVS_DataSource(const StlMesh_Mesh &aMesh,
//!                                const QList<int,QMap<int,int>> maps_localToGlobalNodeID,
//!                                int faceNr)
//! ---------------------------------------------------------------------------------------------
occHandle(StlMesh_Mesh) ExtendedRWStl::ReadExtendedSTLAscii (const QString& thePath,
                                                             std::vector<occHandle(StlMesh_Mesh)> &vecFaceStlMesh_Mesh,
                                                             QList<QMap<int,int>> &maps_localToGlobal_nodeIDs,
                                                             QList<QMap<int,int>> &maps_localToGlobal_elementIDs)
{
    long ipos;
    int nbLines = 0;
    int nbTris = 0;
    int i1,i2,i3;
    occHandle(StlMesh_Mesh) ReadMesh;

    //! --------------
    //! Open the file
    //! --------------
    FILE* file = fopen(thePath.toStdString().c_str(),"r");
    if(file==NULL)
    {
        cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->____error in reading the file____"<<endl;
        return occHandle(StlMesh_Mesh)();
    }
    fseek(file,0L,SEEK_END);
    long filesize = ftell(file);
    rewind(file);

    //! --------------------------
    //! count the number of lines
    //! --------------------------
    for (ipos = 0; ipos < filesize; ++ipos) if(getc(file)=='\n') nbLines++;

    //! ----------------------------
    //! compute number of triangles
    //! ----------------------------
    nbTris = (nbLines / ASCII_LINES_PER_FACET);

    //! -------------------------------------
    //! go back to the beginning of the file
    //! -------------------------------------
    rewind(file);

    //! ----------------
    //! skip the header
    //! ----------------
    while (getc(file) != '\n');

    //! -----------------------------
    //! the overall surface stl mesh
    //! -----------------------------
    ReadMesh = new StlMesh_Mesh();
    ReadMesh->AddDomain();

    //! ----------------------------------------------------
    //! fill the arrays of maps
    //! total number of faces = (faces in the geometry + 1)
    //! ----------------------------------------------------
    for(int i=0; i<vecFaceStlMesh_Mesh.size();i++)
    {
        vecFaceStlMesh_Mesh.at(i)->AddDomain();
        //! an empty map
        QMap<int,int> aMap;
        maps_localToGlobal_elementIDs.append(aMap);
        maps_localToGlobal_nodeIDs.append(aMap);
    }

    //! ------------------------------------------------------
    //! Filter unique vertices to share the nodes of the mesh
    //! ------------------------------------------------------
    BRepBuilderAPI_CellFilter uniqueVertices(Precision::Confusion());
    BRepBuilderAPI_VertexInspector inspector(Precision::Confusion());

    std::vector<BRepBuilderAPI_CellFilter> vecCellFilter;
    std::vector<BRepBuilderAPI_VertexInspector> vecVertexInspector;

    //! ----------------------------------------------------------
    //! identical technique for each face, including the fake "0"
    //! ----------------------------------------------------------
    for(int i=0; i<vecFaceStlMesh_Mesh.size();i++)
    {
        BRepBuilderAPI_CellFilter anUniqueVertices(Precision::Confusion());
        BRepBuilderAPI_VertexInspector anInspector(Precision::Confusion());
        vecCellFilter.push_back(anUniqueVertices);
        vecVertexInspector.push_back(anInspector);
    }

    //! -------------
    //! main reading
    //! -------------
    for(int iTri = 0; iTri < nbTris; iTri++)
    {
        char x[256]="", y[256]="", z[256]="";

        //! ----------------------------------
        //! reading the facet normal
        //! error should be properly reported
        //! ----------------------------------
        if (3 != fscanf(file,"%*s %*s %80s %80s %80s\n", x, y, z)) break;
        gp_XYZ aN (Atof(x), Atof(y), Atof(z));

        //! -------------------------------
        //! skip the keywords "outer loop"
        //! -------------------------------
        if (fscanf(file,"%*s %*s") < 0) break;

        //! -----------------------------------
        //! reading vertex
        //! errors should be properly reported
        //! -----------------------------------
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z)) break;
        gp_XYZ aV1 (Atof(x), Atof(y), Atof(z));
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z)) break;
        gp_XYZ aV2 (Atof(x), Atof(y), Atof(z));
        if (3 != fscanf(file,"%*s %80s %80s %80s\n", x, y, z)) break;
        gp_XYZ aV3 (Atof(x), Atof(y), Atof(z));

        //! ---------------------------------------------------------------
        //! here the facet must be built and put in the mesh datastructure
        //! ---------------------------------------------------------------
        i1 = AddVertex(ReadMesh, uniqueVertices, inspector, aV1);
        i2 = AddVertex(ReadMesh, uniqueVertices, inspector, aV2);
        i3 = AddVertex(ReadMesh, uniqueVertices, inspector, aV3);

        //! global element number: number of the element within the overall surface mesh
        int globalElementNumber = ReadMesh->AddTriangle (i1, i2, i3, aN.X(), aN.Y(), aN.Z());

        //! ----------------------------
        //! skip the keywords "endloop"
        //! ----------------------------
        if (fscanf(file,"%*s") < 0) break;

        //! --------------------------------------------
        //! read the tag
        //! ii1, ii2, ii3 are local (face) node numbers
        //! --------------------------------------------
        unsigned int tag;
        if(fscanf(file,"%u",&tag)==1)
        {
            //cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->____tag (face number): "<<tag<<"____"<<endl;

            int ii1 = AddVertex(vecFaceStlMesh_Mesh[tag], vecCellFilter.at(tag), vecVertexInspector.at(tag), aV1);
            int ii2 = AddVertex(vecFaceStlMesh_Mesh[tag], vecCellFilter.at(tag), vecVertexInspector.at(tag), aV2);
            int ii3 = AddVertex(vecFaceStlMesh_Mesh[tag], vecCellFilter.at(tag), vecVertexInspector.at(tag), aV3);

            //cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->local ____("<<ii1<<", "<<ii2<<", "<<ii3<<")____"<<endl;
            //cout<<"ExtendedRWStl::ReadExtendedSTLAscii()->global____("<<i1<<", "<<i2<<", "<<i3<<")____"<<endl;

            const occHandle(StlMesh_Mesh) &faceSTLMesh = vecFaceStlMesh_Mesh[tag];
            int localElementNumber = faceSTLMesh->AddTriangle(ii1, ii2, ii3, aN.X(), aN.Y(), aN.Z());

            const QMap<int,int> &localToGlobal_nodeID = maps_localToGlobal_nodeIDs.value(tag);
            if(!localToGlobal_nodeID.contains(ii1)) maps_localToGlobal_nodeIDs.operator [](tag).insert(ii1,i1);
            if(!localToGlobal_nodeID.contains(ii2)) maps_localToGlobal_nodeIDs.operator [](tag).insert(ii2,i2);
            if(!localToGlobal_nodeID.contains(ii3)) maps_localToGlobal_nodeIDs.operator [](tag).insert(ii3,i3);

            maps_localToGlobal_elementIDs.operator [](tag).insert(localElementNumber,globalElementNumber);
        }

        //! -----------------------------
        //! skip the keywords "endfacet"
        //! -----------------------------
        if (fscanf(file,"%*s") < 0) break;
    }
    return ReadMesh;
}
