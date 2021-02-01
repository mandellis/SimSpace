#include "geometrydoctor.h"

#include <TopExp_Explorer.hxx>
#include <TopoDS_CompSolid.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>

#include <TopoDS.hxx>

#include <ShapeBuild_ReShape.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <ShapeFix_FixSmallFace.hxx>

#include <BRep_Tool.hxx>
#include <BRepLib.hxx>
#include <BRepOffsetAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepGProp.hxx>

#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

#include <GProp_GProps.hxx>
#include <TopExp.hxx>

geometryDoctor::geometryDoctor(TopoDS_Shape shape)
{
    myShape=shape;
    BuildFMap();
}

// function: Repair
TopoDS_Shape geometryDoctor::Repair(double tolerance,Standard_Boolean fixsmalledges,
                                    Standard_Boolean fixspotstripfaces,
                                    Standard_Boolean sewfaces,
                                    Standard_Boolean makesolids)
{
    int nrc = 0, nrcs = 0,
            nrso = somap.Extent(),
            nrsh = shmap.Extent(),
            nrf = fmap.Extent(),
            nrw = wmap.Extent(),
            nre = emap.Extent(),
            nrv = vmap.Extent();

    cout<<"INITIAL NUMBER OF SOLIDS "<<nrso<<endl;
    cout<<"INITIAL NUMBER OF SHELLS "<<nrsh<<endl;
    cout<<"INITIAL NUMBER OF FACES "<<nrf<<endl;
    cout<<"INITIAL NUMBER OF WIRES "<<nrw<<endl;
    cout<<"INITIAL NUMBER OF EDGES "<<nre<<endl;
    cout<<"INITIAL NUMBER OF VERTEX "<<nrv<<endl;

    TopExp_Explorer exp0;
    TopExp_Explorer exp1;

    for (exp0.Init(myShape,TopAbs_COMPOUND); exp0.More(); exp0.Next())nrc++;
    for (exp0.Init(myShape,TopAbs_COMPSOLID); exp0.More(); exp0.Next())nrcs++;

    cout<<"INITIAL NUMBER OF COMPUNDS "<<nrc<<endl;
    cout<<"INITIAL NUMBER OF COMPSOLIDS "<<nrcs<<endl;

    double surfacecont = 0;

    // MI SFUGGE IL SIGNIFICATO DI QUESTA PARTE DI CODICE
    {
        Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
        rebuild->Apply(myShape);
        for (exp1.Init (myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
        {
            TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
            if (BRep_Tool::Degenerated(edge))rebuild->Remove(edge, false);
       }
       myShape = rebuild->Apply(myShape);
    }
    BuildFMap();

    // CALCOLO AREA TOTALE
    for (exp0.Init (myShape, TopAbs_FACE); exp0.More(); exp0.Next())
    {
        TopoDS_Face face = TopoDS::Face(exp0.Current());
        GProp_GProps system;
        BRepGProp::SurfaceProperties(face, system);
        surfacecont += system.Mass();
    }


    cout << "STARTING HEALING USING TOLERANCE: " << tolerance << endl;
    {
        cout << endl << "Repairing faces" << endl;

        Handle(ShapeFix_Face) sff;
        Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
        rebuild->Apply(myShape);

        for (exp0.Init (myShape, TopAbs_FACE); exp0.More(); exp0.Next())
        {
            TopoDS_Face face = TopoDS::Face (exp0.Current());
            sff = new ShapeFix_Face (face);
            //GB
            sff->FixAddNaturalBoundMode() = Standard_False;
            sff->FixSmallAreaWireMode() = Standard_False;
            //GB
            sff->Perform();

            if(sff->Status(ShapeExtend_DONE1) ||
                    sff->Status(ShapeExtend_DONE2) ||
                    sff->Status(ShapeExtend_DONE3) ||
                    sff->Status(ShapeExtend_DONE4) ||
                    sff->Status(ShapeExtend_DONE5))
                    {
                        cout << "repaired face " << fmap.FindIndex(face) << " ";
                        if(sff->Status(ShapeExtend_DONE1))
                            cout << "(some wires are fixed)" <<endl;
                        else if(sff->Status(ShapeExtend_DONE2))
                            cout << "(orientation of wires fixed)" <<endl;
                        else if(sff->Status(ShapeExtend_DONE3))
                            cout << "(missing seam added)" <<endl;
                        else if(sff->Status(ShapeExtend_DONE4))
                            cout << "(small area wire removed)" <<endl;
                        else if(sff->Status(ShapeExtend_DONE5))
                            cout << "(natural bounds added)" <<endl;
                        TopoDS_Face newface = sff->Face();
                        rebuild->Replace(face, newface, Standard_False);
                    }
       }
       myShape = rebuild->Apply(myShape);
    }

    // MI SFUGGE IL SIGNIFICATO DI QUESTA PARTE DI CODICE
    {
       Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
       rebuild->Apply(myShape);
       for (exp1.Init (myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
       {
           TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
           if (BRep_Tool::Degenerated(edge))rebuild->Remove(edge, false);
       }
       myShape = rebuild->Apply(myShape);
    }

    // FIX SMALL EDGES
    if (fixsmalledges)
    {
       cout << endl << "Fixing small edges" << endl;

       Handle(ShapeFix_Wire) sfw;
       Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
       rebuild->Apply(myShape);

       for (exp0.Init (myShape, TopAbs_FACE); exp0.More(); exp0.Next())
       {
          TopoDS_Face face = TopoDS::Face(exp0.Current());

          for (exp1.Init (face, TopAbs_WIRE); exp1.More(); exp1.Next())
          {
             TopoDS_Wire oldwire = TopoDS::Wire(exp1.Current());
             sfw = new ShapeFix_Wire (oldwire,face,tolerance);
             sfw->ModifyTopologyMode() = Standard_True;
             sfw->ClosedWireMode() = Standard_True;

             bool replace = false;

             replace = sfw->FixReorder() || replace;
             replace = sfw->FixConnected() || replace;

             if (sfw->FixSmall (Standard_False, tolerance) && !(sfw->StatusSmall(ShapeExtend_FAIL1) ||
                sfw->StatusSmall(ShapeExtend_FAIL2) ||
                sfw->StatusSmall(ShapeExtend_FAIL3)))
             {
                 cout << "Fixed small edge in wire " << wmap.FindIndex (oldwire) << endl;
                 replace = true;
             }
             else if (sfw->StatusSmall(ShapeExtend_FAIL1))
                 cerr << "Failed to fix small edge in wire " << wmap.FindIndex(oldwire)
                      << ", edge cannot be checked (no 3d curve and no pcurve)" << endl;
             else if (sfw->StatusSmall(ShapeExtend_FAIL2))
                 cerr << "Failed to fix small edge in wire " << wmap.FindIndex (oldwire)
                       << ", edge is null-length and has different vertives at begin and end, and lockvtx is True or ModifiyTopologyMode is False" << endl;
             else if (sfw->StatusSmall(ShapeExtend_FAIL3))
                 cerr << "Failed to fix small edge in wire " << wmap.FindIndex (oldwire)
                       << ", CheckConnected has failed" << endl;

             replace = sfw->FixEdgeCurves() || replace;
             replace = sfw->FixDegenerated() || replace;
             replace = sfw->FixSelfIntersection() || replace;
             replace = sfw->FixLacking(Standard_True) || replace;

             if(replace)
             {
                 TopoDS_Wire newwire = sfw->Wire();
                 rebuild->Replace(oldwire, newwire, Standard_False);
             }
             //delete sfw; sfw = NULL;
          }
       }
       myShape = rebuild->Apply(myShape);

       {
          BuildFMap();
          Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
          rebuild->Apply(myShape);

          for (exp1.Init (myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
          {
             TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
             if (vmap.FindIndex(TopExp::FirstVertex (edge)) ==
                vmap.FindIndex(TopExp::LastVertex (edge)))
             {
                GProp_GProps system;
                BRepGProp::LinearProperties(edge, system);
                if (system.Mass() < tolerance)
                {
                   cout << "removing degenerated edge " << emap.FindIndex(edge)
                      << " from vertex " << vmap.FindIndex(TopExp::FirstVertex (edge))
                      << " to vertex " << vmap.FindIndex(TopExp::LastVertex (edge)) << endl;
                   rebuild->Remove(edge, false);
                }
             }
          }
          myShape = rebuild->Apply(myShape);
          //delete rebuild; rebuild = NULL;
       }

       {
          Handle_ShapeBuild_ReShape rebuild = new ShapeBuild_ReShape;
          rebuild->Apply(myShape);
          for (exp1.Init (myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
          {
             TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
             if ( BRep_Tool::Degenerated(edge) )
                rebuild->Remove(edge, false);
          }
          myShape = rebuild->Apply(myShape);
       }

       Handle(ShapeFix_Wireframe) sfwf = new ShapeFix_Wireframe;
       sfwf->SetPrecision(tolerance);
       sfwf->Load(myShape);
       sfwf->ModeDropSmallEdges() = Standard_True;

       // MODIFICA GB
       Bnd_Box bb;
       BRepBndLib::Add (myShape, bb);
       double x1,y1,z1,x2,y2,z2,bbDiam;
       bb.Get (x1,y1,z1,x2,y2,z2);
       bbDiam=sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
       sfwf->SetPrecision(bbDiam);
       //

       if (sfwf->FixWireGaps())
       {
          cout << endl << "Fixing wire gaps" << endl;
          if (sfwf->StatusWireGaps(ShapeExtend_OK)) cout << "no gaps found" << endl;
          if (sfwf->StatusWireGaps(ShapeExtend_DONE1)) cout << "some 2D gaps fixed" << endl;
          if (sfwf->StatusWireGaps(ShapeExtend_DONE2)) cout << "some 3D gaps fixed" << endl;
          if (sfwf->StatusWireGaps(ShapeExtend_FAIL1)) cout << "failed to fix some 2D gaps" << endl;
          if (sfwf->StatusWireGaps(ShapeExtend_FAIL2)) cout << "failed to fix some 3D gaps" << endl;
       }
       sfwf->SetPrecision(tolerance);

       // QUESTA PORZIONE DI CODICE E' PER IL CONTROLLO "ESTERNO" DELL'ESECUZIONE
       {
          for (exp1.Init (myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
          {
             TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
             if ( BRep_Tool::Degenerated(edge) )
                cout << "degenerated edge at position 4" << endl;
          }
       }

       if (sfwf->FixSmallEdges())
       {
          cout << endl << "Fixing wire frames" << endl;
          if (sfwf->StatusSmallEdges(ShapeExtend_OK)) cout << "no small edges found" << endl;
          if (sfwf->StatusSmallEdges(ShapeExtend_DONE1)) cout << "some small edges fixed" << endl;
          if (sfwf->StatusSmallEdges(ShapeExtend_FAIL1)) cout << "failed to fix some small edges" << endl;
       }

       myShape = sfwf->Shape();

       //delete sfwf; sfwf = NULL;
       //delete rebuild; rebuild = NULL;
    }

    // QUESTA PORZIONE DI CODICE E' PER IL CONTROLLO "ESTERNO" DELL'ESECUZIONE
    {
       for (exp1.Init (myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
       {
          TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
          if ( BRep_Tool::Degenerated(edge) )
             cout << "degenerated edge at position 5" << endl;
       }
    }

    if (fixspotstripfaces)
    {
       cout << endl << "Fixing spot and strip faces" << endl;
       Handle(ShapeFix_FixSmallFace) sffsm = new ShapeFix_FixSmallFace();
       sffsm -> Init (myShape);
       sffsm -> SetPrecision (tolerance);
       sffsm -> Perform();

       myShape = sffsm -> FixShape();
       //delete sffsm; sffsm = NULL;
    }

    // QUESTA PORZIONE DI CODICE E' PER IL CONTROLLO "ESTERNO" DELL'ESECUZIONE
    {
       for (exp1.Init(myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
       {
          TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
          if ( BRep_Tool::Degenerated(edge) )
             cout << "degenerated edge at position 6" << endl;
       }
    }

    if (sewfaces)
    {
       cout << endl << "Sewing faces" << endl;

       BRepOffsetAPI_Sewing sewedObj(tolerance);

       for (exp0.Init (myShape, TopAbs_FACE); exp0.More(); exp0.Next())
       {
          TopoDS_Face face = TopoDS::Face (exp0.Current());
          sewedObj.Add (face);
       }

       sewedObj.Perform();

       if (!sewedObj.SewedShape().IsNull())
          myShape = sewedObj.SewedShape();
       else
          cout << " not possible";
    }

    // ANCHE QUESTA MI SFUGGE. PROBABILMENTE SOLO PER CONTROLLO
    {
       Handle_ShapeBuild_ReShape rebuild = new ShapeBuild_ReShape;
       rebuild->Apply(myShape);
       for (exp1.Init (myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
       {
          TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
          if (BRep_Tool::Degenerated(edge))
             rebuild->Remove(edge, false);
       }
       myShape = rebuild->Apply(myShape);
    }

    if (makesolids)
    {
       cout << endl << "Making solids" << endl;

       BRepBuilderAPI_MakeSolid ms;
       int count = 0;
       for (exp0.Init(myShape, TopAbs_SHELL); exp0.More(); exp0.Next())
       {
          count++;
          ms.Add (TopoDS::Shell(exp0.Current()));
       }

       if (!count)
       {
          cout << " not possible (no shells)" << endl;
       }
       else
       {
          BRepCheck_Analyzer ba(ms);
          if (ba.IsValid ())
          {
             Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape;
             sfs->Init (ms);
             sfs->SetPrecision(tolerance);
             sfs->SetMaxTolerance(tolerance);
             sfs->Perform();
             myShape = sfs->Shape();

             for (exp0.Init(myShape, TopAbs_SOLID); exp0.More(); exp0.Next())
             {
                TopoDS_Solid solid = TopoDS::Solid(exp0.Current());
                TopoDS_Solid newsolid = solid;
                BRepLib::OrientClosedSolid (newsolid);
                Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
                //rebuild->Apply(myShape);
                rebuild->Replace(solid,newsolid,Standard_False);
                TopoDS_Shape newshape = rebuild->Apply(myShape, TopAbs_COMPSOLID);//, 1);
                //TopoDS_Shape newshape = rebuild->Apply(shape);
                myShape = newshape;
             }
             //delete sfs; sfs = NULL;
          }
          else
             cout << " not possible" << endl;
       }
    }

    // QUESTA PORZIONE DI CODICE E' PER IL CONTROLLO "ESTERNO" DELL'ESECUZIONE
    {
       for (exp1.Init (myShape, TopAbs_EDGE); exp1.More(); exp1.Next())
       {
          TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
          if (BRep_Tool::Degenerated(edge))
             cout << "degenerated edge at position 8" << endl;
       }
    }

    double newsurfacecont = 0;

    for (exp0.Init (myShape, TopAbs_FACE); exp0.More(); exp0.Next())
    {
       TopoDS_Face face = TopoDS::Face(exp0.Current());
       GProp_GProps system;
       BRepGProp::SurfaceProperties(face, system);
       newsurfacecont += system.Mass();
    }

    int nnrc = 0, nnrcs = 0,
       nnrso = somap.Extent(),
       nnrsh = shmap.Extent(),
       nnrf = fmap.Extent(),
       nnrw = wmap.Extent(),
       nnre = emap.Extent(),
       nnrv = vmap.Extent();

    for (exp0.Init(myShape, TopAbs_COMPOUND); exp0.More(); exp0.Next()) nnrc++;
    for (exp0.Init(myShape, TopAbs_COMPSOLID); exp0.More(); exp0.Next()) nnrcs++;

    cout << "Compounds       : " << nnrc << " (" << nrc << ")" << endl;
    cout << "Composite solids: " << nnrcs << " (" << nrcs << ")" << endl;
    cout << "Solids          : " << nnrso << " (" << nrso << ")" << endl;
    cout << "Shells          : " << nnrsh << " (" << nrsh << ")" << endl;
    cout << "Wires           : " << nnrw << " (" << nrw << ")" << endl;
    cout << "Faces           : " << nnrf << " (" << nrf << ")" << endl;
    cout << "Edges           : " << nnre << " (" << nre << ")" << endl;
    cout << "Vertices        : " << nnrv << " (" << nrv << ")" << endl;
    cout << endl;
    cout << "Totol surface area : " << newsurfacecont << " (" << surfacecont << ")" << endl;
    cout << endl;

    return(myShape);
}

// function: BuildFMap
void geometryDoctor::BuildFMap()
{
    shmap.Clear();
    somap.Clear();
    fmap.Clear();
    wmap.Clear();
    emap.Clear();
    vmap.Clear();

    TopExp_Explorer exp0, exp1, exp2, exp3, exp4, exp5;

    for (exp0.Init(myShape,TopAbs_COMPOUND);exp0.More();exp0.Next())
    {
       TopoDS_Compound compound = TopoDS::Compound (exp0.Current());
       int i = 0;
       for(exp1.Init(compound, TopAbs_SHELL);exp1.More(); exp1.Next())
       {
           ++i;
       }
    }

    for (exp0.Init(myShape, TopAbs_SOLID);exp0.More();exp0.Next())
    {
       TopoDS_Solid solid = TopoDS::Solid (exp0.Current());
       if (somap.FindIndex(solid) < 1)
       {
          somap.Add (solid);
          for (exp1.Init(solid, TopAbs_SHELL);
            exp1.More(); exp1.Next())
            {
             TopoDS_Shell shell = TopoDS::Shell (exp1.Current());
             if (shmap.FindIndex(shell) < 1)
             {
                shmap.Add (shell);
                for (exp2.Init(shell, TopAbs_FACE);exp2.More();exp2.Next())
                {
                   TopoDS_Face face = TopoDS::Face(exp2.Current());
                   if (fmap.FindIndex(face) < 1)
                   {
                      fmap.Add (face);
                      for (exp3.Init(exp2.Current(), TopAbs_WIRE);exp3.More(); exp3.Next())
                      {
                         TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
                         if (wmap.FindIndex(wire) < 1)
                         {
                            wmap.Add (wire);
                            for (exp4.Init(exp3.Current(), TopAbs_EDGE);exp4.More(); exp4.Next())
                            {
                               TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
                               if (emap.FindIndex(edge) < 1)
                               {
                                  emap.Add (edge);
                                  for (exp5.Init(exp4.Current(), TopAbs_VERTEX);exp5.More(); exp5.Next())
                                  {
                                     TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                                     if (vmap.FindIndex(vertex) < 1)
                                        vmap.Add (vertex);
                                  }
                               }
                            }
                         }
                      }
                   }
                }
             }
          }
       }
    }

    // Free Shells
    for (exp1.Init(myShape, TopAbs_SHELL, TopAbs_SOLID); exp1.More(); exp1.Next())
    {
       TopoDS_Shell shell = TopoDS::Shell(exp1.Current());
       if (shmap.FindIndex(shell) < 1)
       {
          shmap.Add (shell);

          for (exp2.Init(shell, TopAbs_FACE); exp2.More(); exp2.Next())
          {
             TopoDS_Face face = TopoDS::Face(exp2.Current());
             if (fmap.FindIndex(face) < 1)
             {
                fmap.Add (face);

                for (exp3.Init(face, TopAbs_WIRE); exp3.More(); exp3.Next())
                {
                   TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
                   if (wmap.FindIndex(wire) < 1)
                   {
                      wmap.Add (wire);

                      for (exp4.Init(wire, TopAbs_EDGE); exp4.More(); exp4.Next())
                      {
                         TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
                         if (emap.FindIndex(edge) < 1)
                         {
                            emap.Add (edge);
                            for (exp5.Init(edge, TopAbs_VERTEX); exp5.More(); exp5.Next())
                            {
                               TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                               if (vmap.FindIndex(vertex) < 1)
                                  vmap.Add (vertex);
                            }
                         }
                      }
                   }
                }
             }
          }
       }
    }

    // Free Faces
    for (exp2.Init(myShape, TopAbs_FACE, TopAbs_SHELL); exp2.More(); exp2.Next())
    {
       TopoDS_Face face = TopoDS::Face(exp2.Current());
       if (fmap.FindIndex(face) < 1)
       {
          fmap.Add (face);

          for (exp3.Init(exp2.Current(), TopAbs_WIRE); exp3.More(); exp3.Next())
          {
             TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
             if (wmap.FindIndex(wire) < 1)
             {
                wmap.Add (wire);

                for (exp4.Init(exp3.Current(), TopAbs_EDGE); exp4.More(); exp4.Next())
                {
                   TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
                   if (emap.FindIndex(edge) < 1)
                   {
                      emap.Add (edge);
                      for (exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next())
                      {
                         TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                         if (vmap.FindIndex(vertex) < 1)
                            vmap.Add (vertex);
                      }
                   }
                }
             }
          }
       }
    }

    // Free Wires
    for (exp3.Init(myShape, TopAbs_WIRE, TopAbs_FACE); exp3.More(); exp3.Next())
    {
       TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
       if (wmap.FindIndex(wire) < 1)
       {
          wmap.Add (wire);

          for (exp4.Init(exp3.Current(), TopAbs_EDGE); exp4.More(); exp4.Next())
          {
             TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
             if (emap.FindIndex(edge) < 1)
             {
                emap.Add (edge);
                for (exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next())
                {
                   TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                   if (vmap.FindIndex(vertex) < 1)
                      vmap.Add (vertex);
                }
             }
          }
       }
    }

    // Free Edges
    for (exp4.Init(myShape, TopAbs_EDGE, TopAbs_WIRE); exp4.More(); exp4.Next())
    {
       TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
       if (emap.FindIndex(edge) < 1)
       {
          emap.Add (edge);
          for (exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next())
          {
             TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
             if (vmap.FindIndex(vertex) < 1)
                vmap.Add (vertex);
          }
       }
    }

    // Free Vertices
    for (exp5.Init(myShape, TopAbs_VERTEX, TopAbs_EDGE); exp5.More(); exp5.Next())
    {
       TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
       if (vmap.FindIndex(vertex) < 1)
          vmap.Add (vertex);
    }
}

// function: Analyze
#include <ShapeAnalysis_Shell.hxx>
void geometryDoctor::Analyze()
{
    ShapeAnalysis_Shell sas;
    for(int i=1;i<=shmap.Extent();i++)
    {
        const TopoDS_Shell &shell=TopoDS::Shell(shmap.FindKey(i));
        sas.LoadShells(shell);
        sas.CheckOrientedShells(shell,Standard_True);
        if(sas.HasBadEdges())
            cout<<"Shell is invalid: "<<endl;
        if(sas.HasFreeEdges())
            cout<<"Shell is open: "<<endl;
    }
}
