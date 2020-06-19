//! custom includes
#include "topologytools.h"
#include "ccout.h"

//! OCC
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Builder.hxx>
#include <TopExp_Explorer.hxx>

TopologyTools::TopologyTools()
{
    ;
}

//! -------------------------------
//! function: generateLocationPair
//! details:
//! -------------------------------
QVector<GeometryTag> TopologyTools::generateLocationPairs(geometryDataBase *gDB, const ListOfShape& scope)
{
    //cout<<"TopologyTools::generateLocationPairs()->____function called on "<<scope.Extent()<<" shapes____"<<endl;
    QVector<GeometryTag> vecLocs;
    for(TopTools_ListIteratorOfListOfShape it(scope);it.More();it.Next())
    {
        const TopoDS_Shape &curShape = it.Value();
        curShape.ShapeType();

        GeometryTag loc;
        if(it.Value().ShapeType()==TopAbs_SOLID)
        {
            int bodyIndex = gDB->bodyMap.key(it.Value());
            loc.isParent=true;
            loc.parentShapeNr=bodyIndex;
            loc.subTopNr=bodyIndex;
            loc.subShapeType=TopAbs_SOLID;
        }
        else
        {
            TopAbs_ShapeEnum type;
            int subshapeNumber;
            int parentShapeNumber;
            gDB->getSubShapeNr(it.Value(),parentShapeNumber,subshapeNumber,type);
            loc.isParent=false;
            loc.parentShapeNr=parentShapeNumber;
            loc.subTopNr=subshapeNumber;
            loc.subShapeType=type;
        }
        vecLocs.push_back(loc);
    }
    //cout<<"TopologyTools::generateLocationPairs()->____exiting function____"<<endl;
    return vecLocs;
}

//! ---------------------------
//! function: explore compound
//! details:
//! ---------------------------
TopoDS_Shape TopologyTools::exploreCompound(const TopoDS_Compound &aComp,
                                            QList<TopoDS_Shape> &csolids,
                                            QList<TopoDS_Shape> &solids,
                                            QList<TopoDS_Shape> &shells,
                                            QList<TopoDS_Shape> &faces,
                                            QList<TopoDS_Shape> &wires,
                                            QList<TopoDS_Shape> &edges,
                                            bool load3DBodies,
                                            bool loadSurfaceBodies,
                                            bool loadLineBodies)
{
    int k= 0;
    for(TopoDS_Iterator anIt(aComp); anIt.More(); anIt.Next())
    {
        k++;
        //! ----------------------------
        //! handle the composite solids
        //! ----------------------------
        TopoDS_Shape aShape = anIt.Value();
        if(aShape.ShapeType() == TopAbs_COMPSOLID)
        {
            cout<<"Shape nr: "<<k<<"____found a COMPSOLID____"<<endl;
            csolids<<aShape;
        }

        //! ------------------
        //! handle the solids
        //! ------------------
        if(aShape.ShapeType() == TopAbs_SOLID)
        {
            cout<<"Shape nr: "<<k<<"____found a SOLID____"<<endl;
            solids<<aShape;
        }

        //! -----------------------------------------------------------
        //! the isolated face is inserted into the compound as a SHELL
        //! made of a single FACE
        //! -----------------------------------------------------------
        if(aShape.ShapeType() == TopAbs_SHELL)
        {
            int count = 0;
            for(TopExp_Explorer anExp(aShape,TopAbs_FACE); anExp.More(); anExp.Next())
            {
                count++;
            }
            if(count == 1)
            {
                cout<<"Shape nr: "<<k<<"____found a FACE____"<<endl;
                faces<<aShape;
            }
        }

        //! ---------------------------------------------------------------
        //! the isolated shell is inserted into the compound as a COMPOUND
        //! made of single faces
        //! ---------------------------------------------------------------
        if(aShape.ShapeType() == TopAbs_COMPOUND)
        {
            cout<<"Shape nr: "<<k<<"____found a COMPOUND____"<<endl;
            int edgeCount = 0;
            TopoDS_Builder wireBuilder;
            TopoDS_Wire wire;
            wireBuilder.MakeWire(wire);

            int faceCount = 0;
            TopoDS_Builder shellBuilder;
            TopoDS_Shell shell;
            shellBuilder.MakeShell(shell);
            for(TopoDS_Iterator it(aShape); it.More(); it.Next())
            {
                TopoDS_Shape aShape1 = it.Value();
                if(aShape1.ShapeType() == TopAbs_SHELL)
                {
                    faceCount++;
                    for(TopExp_Explorer anExp(aShape1,TopAbs_FACE); anExp.More(); anExp.Next())
                    {
                        TopoDS_Shape curShape = anExp.Current();
                        shellBuilder.Add(shell,curShape);
                    }
                }

                //! -------------
                //! experimental
                //! -------------
                if(aShape1.ShapeType() == TopAbs_EDGE)
                {
                    edgeCount++;
                    wireBuilder.Add(wire,aShape1);
                }
                //! -----------------
                //! end experimental
                //! -----------------
            }
            if(faceCount>1)
            {
                cout<<"Shape nr: "<<k<<"____ACTUALLY A SHELL____"<<endl;
                shells<<shell;
            }
            if(edgeCount>1)
            {
                cout<<"Shape nr: "<<k<<"____ACTUALLY A WIRE____"<<endl;
                wires<<wire;
            }
            if(edgeCount==1)
            {
                for(TopExp_Explorer anExp(wire,TopAbs_EDGE); anExp.More(); anExp.Next())
                {
                    //! this cycle has one iteration, by definition
                    edges<<anExp.Current();
                    cout<<"Shape nr: "<<k<<"____ACTUALLY AN EDGE____"<<endl;
                }
            }
        }
    }

    cout<<"\n@____Topology tools summary____@"<<endl;
    cout<<"____Nb solids: "<<solids.length()<<"____"<<endl;
    cout<<"____Nb shells: "<<shells.length()<<"____"<<endl;
    cout<<"____Nb faces: "<<faces.length()<<"____"<<endl;
    cout<<"____Nb wires: "<<wires.length()<<"____"<<endl;
    cout<<"____Nb edges: "<<edges.length()<<"____"<<endl;
    cout<<"@____End of the summary________@\n"<<endl;

    TopoDS_Compound compound;
    TopoDS_Builder compoundBuilder;
    compoundBuilder.MakeCompound(compound);

    if(load3DBodies==true)
    {
        for(int i=0; i<csolids.length(); i++) compoundBuilder.Add(compound,csolids.at(i));
        for(int i=0; i<solids.length(); i++) compoundBuilder.Add(compound,solids.at(i));
    }
    if(loadSurfaceBodies==true)
    {
        for(int i=0; i<shells.length(); i++) compoundBuilder.Add(compound,shells.at(i));
        for(int i=0; i<faces.length(); i++) compoundBuilder.Add(compound,faces.at(i));
    }
    if(loadLineBodies==true)
    {
        for(int i=0; i<wires.length(); i++) compoundBuilder.Add(compound,wires.at(i));
        for(int i=0; i<edges.length(); i++) compoundBuilder.Add(compound,edges.at(i));
    }
    return compound;
}


//! -----------------------------------------------------
//! function: makeShellFromoFaces
//! details:  a container for a specialized BRep_Builder
//! -----------------------------------------------------
#include <TopoDS_Builder.hxx>
TopoDS_Shell TopologyTools::makeShellFromFaces(const QList<TopoDS_Face> &aFaceList)
{
    //! -------------------------------
    //! build a face compound and mesh
    //! -------------------------------
    TopoDS_Builder builder;
    TopoDS_Shell aShell;
    builder.MakeShell(aShell);
    for(int i=0; i<aFaceList.length(); i++)
    {
        const TopoDS_Face &F = aFaceList.at(i);
        builder.Add(aShell,F);
    }
    return aShell;
}

