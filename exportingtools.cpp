//! ----------------
//! custom includes
//! ----------------
#include "exportingtools.h"
#include "ng_meshvs_datasource1d.h"
#include "ng_meshvs_datasourceface.h"
#include "meshtools.h"
#include "simulationnodeclass.h"
#include "postobject.h"
#include "geometrytag.h"

//! ---
//! Qt
//! ---
#include <QFileDialog>
#include <QProgressDialog>
#include <QLabel>
#include <QMessageBox>

//! ----
//! OCC
//! ----
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Shape.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TColStd_HArray2OfReal.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Array1OfTriangle.hxx>
#include <NCollection_Array1.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <MeshVS_Mesh.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <MeshVS_EntityType.hxx>

//! ----
//! C++
//! ----
using namespace std;
#include <fstream>

exportingTools::exportingTools()
{
    ;
}

//! ----------------------
//! function: exportCloud
//! details:
//! ----------------------
void exportingTools::exportSTL(const occHandle(AIS_InteractiveContext) &theCTX, meshDataBase *theDB)
{
    //!cout<<"occPreGLWidget::exportSTL()->____function called____"<<endl;
    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(0,"Save as ",USERDATA_DIR,".stl", &selectedFilter,0);

    if(!fileName.isEmpty())
    {
        fileName.append(selectedFilter);

        //! the content of the selection
        TopAbs_ShapeEnum shapeTypeContent;
        std::vector<GeometryTag> vecLoc;
        std::vector<std::pair<int,int>> vecPairs;

        //! --------------------------------------
        //! retrieve the faces from the selection
        //! --------------------------------------
        theCTX->InitSelected();
        if(theCTX->MoreSelected())
        {
            shapeTypeContent=theCTX->SelectedShape().ShapeType();
            for(theCTX->InitSelected();theCTX->MoreSelected();theCTX->NextSelected())
            {
                const TopoDS_Shape &theShape = theCTX->SelectedShape();
                //! topology must be face for exporting an STL file
                switch(shapeTypeContent)
                {
                case TopAbs_FACE:
                {
                    std::pair<int,int> aPair;
                    int parentShapeIndex, childShapeIndex;
                    TopAbs_ShapeEnum type;
                    theDB->getSubShapeNr(theShape,parentShapeIndex,childShapeIndex,type);
                    aPair.first=parentShapeIndex;
                    aPair.second=childShapeIndex;
                    vecPairs.push_back(aPair);
                }
                    break;

                case TopAbs_SOLID:
                {
                    cout<<"function called"<<endl;
                    int bodyIndex = theDB->bodyMap.key(theShape);
                    TopTools_IndexedMapOfShape faceMap = theDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap;
                    for(int faceNr = 1; faceNr<=faceMap.Extent(); faceNr++)
                    {
                        std::pair<int,int> aPair;
                        int parentShapeIndex, childShapeIndex;
                        TopAbs_ShapeEnum type;
                        theDB->getSubShapeNr(faceMap.FindKey(faceNr),parentShapeIndex,childShapeIndex,type);
                        aPair.first=parentShapeIndex;
                        aPair.second=childShapeIndex;
                        vecPairs.push_back(aPair);
                    }
                }
                    break;
                }
            }
        }

        //! ---------------
        //! write the file
        //! ---------------
        //! faces whose stl mesh cannot be written
        std::list<pair<int,int>> invalidFaces;

        if(!vecPairs.empty())
        {
            int N = int(vecPairs.size());

            //! configure a progress dialog
            QProgressDialog progressDiag;
            progressDiag.setCancelButton(0);
            progressDiag.setRange(1,N);
            progressDiag.setWindowModality(Qt::WindowModal);
            progressDiag.setMinimumDuration(0);

            //! initialize the dialog
            int written =  0;
            progressDiag.setValue(written);
            QLabel progress_label;
            progress_label.setText("Saving STL file");
            progressDiag.setLabel(&progress_label);

            ofstream fileout;
            fileout.open(fileName.toStdString());
            int N1,N2,N3;

            char buf1[20], buf2[20], buf3[20];

            //! name ...
            fileout<<"solid"<<endl;

            //! iterate over faces
            for(std::vector<std::pair<int,int>>::iterator it=vecPairs.begin();it!=vecPairs.end();++it)
            {
                //! update the dialog
                written ++;

                std::pair<int,int> aPair = *it;
                //const TopoDS_Face &theFace = TopoDS::Face(theDB->MapOfBodyTopologyMap.alue(aPair.first).faceMap.FindKey(aPair.second));

                const occHandle(Ng_MeshVS_DataSourceFace) &theFaceMesh =
                        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(theDB->ArrayOfMeshDSOnFaces.getValue(aPair.first,aPair.second));

                occHandle(Poly_Triangulation) theTriangulation;
                bool isDone = MeshTools::toPolyTriangulation(theFaceMesh, theTriangulation);
                if(isDone)
                {
                    progressDiag.setLabelText(QString("Converting face nr. %1").arg(aPair.second));
                    progressDiag.setValue(written);

                    Poly_Array1OfTriangle trigs = theTriangulation->Triangles();
                    TColgp_Array1OfPnt nodes = theTriangulation->Nodes();

                    //! triangle data
                    double scale = 1.0;
                    for(int k=trigs.Lower(); k<=trigs.Upper(); k++)
                    {
                        trigs.Value(k).Get(N1,N2,N3);

                        double x_1 = nodes.Value(N1).X()*scale;
                        double y_1 = nodes.Value(N1).Y()*scale;
                        double z_1 = nodes.Value(N1).Z()*scale;
                        double x_2 = nodes.Value(N2).X()*scale;
                        double y_2 = nodes.Value(N2).Y()*scale;
                        double z_2 = nodes.Value(N2).Z()*scale;
                        double x_3 = nodes.Value(N3).X()*scale;
                        double y_3 = nodes.Value(N3).Y()*scale;
                        double z_3 = nodes.Value(N3).Z()*scale;

                        //! write normal
                        fileout<<"facet normal ";

                        double lx = x_2-x_1;
                        double ly = y_2-y_1;
                        double lz = y_2-y_1;
                        double mx = x_3-x_1;
                        double my = y_3-y_1;
                        double mz = y_3-y_1;
                        //! -----------
                        //! i   j   k
                        //! lx  ly  lz
                        //! mx  my  mz
                        //! -----------
                        double ax = ly*mz-lz*my;
                        double ay = lz*mx-lx*mz;
                        double az = lx*my-ly*mx;
                        double L = sqrt(ax*ax+ay*ay+az*az);
                        if(L>1e6)
                        {
                            ax = ax/L; ay = ay/L; az = az/L;
                        }
                        else
                        {
                            ax = ay = az = 0.0;
                        }
                        fileout<<ax<<" "<<ay<<" "<<az<<endl;

                        //std::vector<double> c = MeshTools::triangleGetNormal(gp_Pnt(x_1,y_1,z_1),gp_Pnt(x_2,y_2,z_2),gp_Pnt(x_3,y_3,z_3));
                        //fileout<<c[0]<<" "<<c[1]<<" "<<c[2]<<endl;

                        fileout<<"outer loop"<<endl;
                        sprintf(buf1,"%.6e",x_1);
                        sprintf(buf2,"%.6e",y_1);
                        sprintf(buf3,"%.6e",z_1);
                        fileout<<"vertex "<<buf1<<" "<<buf2<<" "<<buf3<<endl;

                        sprintf(buf1,"%.6e",x_2);
                        sprintf(buf2,"%.6e",y_2);
                        sprintf(buf3,"%.6e",z_2);
                        fileout<<"vertex "<<buf1<<" "<<buf2<<" "<<buf3<<endl;

                        sprintf(buf1,"%.6e",x_3);
                        sprintf(buf2,"%.6e",y_3);
                        sprintf(buf3,"%.6e",z_3);
                        fileout<<"vertex "<<buf1<<" "<<buf2<<" "<<buf3<<endl;
                        fileout<<"endloop"<<endl;
                        fileout<<"endfacet"<<endl;
                    }
                }
                else
                {
                    invalidFaces.push_back(aPair);
                }
            }
            fileout<<"endsolid"<<endl;

            //! close the file
            fileout.close();
        }

        //! ---------------------------------------------
        //! diagnostic messages
        //! ---------------------------------------------
        if(invalidFaces.size()>0)
        {
            QString msg("The following faces cannot be saved: \n");
            for(std::list<pair<int,int>>::iterator it = invalidFaces.begin(); it!=invalidFaces.end(); ++it)
            {
                pair<int,int> aPair;
                std::string bodyName = theDB->MapOfBodyNames.value(aPair.first).toStdString();
                QString submsg = QString("Body: ").append(QString::fromStdString(bodyName).append(QString(" Face nr: %1").arg(aPair.second)));
                msg.append(submsg);
            }
            QMessageBox::warning(0,APPNAME,msg);
        }
    }
}

//! ----------------------
//! function: exportCloud
//! details:
//! ----------------------
void exportingTools::exportCloud(const occHandle(AIS_InteractiveContext) &theCTX, meshDataBase *theDB)
{
    //cout<<"occPreGLWidget::exportCloud()->____function called____"<<endl;
    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(0,"Save as ",USERDATA_DIR,".txt", &selectedFilter,0);

    if(!fileName.isEmpty())
    {
        fileName.append(selectedFilter);
        cout<<fileName.toStdString()<<endl;

        std::vector<std::pair<int,int>> vecPairs;
        std::pair<int,int> aPair;

        //! the content of the selection
        TopAbs_ShapeEnum shapeTypeContent;
        theCTX->InitSelected();
        if(theCTX->MoreSelected())
        {
            shapeTypeContent=theCTX->SelectedShape().ShapeType();
            for(theCTX->InitSelected();theCTX->MoreSelected();theCTX->NextSelected())
            {
                const TopoDS_Shape &theChildShape = theCTX->SelectedShape();
                int parentShapeIndex, childShapeIndex;
                TopAbs_ShapeEnum type;
                theDB->getSubShapeNr(theChildShape,parentShapeIndex,childShapeIndex,type);
                aPair.first=parentShapeIndex;
                aPair.second=childShapeIndex;
                vecPairs.push_back(aPair);
            }
        }

        ofstream fileout;
        fileout.open(fileName.toStdString());

        vector<std::pair<int,int>>::iterator it = vecPairs.begin();
        for(it=vecPairs.begin(); it!=vecPairs.end(); ++it)
        {
            const std::pair<int,int> &curPair = *it;
            int parentShapeIndex = curPair.first;
            int childShapeIndex = curPair.second;

            TColStd_PackedMapOfInteger theMap;
            occHandle(TColStd_HArray2OfReal) coords;
            switch(shapeTypeContent)
            {
            case TopAbs_EDGE:
            {
                const occHandle(Ng_MeshVS_DataSource1D) &theMeshDS = occHandle(Ng_MeshVS_DataSource1D)::
                        DownCast(theDB->ArrayOfMeshDSOnEdges.getValue(parentShapeIndex,childShapeIndex));
                if(!theMeshDS.IsNull())
                {
                    fileout<<"#body nr_"<<parentShapeIndex<<"edge nr_"<<childShapeIndex<<endl;
                    theMap = theMeshDS->GetAllNodes();
                    coords = theMeshDS->getNodesCoords();
                }
            }
                break;

            case TopAbs_FACE:
            {
                const occHandle(Ng_MeshVS_DataSourceFace) &theMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::
                        DownCast(theDB->ArrayOfMeshDSOnFaces.getValue(parentShapeIndex,childShapeIndex));
                if(!theMeshDS.IsNull())
                {
                    fileout<<"#body nr_"<<parentShapeIndex<<" face nr_"<<childShapeIndex<<endl;
                    theMap = theMeshDS->GetAllNodes();
                    coords = theMeshDS->getNodesCoordinates();
                }
            }
                break;
            }

            TColStd_MapIteratorOfPackedMapOfInteger it;
            int i=1;
            for(it.Initialize(theMap);it.More();it.Next())
            {
                //int nodeID = it.Key();
                double x = coords->Value(i,1);
                double y = coords->Value(i,2);
                double z = coords->Value(i,3);
                //cout<<nodeID<<"\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
                QString line = QString("%1 %2 %3").arg(x).arg(y).arg(z);
                fileout<<line.toStdString()<<endl;
                i++;
            }
        }
        fileout.close();
    }
}

//! ----------------------
//! function: exportSTEP
//! details:
//! ----------------------
bool exportingTools::exportSTEP(const TopoDS_Compound &aComp)
{
    cout<<"exportingTools::exportSTEP()->____function called____"<<endl;
    if(aComp.IsNull()) return false;

    QString workingDir = tools::getWorkingDir();
    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(0,"Save STEP file as ",workingDir,".step", &selectedFilter,0);
    if(!fileName.isEmpty())
    {
        //! -------------------------------------
        //! Write resulting compound to the file
        //! -------------------------------------
        STEPControl_Writer aWriter;
        IFSelect_ReturnStatus aStat = aWriter.Transfer(aComp,STEPControl_AsIs);
        fileName = fileName + selectedFilter;
        aStat = aWriter.Write(fileName.toStdString().c_str());
        if(aStat!= IFSelect_RetDone) return true;
        else return false;
    }
    return false;
}

//! -----------------------
//! function: exportResult
//! details:
//! -----------------------
#include <MeshVS_DeformedDataSource.hxx>
void exportingTools::exportNodalResult(SimulationNodeClass *aNode, const std::string &fileName)
{
    if(aNode->isAnalysisResult()==false) return;

    //! -----------------------------------
    //! data and meshes within post object
    //! -----------------------------------
    sharedPostObject aPostObject = aNode->getPropertyValue<sharedPostObject>("Post object");
    const std::map<GeometryTag,std::vector<std::map<int,double>>> &data = aPostObject->getData();
    const std::map<GeometryTag,occHandle(MeshVS_Mesh)> &meshes = aPostObject->getColoredMeshes();

    switch(aNode->getType())
    {
    case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
    {
        //cout<<"____exporting displacements____"<<endl;
        std::ofstream fout;
        fout.open(fileName);
        fout<<"#x\ty\tz\ttot\tdx\tdy\tdz"<<endl;
        for(std::map<GeometryTag,std::vector<std::map<int,double>>>::const_iterator it = data.cbegin(); it!=data.cend(); it++)
        {
            const GeometryTag &aTag = it->first;
            const std::vector<std::map<int,double>> &aRes = it->second;   // contains components

            const occHandle(MeshVS_Mesh) &aMeshVS = meshes.at(aTag);

            const occHandle(MeshVS_DataSource) &aMeshDS= aMeshVS->GetDataSource();
            if(aMeshDS.IsNull()) exit(100);
            MeshVS_DeformedDataSource* aMeshDS_def = (MeshVS_DeformedDataSource*)(aMeshDS.get());
            const occHandle(MeshVS_DataSource) &aMeshDS_nonDef = aMeshDS_def->GetNonDeformedDataSource();

            TColStd_MapIteratorOfPackedMapOfInteger itnodes = aMeshDS_nonDef->GetAllNodes();
            for(int localNodeID = 1; localNodeID<=aMeshDS_nonDef->GetAllNodes().Extent(); localNodeID++, itnodes.Next())
            {
                int globalNodeID = itnodes.Key();
                int NbNodes = -1;
                double b[3];
                TColStd_Array1OfReal coords(*b,1,3);
                MeshVS_EntityType eType;
                if(!aMeshDS_nonDef->GetGeom(globalNodeID,false,coords,NbNodes,eType)) continue;

                //! -----------------------------------------------
                //! write the coordinates and the component values
                //! -----------------------------------------------
                fout<<coords(1)<<"\t"<<coords(2)<<"\t"<<coords(3)<<"\t";
                for(int component = 0; component<aRes.size()-1; component++)
                {
                    double aVal = aRes.at(component).at(globalNodeID);
                    fout<<aVal<<"\t";
                }
                double aVal = aRes.at(aRes.size()-1).at(globalNodeID);
                fout<<aVal<<endl;
            }
        }
        fout.close();
    }
        break;
    default:
    {
        ;
    }
        break;
    }
}
