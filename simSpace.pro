#-------------------------------------------------
#
# Project created by QtCreator 2017-05-30T17:09:38
#
#-------------------------------------------------

QT       += core opengl gui

##############
# QCustomPlot
##############
QT += printsupport

greaterThan(QT_MAJOR_VERSION, 4: QT += widgets)

TARGET = simSpace
TEMPLATE = app
CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS   \
            QCUSTOMPLOT_USE_LIBRARY

##############
# QCustomPlot
##############
#DEFINES += QCUSTOMPLOT_USE_OPENGL

SOURCES += src/mesh/ng_meshvs_datasource1d.cpp \
    src/mesh/ng_meshvs_datasource2d.cpp \
    src/mesh/ng_meshvs_datasource3d.cpp \
    src/mesh/ng_meshvs_datasourceface.cpp \
    src/mesh/tetgenmesher.cpp \
    src/mesh/ng_mesher2.cpp \
    src/viewer/occGLwidget.cpp \
    src/viewer/occPreGLwidget.cpp \
    src/mesh/mesherclass.cpp \
    src/mesh/edgemeshdatasourcebuildesclass.cpp \
    src/viewer/occpostwidget.cpp \
    src/mesh/polygon_triangulate.cpp \
    src/post/frdreader.cpp \
    src/post/occmeshtoccxmesh.cpp \
    src/post/postengine.cpp \
    src/post/posttools.cpp \
    src/post/cubicequation.cpp \
    src/mesh/stldoctor.cpp \
    src/viewer/qprogressindicator.cpp \
    src/viewer/qoccprogressindicator.cpp \
    src/mesh/extendedrwstl.cpp \
    src/mesh/ng_meshvs_deformeddatasource2d.cpp \
    src/mesh/occmesher.cpp \
    src/viewer/scaleselector.cpp \
    src/mesh/simplifymesh.cpp \
    src/mesh/curvature.cpp \
    src/mesh/prismaticlayer.cpp \
    src/gui/itemselector.cpp \
    src/mesh/netgentools.cpp \
    src/gui/qfileselect.cpp \
    src/memory/memoryprofiler.cpp \
    src/geometry/geometryhealing.cpp \
    src/mesh/igtools.cpp \
    src/mesh/tetwildmesher.cpp \
    src/mesh/mshconvert.cpp \
    src/mesh/meshuvprojection.cpp \
    src/mesh/customMesher/custommesher.cpp \
    src/mesh/surfacemeshcutter.cpp \
    src/mesh/datasourcebuildercontroller.cpp \
    src/mesh/surfacemeshsmoother.cpp \
    src/mesh/nettgentool2.cpp \
    src/mesh/qmorph.cpp \
    src/mesh/triangleaccessor.cpp \
    src/geometry/TriDist.cpp \
    src/gui/usermessagesmodel.cpp \
    src/geometry/geometryface.cpp \
    src/geometry/OCCface.cpp \
    src/geometry/meshface.cpp \
    src/mesh/tetsplitter.cpp \
    src/mesh/tethex.cpp \
    src/mesh/meshselectnlayers.cpp \
    src/viewer/wbtrihedron.cpp \
    src/mesh/meshslicer.cpp \
    src/mesh/silcemeshvs_mesh.cpp \
    src/mesh/tetqualityclass.cpp \
    src/viewer/qhistogram.cpp \
    src/electrostatic/poissonsolver.cpp \
    src/electrostatic/particlesinfieldssolver.cpp \
    src/electrostatic/particlesemitter.cpp \
    src/mesh/meshnodesrenumberingtool.cpp \
    src/mesh/isostripbuilder.cpp \
    src/mesh/rayintersectmesh.cpp \
    src/post/rainflow.cpp \
    src/mesh/datasourcebuilder.cpp \
    src/mesh/smoothingtools.cpp \
    src/mesh/pointtomeshdistance.cpp \
    src/mesh/surfacemeshtofacemeshes.cpp \
    ext/occ_extended/ais_arrowmarker.cpp \
    ext/occ_extended/ais_colorscaleextended.cpp \
    ext/occ_extended/ais_cormarker.cpp \
    ext/occ_extended/ais_curvedarrowmarker.cpp \
    ext/occ_extended/ais_customtrihedron.cpp \
    ext/occ_extended/ais_doublearrowmarker.cpp \
    ext/occ_extended/ais_extendedshape.cpp \
    ext/occ_extended/ais_meshsegmentmarker.cpp \
    ext/occ_extended/ais_spheremarker.cpp \
    ext/occ_extended/arrayofcolors.cpp \
    ext/pymesh/MshLoader.cpp \
    ext/pymesh/MshSaver.cpp \
    src/controller/solutionworker.cpp \
    src/connections/contactfinder.cpp \
    src/connections/contactparameters.cpp \
    src/databases/geometrydatabase.cpp \
    src/databases/meshdatabase.cpp \
    src/databases/simulationdatabase.cpp \
    src/main/main.cpp \
    src/main/maintreetools.cpp \
    src/main/mainwindow.cpp \
    src/mesh/isostrip.cpp \
    src/post/postobject.cpp \
    src/utils/ccout.cpp \
    src/utils/parser.cpp \
    src/viewer/dockableviewport.cpp \
    src/viewer/shapeselector.cpp \
    src/viewer/shapeselectorbox.cpp \
    src/controller/meshingserver.cpp \
    src/controller/meshworker.cpp \
    src/controller/modelloadercontroller.cpp \
    src/gui/optionsWidget/colorselector.cpp \
    src/gui/optionsWidget/optionswidget.cpp \
    src/gui/optionsWidget/qtablestandarditem.cpp \
    src/gui/optionsWidget/tabledelegate.cpp \
    src/gui/boundaryvaluemanager.cpp \
    src/gui/contextmenubuilder.cpp \
    src/gui/deserializerclass.cpp \
    src/gui/detailviewer.cpp \
    src/gui/directionselector.cpp \
    src/gui/generaldelegate.cpp \
    src/gui/geometrynodebuilder.cpp \
    src/gui/lineedit.cpp \
    src/gui/load.cpp \
    src/gui/markerbuilder.cpp \
    src/gui/markers.cpp \
    src/gui/meshselector.cpp \
    src/gui/meshtoolbar.cpp \
    src/gui/nodefactory.cpp \
    src/gui/property.cpp \
    src/gui/qextendedstandarditem.cpp \
    src/gui/qpushbuttonextended.cpp \
    src/gui/resultstoolbar.cpp \
    src/gui/serializerclass.cpp \
    src/gui/simulationnodeclass.cpp \
    src/gui/stepimporter.cpp \
    src/gui/stlapiwriter.cpp \
    src/gui/writelabelclass.cpp \
    src/main/simulationmanager.cpp \
    src/mapper/mapper3dclass.cpp \
    src/mapper/openfoamreader.cpp \
    src/utils/cliptool/clipTool.cpp \
    src/utils/cliptool/cliptooldelegate.cpp \
    src/utils/feaTool/bolttool.cpp \
    src/utils/exportingtools.cpp \
    src/utils/geomtoolsclass.cpp \
    src/utils/graphicstools.cpp \
    src/utils/meshtools.cpp \
    src/utils/modelloader.cpp \
    src/utils/steptools.cpp \
    src/utils/testtools.cpp \
    src/utils/textviewer.cpp \
    src/utils/tools.cpp \
    src/utils/topologytools.cpp \
    src/utils/vectortool.cpp \
    src/gui/systemConsole/SimulationMonitor.cpp \
    src/gui/systemConsole/systemconsole.cpp \
    src/gui/StepGenerationWidget/stepgenerationwidget.cpp \
    src/gui/tabularData/customtablemodel.cpp \
    src/gui/tabularData/tableviewclass.cpp \
    src/gui/tabularData/tableviewclassitemdelegate.cpp \
    src/gui/tabularData/tablewidget.cpp \
    ext/StlMesh/StlMesh.cxx \
    ext/StlMesh/StlMesh_Mesh.cxx \
    ext/StlMesh/StlMesh_MeshDomain.cxx \
    ext/StlMesh/StlMesh_MeshExplorer.cxx \
    ext/StlMesh/StlMesh_MeshTriangle.cxx \
    ext/StlMesh/StlTransfer.cxx \
    src/gui/simulationmanagerdelegate.cpp \
    src/controller/interpolatorcontroller.cpp \
    src/controller/openfoamcontroller.cpp \
    src/controller/ccxsolvermanager.cpp \
    src/post/convergencedatachart.cpp \
    src/solver/ccx/ccxconsoletofile.cpp \
    src/solver/ccx/consolereader.cpp \
    src/solver/ccx/writesolverfileclass.cpp \
    src/solver/OF/ofwrite.cpp \
    src/solver/inputfilegenerator.cpp \
    src/gui/tabularData/tabulardataviewerclass.cpp

HEADERS  += src/mesh/ng_meshvs_datasource1d.h \
    src/mesh/ng_meshvs_datasource2d.h \
    src/mesh/ng_meshvs_datasource3d.h \
    src/mesh/ng_meshvs_datasourceface.h \
    src/mesh/tetgenmesher.h \
    src/mesh/ng_mesher2.h \
    src/viewer/occGLwidget.h \
    src/viewer/occPreGLwidget.h \
    src/mesh/elementtypes.h \
    src/mesh/mesherclass.h \
    src/mesh/edgemeshdatasourcebuildesclass.h \
    src/viewer/occpostwidget.h \
    src/mesh/polygon_triangulate.h \
    src/post/frdreader.h \
    src/post/occmeshtoccxmesh.h \
    src/post/postengine.h \
    src/post/posttools.h \
    src/post/cubicequation.h \
    src/post/runterminationdata.h \
    src/mesh/stldoctor.h \
    src/viewer/qprogressindicator.h \
    src/viewer/qoccprogressindicator.h \
    src/mesh/extendedrwstl.h \
    src/mesh/ng_meshvs_deformeddatasource2d.h \
    src/mesh/occmesher.h \
    src/viewer/scaleselector.h \
    src/mesh/simplifymesh.h \
    src/mesh/curvature.h \
    src/mesh/prismaticlayer.h \
    src/mesh/meshpoint.h \
    src/mesh/mesh.h \
    src/mesh/meshelement2d.h \
    src/mesh/prismaticlayerparameters.h \
    src/gui/itemselector.h \
    src/mesh/netgentools.h \
    src/gui/qfileselect.h \
    src/memory/memoryprofiler.h \
    src/geometry/geometryhealing.h \
    src/mesh/igtools.h \
    src/geometry/geometrytag.h \
    src/geometry/polygon.h \
    src/mesh/tetwildmesher.h \
    src/mesh/mshconvert.h \
    src/mesh/mshconvert.h \
    src/mesh/meshuvprojection.h \
    src/mesh/indexedmapofmeshdatasources.h \
    src/mesh/customMesher/custommesher.h \
    src/viewer/resultpresentation.h \
    src/mesh/surfacemeshcutter.h \
    src/mesh/simulationdata.h \
    src/geometry/polyhedron.h \
    src/mesh/datasourcebuildercontroller.h \
    src/mesh/surfacemeshsmoother.h \
    src/mesh/nettgentool2.h \
    src/mesh/qmorph.h \
    src/mesh/triangleaccessor.h \
    src/geometry/triangletotriangleintersection.h \
    src/geometry/Tri.h \
    src/geometry/TriDist.h \
    src/geometry/PQP_Compile.h \
    src/geometry/MatVec.h \
    src/gui/usermessagesmodel.h \
    src/mesh/userMessage.h \
    src/geometry/geometryface.h \
    src/geometry/OCCface.h \
    src/geometry/meshface.h \
    src/mesh/tetsplitter.h \
    src/mesh/tethex.h \
    src/mesh/meshelementbycoords.h \
    src/mesh/meshselectnlayers.h \
    src/geometry/shapecomparison.h \
    src/viewer/wbtrihedron.h \
    src/mesh/meshslicer.h \
    src/mesh/slicemeshvs_mesh.h \
    src/mesh/tetqualityclass.h \
    src/viewer/qhistogramdata.h \
    src/viewer/qhistogram.h \
    src/electrostatic/poissonsolver.h \
    src/electrostatic/particle.h \
    src/electrostatic/particlesinfieldssolver.h \
    src/electrostatic/particlesemitter.h \
    src/mesh/renumberingtool.h \
    src/mesh/meshnodesrenumberingtool.h \
    src/mesh/isostripbuilder.h \
    src/mesh/rayintersectmesh.h \
    src/post/rainflow.h \
    src/mesh/datasourcebuilder.h \
    src/mesh/smoothingtools.h \
    src/mesh/pointtomeshdistance.h \
    src/mesh/surfacemeshtofacemeshes.h \
    src/connections/connectionpairgenerationoptions.h \
    src/connections/contactfinder.h \
    src/connections/contactparameters.h \
    src/databases/geometrydatabase.h \
    src/databases/meshdatabase.h \
    src/databases/simulationdatabase.h \
    src/geometry/PQP_Internal.h \
    src/main/maintreetools.h \
    src/main/mainwindow.h \
    src/main/mydefines.h \
    src/mesh/cgal_tools.h \
    src/mesh/isostrip.h \
    src/post/postobject.h \
    src/registeredMetatypes/meshvs_mesh_handle_reg.h \
    src/utils/ccout.h \
    src/utils/mathtools.h \
    src/utils/myconstant.h \
    src/utils/myenumvariables.h \
    src/utils/mystdlib.h \
    src/utils/nr3.h \
    src/utils/parser.h \
    src/viewer/dockableviewport.h \
    src/viewer/shapeselector.h \
    src/viewer/shapeselectorbox.h \
    ext/occ_extended/ais_arrowmarker.h \
    ext/occ_extended/ais_colorscaleextended.h \
    ext/occ_extended/ais_cormarker.h \
    ext/occ_extended/ais_curvedarrowmarker.h \
    ext/occ_extended/ais_customsignatures.h \
    ext/occ_extended/ais_customtrihedron.h \
    ext/occ_extended/ais_doublearrowmarker.h \
    ext/occ_extended/ais_extendedshape.h \
    ext/occ_extended/ais_meshsegmentmarker.h \
    ext/occ_extended/ais_spheremarker.h \
    ext/occ_extended/arrayofcolors.h \
    ext/occ_extended/handle_ais_arrowmarker_reg.h \
    ext/occ_extended/handle_ais_coloredshape_reg.h \
    ext/occ_extended/handle_ais_curvedarrowmarker_reg.h \
    ext/occ_extended/handle_ais_customtrihedron_reg.h \
    ext/occ_extended/handle_ais_doublearrowmarker_reg.h \
    ext/occ_extended/handle_ais_spheremarker_reg.h \
    ext/occ_extended/handle_ais_trihedron_reg.h \
    ext/occ_extended/handle_meshvs_datasource_reg.h \
    ext/occ_extended/handle_qoccprogressindicator_reg.h \
    src/connections/prebuiltcontactoptions.h \
    src/controller/meshingserver.h \
    src/controller/meshworker.h \
    src/controller/modelloadercontroller.h \
    src/gui/optionsWidget/celldata.h \
    src/gui/optionsWidget/colorselector.h \
    src/gui/optionsWidget/optionswidget.h \
    src/gui/optionsWidget/qtablestandarditem.h \
    src/gui/optionsWidget/tabledelegate.h \
    src/gui/optionsWidget/viewoptions.h \
    src/gui/actions3d.h \
    src/gui/boundaryvaluemanager.h \
    src/gui/contextmenubuilder.h \
    src/gui/deserializerclass.h \
    src/gui/detailviewer.h \
    src/gui/directionselector.h \
    src/gui/generaldelegate.h \
    src/gui/geometrynodebuilder.h \
    src/gui/lineedit.h \
    src/gui/load.h \
    src/gui/markerbuilder.h \
    src/gui/markers.h \
    src/gui/meshselector.h \
    src/gui/meshtoolbar.h \
    src/gui/nodefactory.h \
    src/gui/noteditabledelegate.h \
    src/gui/property.h \
    src/gui/qbackgroundevent.h \
    src/gui/qconsoleevent.h \
    src/gui/qextendedstandarditem.h \
    src/gui/qprogressevent.h \
    src/gui/qpushbuttonextended.h \
    src/gui/qsimulationstatusevent.h \
    src/gui/qtabwidgetextended.h \
    src/gui/resultstoolbar.h \
    src/gui/se_exception.h \
    src/gui/selectionmodes.h \
    src/gui/serializerclass.h \
    src/gui/simulationnodeclass.h \
    src/gui/stepimporter.h \
    src/gui/stlapiwriter.h \
    src/gui/workingmode.h \
    src/gui/writelabelclass.h \
    src/main/simulationmanager.h \
    src/main/ui_mainwindow.h \
    src/mapper/mapper3dclass.h \
    src/mapper/openfoamreader.h \
    src/registeredMetatypes/listofmesh.h \
    src/registeredMetatypes/listofshape.h \
    src/registeredMetatypes/mapofmeshdatasources.h \
    src/registeredMetatypes/topods_shape_reg.h \
    src/utils/cliptool/clipTool.h \
    src/utils/cliptool/cliptooldelegate.h \
    src/utils/feaTool/bolttool.h \
    src/utils/exportingtools.h \
    src/utils/geomtoolsclass.h \
    src/utils/global.h \
    src/utils/graphicstools.h \
    src/utils/hash_c.h \
    src/utils/meshtools.h \
    src/utils/modelloader.h \
    src/utils/steptools.h \
    src/utils/testtools.h \
    src/utils/textviewer.h \
    src/utils/tools.h \
    src/utils/topologytools.h \
    src/utils/vectortool.h \
    src/viewer/displaymode.h \
    src/viewer/displayquality.h \
    src/viewer/occhandle.h \
    src/gui/systemConsole/SimulationMonitor.h \
    src/gui/systemConsole/systemconsole.h \
    src/gui/StepGenerationWidget/stepgenerationwidget.h \
    src/gui/tabularData/customtablemodel.h \
    src/gui/tabularData/tableviewclass.h \
    src/gui/tabularData/tableviewclassitemdelegate.h \
    src/gui/tabularData/tablewidget.h \
    src/gui/tabularData/tabulardatacolumns.h \
    src/utils/ccxtools.h \
    ext/StlMesh/NCollection_StlIterator.hxx \
    ext/StlMesh/StlAPI.hxx \
    ext/StlMesh/StlAPI_ErrorStatus.hxx \
    ext/StlMesh/StlAPI_Reader.hxx \
    ext/StlMesh/StlAPI_Writer.hxx \
    ext/StlMesh/StlMesh.hxx \
    ext/StlMesh/StlMesh_Mesh.hxx \
    ext/StlMesh/StlMesh_MeshDomain.hxx \
    ext/StlMesh/StlMesh_MeshExplorer.hxx \
    ext/StlMesh/StlMesh_MeshTriangle.hxx \
    ext/StlMesh/StlMesh_SequenceOfMesh.hxx \
    ext/StlMesh/StlMesh_SequenceOfMeshDomain.hxx \
    ext/StlMesh/StlMesh_SequenceOfMeshTriangle.hxx \
    ext/StlMesh/StlTransfer.hxx \
    ext/pymesh/Exception.h \
    ext/pymesh/MshLoader.h \
    ext/pymesh/MshSaver.h \
    src/gui/simulationmanagerdelegate.h \
    src/controller/interpolatorcontroller.h \
    src/controller/openfoamcontroller.h \
    src/controller/ccxsolvermanager.h \
    src/controller/solutionworker.h \
    src/post/convergencedatachart.h \
    src/solver/ccx/ccxconsoletofile.h \
    src/solver/ccx/ccxsolvermessage.h \
    src/solver/ccx/consolereader.h \
    src/solver/ccx/qccxsolvermessageevent.h \
    src/solver/ccx/solutioninfo.h \
    src/solver/ccx/writesolverfileclass.h \
    src/solver/OF/ofwrite.h \
    src/solver/inputfilegenerator.h \
    src/gui/tabularData/tabulardataviewerclass.h

FORMS    += src/main/mainwindow.ui

DEFINES += WNT  \
           OCCGEOMETRY \
           DISPLAYBOTTOMTABLES  \
           #MESH_DIAGNOSTICS    \
           _TURNONFPES_ \
           TETLIBRARY  \
           DEBUG_VERSION   \
           #COSTAMP_VERSION \
           GENERATE_FACE_MESH_DATASOURCES   \
           #USE_MESHFIX

DEFINES += QCUSTOMPLOT_USE_LIBRARY

#DEFINES += ONLY_MESHER
DEFINES += NEW_HASH

INCLUDEPATH = $$PWD/src/geometry \
              C:/OpenCASCADE7.3.0-vc14-64/opencascade-7.3.0/inc \
              $$PWD/ext/StlMesh    \
              "D:/Work/Costamp/OCC lib/EMESH_7.3.0_binaries_win64vc14/inc"     \
              "D:/Work/Costamp/OCC lib/OMF_7.3.0_binaries_win64vc14/inc"   \
              $$PWD/src/gui/optionsWidget    \
              $$PWD/src/connections \
              $$PWD/src/controller \
              $$PWD/src/geometry \
              $$PWD/src/main \
              $$PWD/src/utils \
              $$PWD/src/gui/tabularData    \
              $$PWD/src/gui/systemConsole    \
              $$PWD/src/gui/stepGenerationWidget    \
              $$PWD/src/solver    \
              $$PWD/src/solver/ccx    \
              $$PWD/src/solver/OF    \
              $$PWD/src/mesh   \
              $$PWD/src/viewer   \
              $$PWD/src/post \
              $$PWD/src/testTools    \
              $$PWD/src/mapper    \
              $$PWD/src/databases    \
              $$PWD/src/registeredMetatypes  \
              $$PWD/src/gui  \
              $$PWD/src/utils  \
              $$PWD/ext/eigen    \
              $$PWD    \
              C:/CGAL-4.14/include    \
              C:/local/boost_1_69_0   \                 # needed by cgal
              C:/CGAL-4.14/auxiliary/gmp/include    \   # needed by cgal
              $$PWD/ext/libigl   \
              $$PWD/ext/pymesh    \
              $$PWD/ext/occ_extended

LIBS += \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKBin   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKBinL   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKBinTObj   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKBinXCAF   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKBO   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKBool   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKBRep   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKCAF   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKCDF   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKD3DHost   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKDCAF   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKDraw   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKernel   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKFeat   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKFillet   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKG2d   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKG3d   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKGeomAlgo   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKGeomBase   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKHLR   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKIGES   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKIVtk   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKIVtkDraw   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKLCAF   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKMath   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKMesh   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKMeshVS   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKOffset   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKOpenGl   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKPrim   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKQADraw   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKService   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKShHealing   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKStd   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKStdL   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKSTEP   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKSTEP209   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKSTEPAttr   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKSTEPBase   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKSTL   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKTObj   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKTObjDRAW   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKTopAlgo   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKTopTest   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKV3d   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKVCAF   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKViewerTest   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKVRML   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXCAF   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXDEDRAW   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXDEIGES   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXDESTEP   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXMesh   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXml   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXmlL   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXmlTObj   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXmlXCAF   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXSBase   \
-lC:\OpenCASCADE7.3.0-vc14-64\opencascade-7.3.0\win64\vc14\lib\TKXSDRAW   \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOCCLicense"    \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMF"    \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFBase"    \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFCAF" \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFCAM" \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFCAMTest" \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFQM"  \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFQMTest"  \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFTest"    \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFVS"  \
-l"D:\Work\Costamp\OCC lib\OMF_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOMFXCAF" \
-l"D:\Work\Costamp\OCC lib\EMESH_7.3.0_binaries_win64vc14\win64\vc14\lib\TKOCCLicense"    \
-l"D:\Work\Costamp\OCC lib\EMESH_7.3.0_binaries_win64vc14\win64\vc14\lib\TKQMesh"    \
-l"D:\Work\Costamp\OCC lib\EMESH_7.3.0_binaries_win64vc14\win64\vc14\lib\TKEMeshTest"    \
-l"D:\Work\Costamp\OCC lib\EMESH_7.3.0_binaries_win64vc14\win64\vc14\lib\TKEMesh"    \
-l"D:\quazip-0.7.3\build-quazip-Desktop_Qt_5_8_0_MSVC2013_64bit-Release\release\quazip" \
-l"C:\local\boost_1_69_0\lib64-msvc-14.0\libboost_thread-vc140-mt-x64-1_69" \
-l"C:\CGAL-4.14\build_MSVC14.0\lib\Release\CGAL-vc140-mt-4.14"   \
-l"C:\CGAL-4.14\build_MSVC14.0\lib\Release\CGAL_Core-vc140-mt-4.14" \
#-l"C:\CGAL-4.14\build_MSVC14.0\lib\Debug\CGAL-vc140-mt-gd-4.14"   \
#-l"C:\CGAL-4.14\build_MSVC14.0\lib\Debug\CGAL_Core-vc140-mt-gd-4.14" \

RESOURCES += \
    resources/cursors.qrc \
    resources/icons.qrc \
    resources/stylesheets.qrc

#----------------------
# unhandled exceptions
#----------------------
QMAKE_CXXFLAGS_EXCEPTIONS_ON = /EHa
QMAKE_CXXFLAGS_STL_ON = /EHa
QMAKE_CXXFLAGS += /EHa

#----------------
# quazip library
#----------------
win32: LIBS += -L$$PWD/ext/quazip/bin/ -lquazip

INCLUDEPATH += $$PWD/ext/quazip/inc
DEPENDPATH += $$PWD/ext/quazip/inc

#-----------------------------------------------
# zlib - only need the .h file (see Quazip doc)
#-----------------------------------------------
win32: LIBS += -L$$PWD/ext/zlib/x64/lib/ -lzlib1

INCLUDEPATH += $$PWD/ext/zlib/include
DEPENDPATH += $$PWD/ext/zlib/include

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/ext/zlib/x64/lib/zlib1.lib
else:win32-g++: PRE_TARGETDEPS += $$PWD/ext/zlib/x64/lib/libzlib1.a

#-------
# nglib
#-------
win32: LIBS += -L$$PWD/ext/nglib/lib/ -lnglib

INCLUDEPATH += $$PWD/ext/nglib/inc
DEPENDPATH += $$PWD/ext/nglib/inc

# --------------
# Tetgen define
# --------------
win32: LIBS += -L$$PWD/ext/tetgen/bin/ -ltet

#----------------
# tetgen library
#----------------
INCLUDEPATH += $$PWD/ext/tetgen/inc
DEPENDPATH += $$PWD/ext/tetgen/inc

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/ext/tetgen/bin/tet.lib
else:win32-g++: PRE_TARGETDEPS += $$PWD/ext/tetgen/bin/libtet.a

win32: LIBS += -LC:/CGAL-4.14/auxiliary/gmp/lib/ -llibgmp-10

#----------------------------------------------------
# CGAL dependencies - this feature is not documented
#----------------------------------------------------
INCLUDEPATH += C:/CGAL-4.14/auxiliary/gmp/include
DEPENDPATH += C:/CGAL-4.14/auxiliary/gmp/include

win32:!win32-g++: PRE_TARGETDEPS += C:/CGAL-4.14/auxiliary/gmp/lib/libgmp-10.lib
else:win32-g++: PRE_TARGETDEPS += C:/CGAL-4.14/auxiliary/gmp/lib/liblibgmp-10.a

win32: LIBS += -LC:/CGAL-4.14/auxiliary/gmp/lib/ -llibmpfr-4

INCLUDEPATH += C:/CGAL-4.14/auxiliary/gmp/include
DEPENDPATH += C:/CGAL-4.14/auxiliary/gmp/include

win32:!win32-g++: PRE_TARGETDEPS += C:/CGAL-4.14/auxiliary/gmp/lib/libmpfr-4.lib
else:win32-g++: PRE_TARGETDEPS += C:/CGAL-4.14/auxiliary/gmp/lib/liblibmpfr-4.a

# ---------------
# TetWild mesher
# ---------------
win32: LIBS += -L$$PWD/ext/tetwild/lib/ -ltetwild
INCLUDEPATH += $$PWD/ext/tetwild/inc
DEPENDPATH += $$PWD/ext/tetwild/inc

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/ext/tetwild/lib/tetwild.lib
else:win32-g++: PRE_TARGETDEPS += $$PWD/ext/tetwild/lib/libtetwild.a

# ---------------------------------------------
# this is for generating a valid debug version
# ---------------------------------------------
QMAKE_CXXFLAGS += -bigobj

# --------------------
# QCustomPlot library
# --------------------
LIBS += -L$$PWD/ext/QCustomPlot/qcp/ -lqcustomplot2
#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/QCustomPlot/qcp/ -lqcustomplot2
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/QCustomPlot/qcp/ -lqcustomplot2d

INCLUDEPATH += $$PWD/ext/QCustomPlot/qcp
DEPENDPATH += $$PWD/ext/QCustomPlot/qcp

# -------
# embree
# -------
win32: LIBS += -L$$PWD/ext/Embree/embree-3.11.0.x64.vc14.windows/lib/ -lembree3
INCLUDEPATH += $$PWD/ext/Embree/embree-3.11.0.x64.vc14.windows/include/embree3
DEPENDPATH += $$PWD/ext/Embree/embree-3.11.0.x64.vc14.windows/include/embree3


