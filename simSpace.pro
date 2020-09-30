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

SOURCES += main.cpp\
        mainwindow.cpp \
    ais_extendedshape.cpp \
    arrayofcolors.cpp \
    textviewer.cpp \
    simulationnodeclass.cpp \
    shapeselectorbox.cpp \
    property.cpp \
    generaldelegate.cpp \
    detailviewer.cpp \
    simulationmanager.cpp \
    qextendedstandarditem.cpp \
    shapeselector.cpp \
    writesolverfileclass.cpp \
    geomtoolsclass.cpp \
    load.cpp \
    customtablemodel.cpp \
    tableviewclass.cpp \
    tableviewclassitemdelegate.cpp \
    tablewidget.cpp \
    lineedit.cpp \
    directionselector.cpp \
    markers.cpp \
    vectortool.cpp \
    writelabelclass.cpp \
    serializerclass.cpp \
    deserializerclass.cpp \
    tools.cpp \
    nodefactory.cpp \
    qpushbuttonextended.cpp \
    exportingtools.cpp \
    topologytools.cpp \
    graphicstools.cpp \
    mapper3dclass.cpp \
    stepimporter.cpp \
    src/mesh/ng_meshvs_datasource1d.cpp \
    src/mesh/ng_meshvs_datasource2d.cpp \
    src/mesh/ng_meshvs_datasource3d.cpp \
    src/mesh/ng_meshvs_datasourceface.cpp \
    src/mesh/tetgenmesher.cpp \
    src/mesh/ng_mesher2.cpp \
    src/viewer/occGLwidget.cpp \
    src/viewer/occPreGLwidget.cpp \
    testtools.cpp \
    src/mesh/mesherclass.cpp \
    src/mesh/edgemeshdatasourcebuildesclass.cpp \
    ais_colorscaleextended.cpp \
    src/viewer/occpostwidget.cpp \
    stlapiwriter.cpp \
    openfoamreader.cpp \
    src/mesh/polygon_triangulate.cpp \
    contextmenubuilder.cpp \
    postobject.cpp \
    optionsWidget/colorselector.cpp \
    optionsWidget/optionswidget.cpp \
    optionsWidget/qtablestandarditem.cpp \
    optionsWidget/tabledelegate.cpp \
    meshingserver.cpp \
    meshworker.cpp \
    systemConsole/systemconsole.cpp \
    ccout.cpp \
    solutionworker.cpp \
    src/post/frdreader.cpp \
    src/post/occmeshtoccxmesh.cpp \
    ccxconsoletofile.cpp \
    src/post/postengine.cpp \
    src/post/posttools.cpp \
    src/post/cubicequation.cpp \
    parser.cpp \
    src/mesh/stldoctor.cpp \
    modelloader.cpp \
    src/viewer/qprogressindicator.cpp \
    src/viewer/qoccprogressindicator.cpp \
    modelloadercontroller.cpp \
    src/mesh/extendedrwstl.cpp \
    src/mesh/ng_meshvs_deformeddatasource2d.cpp \
    ais_cormarker.cpp \
    ais_meshsegmentmarker.cpp \
    systemConsole/SimulationMonitor.cpp \
    src/mesh/occmesher.cpp \
    resultstoolbar.cpp \
    src/viewer/scaleselector.cpp \
    dockableviewport.cpp \
    contactparameters.cpp \
    ais_spheremarker.cpp \
    ais_arrowmarker.cpp \
    src/cliptool/clipTool.cpp \
    src/cliptool/cliptooldelegate.cpp \
    src/mesh/simplifymesh.cpp \
    meshtoolbar.cpp \
    src/mesh/curvature.cpp \
    src/mesh/prismaticlayer.cpp \
    bolttool.cpp \
    ais_doublearrowmarker.cpp \
    markerbuilder.cpp \
    ais_customtrihedron.cpp \
    maintreetools.cpp \
    ais_curvedarrowmarker.cpp \
    databases/geometrydatabase.cpp \
    databases/meshdatabase.cpp \
    databases/simulationdatabase.cpp \
    steptools.cpp \
    geometrynodebuilder.cpp \
    simulationmanagerdelegate.cpp \
    src/gui/itemselector.cpp \
    meshtools.cpp \
    src/mesh/netgentools.cpp \
    src/gui/qfileselect.cpp \
    src/memory/memoryprofiler.cpp \
    src/geometry/geometryhealing.cpp \
    src/mesh/igtools.cpp \
    compatibility/StlMesh/StlMesh.cxx \
    compatibility/StlMesh/StlMesh_Mesh.cxx \
    compatibility/StlMesh/StlMesh_MeshDomain.cxx \
    compatibility/StlMesh/StlMesh_MeshExplorer.cxx \
    compatibility/StlMesh/StlMesh_MeshTriangle.cxx \
    src/mesh/tetwildmesher.cpp \
    src/mesh/mshconvert.cpp \
    src/mesh/backgroundmeshbuilder.cpp \
    src/mesh/meshuvprojection.cpp \
    src/mesh/facedatasourcebuilder.cpp \
    src/mesh/customMesher/custommesher.cpp \
    src/mesh/surfacemeshcutter.cpp \
    openfoamcontroller.cpp \
    interpolatorcontroller.cpp \
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
    contactfinder.cpp \
    src/mesh/tetqualityclass.cpp \
    src/viewer/qhistogram.cpp \
    tabulardataviewerclass1.cpp \
    src/post/convergencedatachart1.cpp \
    src/electrostatic/poissonsolver.cpp \
    src/electrostatic/particlesinfieldssolver.cpp \
    src/electrostatic/particlesemitter.cpp \
    src/mesh/meshnodesrenumberingtool.cpp \
    inputfilegenerator.cpp \
    ccxsolvermanager1.cpp \
    src/mesh/isostripbuilder.cpp \
    src/mesh/rayintersectmesh.cpp \
    meshselector.cpp \
    src/post/rainflow.cpp

HEADERS  += mainwindow.h \
    actions3d.h \
    selectionmodes.h \
    ais_extendedshape.h \
    displaymode.h \
    displayquality.h \
    workingmode.h \
    mydefines.h \
    myenumvariables.h \
    textviewer.h \
    simulationnodeclass.h \
    shapeselectorbox.h \
    property.h \
    listofshape.h \
    generaldelegate.h \
    detailviewer.h \
    simulationmanager.h \
    qextendedstandarditem.h \
    shapeselector.h \
    writesolverfileclass.h \
    geomtoolsclass.h \
    load.h \
    customtablemodel.h \
    tableviewclass.h \
    tableviewclassitemdelegate.h \
    tablewidget.h \
    lineedit.h \
    directionselector.h \
    markers.h \
    vectortool.h \
    writelabelclass.h \
    serializerclass.h \
    deserializerclass.h \
    tools.h \
    nodefactory.h \
    meshtools.h \
    qpushbuttonextended.h \
    exportingtools.h \
    topologytools.h \
    graphicstools.h \
    mapper3dclass.h \
    arrayofcolors.h \
    mathtools.h \
    nr3.h \
    stepimporter.h \
    se_exception.h \
    src/mesh/ng_meshvs_datasource1d.h \
    src/mesh/ng_meshvs_datasource2d.h \
    src/mesh/ng_meshvs_datasource3d.h \
    src/mesh/ng_meshvs_datasourceface.h \
    src/mesh/tetgenmesher.h \
    src/mesh/ng_mesher2.h \
    src/viewer/occGLwidget.h \
    src/viewer/occPreGLwidget.h \
    testtools.h \
    src/mesh/elementtypes.h \
    src/mesh/mesherclass.h \
    src/mesh/edgemeshdatasourcebuildesclass.h \
    ais_colorscaleextended.h \
    src/viewer/occpostwidget.h \
    listofmesh.h \
    stlapiwriter.h \
    openfoamreader.h \
    src/mesh/polygon_triangulate.h \
    contextmenubuilder.h \
    postobject.h \
    mapofmeshdatasources.h \
    optionsWidget/celldata.h \
    optionsWidget/colorselector.h \
    optionsWidget/optionswidget.h \
    optionsWidget/qtablestandarditem.h \
    optionsWidget/tabledelegate.h \
    optionsWidget/viewoptions.h \
    qprogressevent.h \
    qbackgroundevent.h \
    meshingserver.h \
    meshworker.h \
    systemConsole/systemconsole.h \
    ccout.h \
    qconsoleevent.h \
    solutionworker.h \
    src/post/frdreader.h \
    src/post/occmeshtoccxmesh.h \
    ccxconsoletofile.h \
    qtabwidgetextended.h \
    solutioninfo.h \
    src/post/postengine.h \
    src/post/posttools.h \
    src/post/cubicequation.h \
    qsimulationstatusevent.h \
    src/post/runterminationdata.h \
    parser.h \
    prebuiltcontactoptions.h \
    src/mesh/stldoctor.h \
    modelloader.h \
    src/viewer/qprogressindicator.h \
    src/viewer/qoccprogressindicator.h \
    modelloadercontroller.h \
    topods_shape_reg.h \
    handle_qoccprogressindicator_reg.h \
    src/mesh/extendedrwstl.h \
    src/mesh/ng_meshvs_deformeddatasource2d.h \
    ais_cormarker.h \
    ais_meshsegmentmarker.h \
    ais_extendedshape.h \
    prebuiltcontactoptions.h \
    systemConsole/SimulationMonitor.h \
    noteditabledelegate.h \
    src/mesh/occmesher.h \
    resultstoolbar.h \
    src/viewer/scaleselector.h \
    dockableviewport.h \
    tabulardatacolumns.h \
    contactparameters.h \
    ais_spheremarker.h \
    ais_arrowmarker.h \
    ais_customsignatures.h \
    src/cliptool/clipTool.h \
    src/cliptool/cliptooldelegate.h \
    src/mesh/simplifymesh.h \
    meshtoolbar.h \
    src/mesh/curvature.h \
    src/mesh/prismaticlayer.h \
    src/mesh/meshpoint.h \
    src/mesh/mesh.h \
    bolttool.h \
    src/mesh/meshelement2d.h \
    ais_doublearrowmarker.h \
    markerbuilder.h \
    ais_customtrihedron.h \
    myconstant.h \
    handle_ais_trihedron_reg.h \
    handle_ais_customtrihedron_reg.h \
    handle_ais_doublearrowmarker_reg.h \
    handle_ais_arrowmarker_reg.h \
    maintreetools.h \
    handle_ais_curvedarrowmarker_reg.h \
    ais_curvedarrowmarker.h \
    src/mesh/prismaticlayerparameters.h \
    databases/geometrydatabase.h \
    databases/meshdatabase.h \
    databases/simulationdatabase.h \
    meshvs_mesh_handle_reg.h \
    handle_ais_coloredshape_reg.h \
    handle_ais_spheremarker_reg.h \
    steptools.h \
    geometrynodebuilder.h \
    simulationmanagerdelegate.h \
    src/gui/itemselector.h \
    src/mesh/netgentools.h \
    src/gui/qfileselect.h \
    src/memory/memoryprofiler.h \
    src/geometry/geometryhealing.h \
    src/mesh/igtools.h \
    hash_c.h \
    occhandle.h \
    compatibility/StlMesh/NCollection_StlIterator.hxx \
    compatibility/StlMesh/StlAPI.hxx \
    compatibility/StlMesh/StlAPI_ErrorStatus.hxx \
    compatibility/StlMesh/StlAPI_Reader.hxx \
    compatibility/StlMesh/StlAPI_Writer.hxx \
    compatibility/StlMesh/StlMesh.hxx \
    compatibility/StlMesh/StlMesh_Mesh.hxx \
    compatibility/StlMesh/StlMesh_MeshDomain.hxx \
    compatibility/StlMesh/StlMesh_MeshExplorer.hxx \
    compatibility/StlMesh/StlMesh_MeshTriangle.hxx \
    compatibility/StlMesh/StlMesh_SequenceOfMesh.hxx \
    compatibility/StlMesh/StlMesh_SequenceOfMeshDomain.hxx \
    compatibility/StlMesh/StlMesh_SequenceOfMeshTriangle.hxx \
    src/mesh/cgal_tools.h \
    src/geometry/geometrytag.h \
    src/geometry/polygon.h \
    src/mesh/tetwildmesher.h \
    src/mesh/mshconvert.h \
    src/mesh/mshconvert.h \
    src/mesh/backgroundmeshbuilder.h \
    libigl/include/igl/copyleft/cgal/test.h \
    src/mesh/meshuvprojection.h \
    src/mesh/facedatasourcebuilder.h \
    src/mesh/indexedmapofmeshdatasources.h \
    src/mesh/customMesher/custommesher.h \
    connectionpairgenerationoptions.h \
    src/viewer/resultpresentation.h \
    src/mesh/surfacemeshcutter.h \
    openfoamcontroller.h \
    src/mesh/simulationdata.h \
    interpolatorcontroller.h \
    src/geometry/polyhedron.h \
    src/mesh/datasourcebuildercontroller.h \
    src/mesh/surfacemeshsmoother.h \
    src/mesh/nettgentool2.h \
    global.h \
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
    mystdlib.h \
    src/geometry/shapecomparison.h \
    src/viewer/wbtrihedron.h \
    src/mesh/meshslicer.h \
    src/mesh/slicemeshvs_mesh.h \
    contactfinder.h \
    src/mesh/tetqualityclass.h \
    src/viewer/qhistogramdata.h \
    src/viewer/qhistogram.h \
    qccxsolvermessageevent.h \
    ccxsolvermessage.h \
    ccxtools.h \
    tabulardataviewerclass1.h \
    src/post/convergencedatachart1.h \
    src/electrostatic/poissonsolver.h \
    src/electrostatic/particle.h \
    src/electrostatic/particlesinfieldssolver.h \
    src/electrostatic/particlesemitter.h \
    src/mesh/renumberingtool.h \
    src/mesh/meshnodesrenumberingtool.h \
    inputfilegenerator.h \
    ccxsolvermanager1.h \
    src/mesh/isostripbuilder.h \
    src/mesh/rayintersectmesh.h \
    meshselector.h \
    src/post/rainflow.h

FORMS    += mainwindow.ui

DEFINES += WNT  \
           OCCGEOMETRY \
           DISPLAYBOTTOMTABLES  \
           #MESH_DIAGNOSTICS    \
           _TURNONFPES_ \
           TETLIBRARY  \
           DEBUG_VERSION   \
           COSTAMP_VERSION \
           GENERATE_FACE_MESH_DATASOURCES   \
           #USE_MESHFIX

DEFINES += QCUSTOMPLOT_USE_LIBRARY

#DEFINES += ONLY_MESHER
DEFINES += NEW_HASH

INCLUDEPATH = D:/Work/Qt/SimSpace/src/geometry \
              C:/OpenCASCADE7.3.0-vc14-64/opencascade-7.3.0/inc \
              D:/Work/Qt/SimSpace/compatibility/StlMesh    \
              "D:/Work/Costamp/OCC lib/EMESH_7.3.0_binaries_win64vc14/inc"     \
              "D:/Work/Costamp/OCC lib/OMF_7.3.0_binaries_win64vc14/inc"   \
              D:/Work/Qt/SimSpace/src/mesh   \
              D:/Work/Qt/SimSpace/src/viewer   \
              D:/Work/Qt/SimSpace/optionsWidget    \
              D:/Work/Qt/SimSpace/src/post \
              D:/Work/Qt/SimSpace/src/testTools    \
              D:/Work/Qt/SimSpace/databases    \
              D:/Work/Qt/SimSpace/registeredMetatypes  \
              D:/Work/Qt/SimSpace/src/gui  \
              D:/Work/Qt/SimSpace/eigen    \
              D:/Work/Qt/SimSpace    \
              C:/CGAL-4.14/include    \
              C:/local/boost_1_69_0   \                 # needed by cgal
              C:/CGAL-4.14/auxiliary/gmp/include    \   # needed by cgal
              D:/Work/Qt/SimSpace/libigl   \
              D:/Work/Qt/SimSpace/src/pymesh

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
    icons.qrc \
    cursors.qrc \
    stylesheets.qrc \
    src/TimeStepBuilder/resources.qrc

#----------------------
# unhandled exceptions
#----------------------
QMAKE_CXXFLAGS_EXCEPTIONS_ON = /EHa
QMAKE_CXXFLAGS_STL_ON = /EHa
QMAKE_CXXFLAGS += /EHa

#----------------
# quazip library
#----------------
win32: LIBS += -L$$PWD/quazip/bin/ -lquazip

INCLUDEPATH += $$PWD/quazip/inc
DEPENDPATH += $$PWD/quazip/inc

#-----------------------------------------------
# zlib - only need the .h file (see Quazip doc)
#-----------------------------------------------
win32: LIBS += -L$$PWD/zlib/x64/lib/ -lzlib1

INCLUDEPATH += $$PWD/zlib/include
DEPENDPATH += $$PWD/zlib/include

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/zlib/x64/lib/zlib1.lib
else:win32-g++: PRE_TARGETDEPS += $$PWD/zlib/x64/lib/libzlib1.a

#-------
# nglib
#-------
win32: LIBS += -L$$PWD/nglib/lib/ -lnglib

INCLUDEPATH += $$PWD/nglib/inc
DEPENDPATH += $$PWD/nglib/inc

# --------------
# Tetgen define
# --------------
win32: LIBS += -L$$PWD/tetgen/bin/ -ltet

#----------------
# tetgen library
#----------------
INCLUDEPATH += $$PWD/tetgen/inc
DEPENDPATH += $$PWD/tetgen/inc

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/tetgen/bin/tet.lib
else:win32-g++: PRE_TARGETDEPS += $$PWD/tetgen/bin/libtet.a

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

DISTFILES += \
    compatibility/StlMesh/StlMesh_Mesh.lxx \
    compatibility/StlMesh/StlMesh_MeshDomain.lxx \
    compatibility/StlMesh/StlMesh_MeshExplorer.lxx

# ---------------
# TetWild mesher
# ---------------
win32: LIBS += -L$$PWD/tetwild/lib/ -ltetwild
INCLUDEPATH += $$PWD/tetwild/inc
DEPENDPATH += $$PWD/tetwild/inc

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/tetwild/lib/tetwild.lib
else:win32-g++: PRE_TARGETDEPS += $$PWD/tetwild/lib/libtetwild.a

# ---------------------------------------------
# this is for generating a valid debug version
# ---------------------------------------------
QMAKE_CXXFLAGS += -bigobj

# --------------------
# QCustomPlot library
# --------------------
LIBS += -L$$PWD/QCustomPlot/qcp/ -lqcustomplot2
#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/QCustomPlot/qcp/ -lqcustomplot2
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/QCustomPlot/qcp/ -lqcustomplot2d

INCLUDEPATH += $$PWD/QCustomPlot/qcp
DEPENDPATH += $$PWD/QCustomPlot/qcp

# --------------
# manifest file
# --------------
#win32:
#{
#QMAKE_POST_LINK += mt -nologo -manifest $$PWD/manifest.xml -outputresource:$$OUT_PWD/$$TARGET”.exe” $$escape_expand(\n\t)
#}
