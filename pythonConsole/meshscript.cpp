#include "meshscript.h"
#include <QString>
#include <QObject>
#include "simulationmanager.h"
#include "generaldelegate.h"
#include "property.h"
#include "pythonconsole.h"
//#include "src/geometry/geometrytag.h"

MeshScript::MeshScript(QObject *parent) : QObject(parent)
{
    mw = static_cast<MainWindow*>(this->parent());
}

//Item & node creation

void MeshScript::generateMesh()
{
    mw->getSimulationManager()->buildVolumeMesh();
}

void MeshScript::generateSurfaceMesh()
{
    mw->getSimulationManager()->buildSurfaceMesh();
}

void MeshScript::insertMethod(QString name)
{
    if (name.isEmpty()) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in insertMethod: name not valid\n");
        return;
    }
    mw->getSimulationManager()->handleItem(74);
    QExtendedStandardItem * item = this->getMeshItem("Method");
    QVariant data;
    data.setValue(name);
    item->setData(data,Qt::DisplayRole);
}

void MeshScript::insertMeshType(QString name)
{
    if (name.isEmpty()) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in insertMeshType: name not valid\n");
        return;
    }
    mw->getSimulationManager()->handleItem(79);
    QExtendedStandardItem * item = this->getMeshItem("Mesh type");
    QVariant data;
    data.setValue(name);
    item->setData(data,Qt::DisplayRole);
}

void MeshScript::insertBodySizing(QString name)
{
    if (name.isEmpty()) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in insertBodySizing: name not valid\n");
        return;
    }
    mw->getSimulationManager()->handleItem(51);
    QExtendedStandardItem * item = this->getMeshItem("Body sizing");
    QVariant data;
    data.setValue(name);
    item->setData(data,Qt::DisplayRole);
}

void MeshScript::insertFaceSizing(QString name)
{
    if (name.isEmpty()) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in insertFaceSizing: name not valid\n");
        return;
    }
    mw->getSimulationManager()->handleItem(72);
    QExtendedStandardItem * item = this->getMeshItem("Face sizing");
    QVariant data;
    data.setValue(name);
    item->setData(data,Qt::DisplayRole);
}

void MeshScript::insertEdgeSizing(QString name)
{
    if (name.isEmpty()) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in insertEdgeSizing: name not valid\n");
        return;
    }
    mw->getSimulationManager()->handleItem(68);
    QExtendedStandardItem * item = this->getMeshItem("Edge sizing");
    QVariant data;
    data.setValue(name);
    item->setData(data,Qt::DisplayRole);
}

void MeshScript::insertVertexSizing(QString name)
{
    if (name.isEmpty()) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in insertVertexSizing: name not valid\n");
        return;
    }
    mw->getSimulationManager()->handleItem(69);
    QExtendedStandardItem * item = this->getMeshItem("Vertex sizing");
    QVariant data;
    data.setValue(name);
    item->setData(data,Qt::DisplayRole);
}

void MeshScript::insertPrismaticLayer(QString name)
{
    if (name.isEmpty()) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in insertPrismaticLayer: name not valid\n");
        return;
    }
    mw->getSimulationManager()->handleItem(75);
    QExtendedStandardItem * item = this->getMeshItem("Prismatic layer");
    QVariant data;
    data.setValue(name);
    item->setData(data,Qt::DisplayRole);
}

void MeshScript::insertMeshMetric(QString name)
{
    if (name.isEmpty()) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in insertMeshMetric: name not valid\n");
        return;
    }
    mw->getSimulationManager()->handleItem(89);
    QExtendedStandardItem * item = this->getMeshItem("Mesh metric");
    QVariant data;
    data.setValue(name);
    item->setData(data,Qt::DisplayRole);
}

void MeshScript::setShowMeshNodes(QString value)
{
    SimulationNodeClass *aNode = mw->getSimulationManager()->Mesh_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    if (value.compare("Active") == 0) data.setValue(1);
    else if (value.compare("Off") == 0) data.setValue(0);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setShowMeshNodes: invalid value\n");
        return;
    }
    Property property_showMeshNodes("Show mesh nodes",data,Property::PropertyGroup_MeshViewOptions);
    aNode->replaceProperty("Show mesh nodes", property_showMeshNodes);
}

void MeshScript::setRelevance(int value)
{
    SimulationNodeClass *aNode = mw->getSimulationManager()->Mesh_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    if (value < -100 || value > 100) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setRelevance: invalid value\n");
        return;
    }
    Property property_relevance("Relevance",data,Property::PropertyGroup_Sizing);
    aNode->replaceProperty("Relevance", property_relevance);
}

void MeshScript::setInitialSizeSeed(QString value)
{
    SimulationNodeClass *aNode = mw->getSimulationManager()->Mesh_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    if (value.compare("Program controlled") == 0) data.setValue(0);
    else if (value.compare("Assembly") == 0) data.setValue(1);
    else if (value.compare("Part") == 0) data.setValue(2);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setInitialSizeSeed: invalid value\n");
        return;
    }
    Property property_initialSizeSeed("Initial size seed",data,Property::PropertyGroup_Sizing);
    aNode->replaceProperty("Initial size seed", property_initialSizeSeed);
}

void MeshScript::setSmoothing(QString value)
{
    SimulationNodeClass *aNode = mw->getSimulationManager()->Mesh_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    if (value.compare("Off") == 0) data.setValue(0);
    else if (value.compare("Low") == 0) data.setValue(1);
    else if (value.compare("Medium") == 0) data.setValue(2);
    else if (value.compare("High") == 0) data.setValue(3);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setSmoothing: invalid value\n");
        return;
    }
    Property property_smoothing("Smoothing",data,Property::PropertyGroup_Sizing);
    aNode->replaceProperty("Smoothing", property_smoothing);
}

void MeshScript::setElementMidsideNodes(QString value)
{
    SimulationNodeClass *aNode = mw->getSimulationManager()->Mesh_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    if (value.compare("Program controlled") == 0) data.setValue(0);
    else if (value.compare("Dropped") == 0) data.setValue(1);
    else if (value.compare("Kept") == 0) data.setValue(2);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setElementMidsideNodes: invalid value\n");
        return;
    }
     Property property_elementMidsideNodes("Element midside nodes", data, Property::PropertyGroup_Sizing);
    aNode->replaceProperty("Element midside nodes", property_elementMidsideNodes);
}

void MeshScript::setSubmeshes(QString value)
{
    SimulationNodeClass *aNode = mw->getSimulationManager()->Mesh_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    if (value.compare("Deferred") == 0) data.setValue(false);
    else if (value.compare("At meshing time") == 0) data.setValue(true);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setSubmeshes: invalid value\n");
        return;
    }
    Property property_submeshes("Submeshes",data,Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Submeshes", property_submeshes);
}

void MeshScript::setStraightSidedElements(QString value)
{
    SimulationNodeClass *aNode = mw->getSimulationManager()->Mesh_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    if (value.compare("Yes") == 0) data.setValue(true);
    else if (value.compare("No") == 0) data.setValue(false);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setStraightSidedElements: invalid value\n");
        return;
    }
    Property property_curvedElements("Straight sided elements",data, Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Straight sided elements", property_curvedElements);
}

//General

void MeshScript::setScopingMethod(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setScopingMethod: item not exist\n");
        return;
    }
    QVariant data;
    if (value.compare("Named selection") == 0) data.setValue(Property::ScopingMethod_NamedSelection);
    else if (value.compare("Geometry selection") == 0) data.setValue(Property::ScopingMethod_GeometrySelection);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setScopingMethod: invalid value\n");
        return;
    }
    Property prop_scopingMethod("Scoping method",data,Property::PropertyGroup_Scope);
    aNode->replaceProperty("Scoping method", prop_scopingMethod);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->scopingMethodChanged();
}

void MeshScript::setVolumeScope(QString nodeName, int pSNr)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setVolumeScope: item not exist\n");
        return;
    }
    Property::ScopingMethod verify = aNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
    if(verify != Property::ScopingMethod_GeometrySelection) return;
    GeometryTag loc;
    loc.isParent = true;
    loc.subShapeType = TopAbs_SOLID;
    loc.parentShapeNr = pSNr;
    loc.subTopNr = pSNr;
    std::vector<GeometryTag> vecLocs = aNode->getPropertyValue<std::vector<GeometryTag>>("Geometry");
    vecLocs.push_back(loc);
    QVariant data;
    data.setValue(vecLocs);
    Property prop_scope("Geometry",data,Property::PropertyGroup_Scope);
    Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
    aNode->replaceProperty("Geometry", prop_scope);
    aNode->replaceProperty("Tags", prop_tags);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->scopeChanged();
}

void MeshScript::setFaceScope(QString nodeName, int pSNr, int sTNr)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setFaceScope: item not exist\n");
        return;
    }
    Property::ScopingMethod verify = aNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
    if(verify != Property::ScopingMethod_GeometrySelection) return;
    GeometryTag loc;
    loc.isParent = false;
    loc.subShapeType = TopAbs_FACE;
    loc.parentShapeNr = pSNr;
    loc.subTopNr = sTNr;
    std::vector<GeometryTag> vecLocs = aNode->getPropertyValue<std::vector<GeometryTag>>("Geometry");
    vecLocs.push_back(loc);
    QVariant data;
    data.setValue(vecLocs);
    Property prop_scope("Geometry",data,Property::PropertyGroup_Scope);
    Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
    aNode->replaceProperty("Geometry", prop_scope);
    aNode->replaceProperty("Tags", prop_tags);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->scopeChanged();
}

void MeshScript::setEdgeScope(QString nodeName, int pSNr, int sTNr)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setEdgeScope: item not exist\n");
        return;
    }
    Property::ScopingMethod verify = aNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
    if(verify != Property::ScopingMethod_GeometrySelection) return;
    GeometryTag loc;
    loc.isParent = false;
    loc.subShapeType = TopAbs_EDGE;
    loc.parentShapeNr = pSNr;
    loc.subTopNr = sTNr;
    std::vector<GeometryTag> vecLocs = aNode->getPropertyValue<std::vector<GeometryTag>>("Geometry");
    vecLocs.push_back(loc);
    QVariant data;
    data.setValue(vecLocs);
    Property prop_scope("Geometry",data,Property::PropertyGroup_Scope);
    Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
    aNode->replaceProperty("Geometry", prop_scope);
    aNode->replaceProperty("Tags", prop_tags);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->scopeChanged();
}

void MeshScript::setNamedSelection(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setNamedSelection: item not exist\n");
        return;
    }
    Property::ScopingMethod verify = aNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
    if(verify != Property::ScopingMethod_NamedSelection) return;
    QExtendedStandardItem *nSItem = this->getNSItem(value);
    if (nSItem) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setNamedSelection: named selection not exist\n");
        return;
    }
    if(nSItem == nullptr) return;
    void *p = (void*)nSItem;
    QVariant data;
    data.setValue(p);
    Property prop_namedSelection("Named selection",data,Property::PropertyGroup_Scope);
    aNode->replaceProperty("Named selection", prop_namedSelection);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->scopeChanged();
}

void MeshScript::setSuppressed(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setSuppressed: item not exist\n");
        return;
    }
    QVariant data;
    bool state = false;
    if (value.compare("Active") == 0) {
        data.setValue(Property::SuppressionStatus_Active);
        state = true;
    }
    else if (value.compare("Suppressed") == 0) data.setValue(Property::SuppressionStatus_Suppressed);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setSuppressed: invalid value\n");
        return;
    }
    Property prop_suppressed("Suppressed",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Suppressed", prop_suppressed);
    if (state) emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->suppressionChanged(Property::SuppressionStatus_Active);
    else emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->suppressionChanged(Property::SuppressionStatus_Suppressed);
}

//Mesh method

void MeshScript::setPatchConforming(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setPatchConforming: item not recoverable\n");
        return;
    }
    bool var;
    if (value.compare("On") == 0) var = true;
    else if (value.compare("Off") == 0) var = false;
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setPatchConforming: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(var);
    Property prop_useBRep("Patch conforming",data,Property::PropertyGroup_Method);
    aNode->replaceProperty("Patch conforming", prop_useBRep);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->BRepFlagChanged();
}

void MeshScript::setTessellator(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setTessellator: item not recoverable\n");
        return;
    }
    bool verify = aNode->getPropertyValue<bool>("Patch conforming");
    if(verify) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setTessellator: action not possible\n");
        return;
    }
    QVariant data;
    if (value.compare("Express mesh") == 0) data.setValue(Property::meshEngine2D_OCC_ExpressMesh);
    else if (value.compare("Standard STL") == 0) data.setValue(Property::meshEngine2D_OCC_STL);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setTessellator: invalid value\n");
        return;
    }
    Property prop_surfaceDiscretizer("Tessellator",data,Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Tessellator", prop_surfaceDiscretizer);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->TessellatorChanged();
}

void MeshScript::setAngularDeflection(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setAngularDeflection: item not recoverable\n");
        return;
    }
    bool verify = aNode->getPropertyValue<bool>("Patch conforming");
    if (value < 0.001 || value > 0.7854) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setAngularDeflection: invalid value\n");
        return;
    }
    if (verify) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setAngularDeflection: action not possible\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_angularDeflection("Angular deflection",data,Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Angular deflection", prop_angularDeflection);
}

void MeshScript::setLinearDeflection(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setLinearDeflection: item not recoverable\n");
        return;
    }
    bool verify = aNode->getPropertyValue<bool>("Patch conforming");
    if (value < 0.001 || value > 100) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setLinearDeflection: invalid value\n");
        return;
    }
    if (verify) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setLinearDeflection: action not possible\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_linearDeflection("Linear deflection",data,Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Linear deflection", prop_linearDeflection);
}

void MeshScript::setMinFaceSize(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMinFaceSize: item not recoverable\n");
        return;
    }
    bool verify = aNode->getPropertyValue<bool>("Patch conforming");
    double max = aNode->getPropertyValue<double>("Max face size");
    Property::meshEngine2D Tessellator = aNode->getPropertyValue<Property::meshEngine2D>("Tessellator");
    if(value < 1e-6 || value > max) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMinFaceSize: invalid value\n");
        return;
    }
    if(verify || Tessellator != Property::meshEngine2D_OCC_ExpressMesh) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMinFaceSize: action not possible\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_minFaceSize("Min face size",data,Property::PropertyGroup_Method);
    aNode->replaceProperty("Min face size", prop_minFaceSize);
}

void MeshScript::setMaxFaceSize(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMaxFaceSize: item not recoverable\n");
        return;
    }
    bool verify = aNode->getPropertyValue<bool>("Patch conforming");
    double min = aNode->getPropertyValue<double>("Min face size");
    Property::meshEngine2D tessellator = aNode->getPropertyValue<Property::meshEngine2D>("Tessellator");
    if(value < min) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMaxFaceSize: invalid value\n");
        return;
    }
    if(verify || tessellator != Property::meshEngine2D_OCC_ExpressMesh) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMaxFaceSize: action not possible\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_maxFaceSize("Max face size",data,Property::PropertyGroup_Method);
    aNode->replaceProperty("Max face size", prop_maxFaceSize);
}

void MeshScript::setMeshOrder(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMeshOrder: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("First") == 0) data.setValue(Property::meshOrder_First);
    else if (value.compare("Second") == 0) data.setValue(Property::meshOrder_Second);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMeshOrder: invalid value\n");
        return;
    }
    Property prop_meshOrder("Mesh order",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Mesh order", prop_meshOrder);
}

void MeshScript::setSurfaceMesher(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setSurfaceMesher: item not recoverable\n");
        return;
    }
    bool verify = aNode->getPropertyValue<bool>("Patch conforming");
    Property::meshEngine3D volumeMesher = aNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");
    QVariant data;
    if (volumeMesher == Property::meshEngine3D_Tetgen_BR || volumeMesher == Property::meshEngine3D_TetWild) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setSurfaceMesher: action not possible\n");
        return;
    }
    if(verify) {
        if (value.compare("Netgen") == 0) data.setValue(Property::meshEngine2D_Netgen);
        else if (value.compare("Express mesh") == 0) data.setValue(Property::meshEngine2D_OCC_ExpressMesh);
        else {
            mw->myPythonConsole->pyQtScrCons->stdErr("error in setSurfaceMesher: invalid value\n");
            return;
        }
    } else {
        if (value.compare("Netgen STL") == 0) data.setValue(Property::meshEngine2D_Netgen_STL);
        else {
            mw->myPythonConsole->pyQtScrCons->stdErr("error in setSurfaceMesher: invalid value\n");
            return;
        }
    }
    Property prop_surfaceMesher("Surface mesher",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Surface mesher", prop_surfaceMesher);
}

void MeshScript::setVolumeMesher(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setVolumeMesher: item not recoverable\n");
        return;
    }
    bool verify = aNode->getPropertyValue<bool>("Patch conforming");
    QVariant data;
    if(verify) {
        if (value.compare("Netgen") == 0) data.setValue(Property::meshEngine3D_Netgen);
        else if (value.compare("Tetgen") == 0) data.setValue(Property::meshEngine3D_Tetgen);
        else {
            mw->myPythonConsole->pyQtScrCons->stdErr("error in setVolumeMesher: invalid value\n");
            return;
        }
    } else {
        if (value.compare("Netgen STL") == 0) data.setValue(Property::meshEngine3D_Netgen_STL);
        else if (value.compare("Tetgen") == 0) data.setValue(Property::meshEngine3D_Tetgen);
        else if (value.compare("Tetgen BR") == 0) data.setValue(Property::meshEngine3D_Tetgen_BR);
        else if (value.compare("Experimental mesher") == 0) data.setValue(Property::meshEngine3D_TetWild);
        else {
            mw->myPythonConsole->pyQtScrCons->stdErr("error in setVolumeMesher: invalid value\n");
            return;
        }
    }
    Property prop_volumeMesher("Volume mesher",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Volume mesher", prop_volumeMesher);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->volumeMesherChanged();
}

void MeshScript::setRunInMemory(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMethod) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setRunInMemory: item not recoverable\n");
        return;
    }
    Property::meshEngine3D volumeMesher = aNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");
    QVariant data;
    if(volumeMesher == Property::meshEngine3D_Tetgen || volumeMesher == Property::meshEngine3D_Tetgen_BR ||
            volumeMesher == Property::meshEngine3D_TetWild) {
        if (value.compare("Active") == 0) data.setValue(true);
        else if (value.compare("Off") == 0) data.setValue(false);
        else {
            mw->myPythonConsole->pyQtScrCons->stdErr("error in setRunInMemory: invalid value\n");
            return;
        }
    } else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setRunInMemory: action not possible\n");
        return;
    }
    Property prop_runInMemory("Run in memory",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Run in memory", prop_runInMemory);
}

//Mesh type

void MeshScript::setMeshType(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMeshType) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMeshType: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("Full tri-tet") == 0) data.setValue(0);
    else if (value.compare("Full quad-hexa") == 0) data.setValue(1);
    else if (value.compare("Quad-hexa dominant") == 0) data.setValue(2);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMeshType: invalid value\n");
        return;
    }
    Property prop_meshType("Mesh type",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Mesh type", prop_meshType);
}

//Body sizing

void MeshScript::setMinElementSize(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshBodyMeshControl) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMinElementSize: item not recoverable\n");
        return;
    }
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMinElementSize: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_minElSize("Min element size",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Min element size", prop_minElSize);
}

void MeshScript::setMaxElementSize(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshBodyMeshControl) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMaxElementSize: item not recoverable\n");
        return;
    }
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMaxElementSize: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_maxElSize("Max element size",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Max element size", prop_maxElSize);
}

void MeshScript::setGrading(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshBodyMeshControl) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setGrading: item not recoverable\n");
        return;
    }
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setGrading: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_grading("Grading",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Grading", prop_grading);
}

//Face sizing

void MeshScript::setFaceSizing(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshFaceSize) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setFaceSizing: item not recoverable\n");
        return;
    }
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setFaceSizing: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_faceMeshSizing("Face sizing",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Face sizing", prop_faceMeshSizing);
}

//Edge sizing

void MeshScript::setSizingType(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshEdgeSize) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setSizingType: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("Element size") == 0) data.setValue(0);
    else if (value.compare("Number of divisions") == 0) data.setValue(1);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setSizingType: invalid value\n");
        return;
    }
    Property prop_edgeSizingType("Sizing type",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Sizing type", prop_edgeSizingType);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->typeOfSizingChanged();

}

void MeshScript::setEdgeElementSize(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshEdgeSize) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setEdgeElementSize: item not recoverable\n");
        return;
    }
    int type = aNode->getPropertyValue<int>("Sizing type");
    if(type != 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setEdgeElementSize: action not possible\n");
        return;
    }
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setEdgeElementSize: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_edgeMeshSizing("Element size",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Element size", prop_edgeMeshSizing);
}

void MeshScript::setNumberOfDivisions(QString nodeName, int value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshEdgeSize) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setNumberOfDivisions: item not recoverable\n");
        return;
    }
    int type = aNode->getPropertyValue<int>("Sizing type");
    if(value < 1) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setNumberOfDivisions: invalid value\n");
        return;
    }
    if(type != 1) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setNumberOfDivisions: action not possible\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_NbDivisions("Number of divisions",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Number of divisions", prop_NbDivisions);
}

void MeshScript::setPinball(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshEdgeSize) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setPinball: item not recoverable\n");
        return;
    }
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setPinball: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_edgeSizingType("Pinball",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Pinball", prop_edgeSizingType);
}

void MeshScript::setVertexElementSize(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshEdgeSize) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setVertexElementSize: item not recoverable\n");
        return;
    }
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setVertexElementSize: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_edgeMeshSizing("Element size",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Element size", prop_edgeMeshSizing);
}

//Prismatic layer

void MeshScript::setBoundaryScopingMethod(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setBoundaryScopingMethod: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("Named selection") == 0) data.setValue(Property::ScopingMethod_NamedSelection);
    else if (value.compare("Geometry selection") == 0) data.setValue(Property::ScopingMethod_GeometrySelection);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setBoundaryScopingMethod: invalid value\n");
        return;
    }
    Property prop_boundaryScopingMethod("Boundary scoping method",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Boundary scoping method", prop_boundaryScopingMethod);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->boundaryScopingMethodChanged();
}

void MeshScript::setBoundary(QString nodeName, int pSNr, int sTNr)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setBoundary: item not recoverable\n");
        return;
    }
    Property::ScopingMethod verify = aNode->getPropertyValue<Property::ScopingMethod>("Boundary scoping method");
    if(verify != Property::ScopingMethod_GeometrySelection) return;
    GeometryTag loc;
    loc.isParent = false;
    loc.subShapeType = TopAbs_FACE;
    loc.parentShapeNr = pSNr;
    loc.subTopNr = sTNr;
    std::vector<GeometryTag> vecLocs = aNode->getPropertyValue<std::vector<GeometryTag>>("Geometry");
    vecLocs.push_back(loc);
    QVariant data;
    data.setValue(vecLocs);
    Property prop_boundaryScope("Boundary",data,Property::PropertyGroup_Definition);
    Property prop_boundaryTags("Boundary tags",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Boundary", prop_boundaryScope);
    aNode->replaceProperty("Boundary tags", prop_boundaryTags);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->boundaryScopeChanged();
}

void MeshScript::setBoundaryNamedSelection(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setBoundaryNamedSelection: item not recoverable\n");
        return;
    }
    Property::ScopingMethod verify = aNode->getPropertyValue<Property::ScopingMethod>("Boundary scoping method");
    if(verify != Property::ScopingMethod_NamedSelection) return;
    QExtendedStandardItem *nSItem = this->getNSItem(value);
    if (nSItem) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setBoundaryNamedSelection: named selection not exist\n");
        return;
    }
    void *p = (void*)nSItem;
    QVariant data;
    data.setValue(p);
    Property prop_namedSelection("Boundary named selection",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Boundary named selection", prop_namedSelection);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->boundaryScopeChanged();
}

void MeshScript::setOptions(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setOptions: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("First layer thickness") == 0) data.setValue(0);
    else if (value.compare("Total thickness") == 0) data.setValue(1);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setOptions: invalid value\n");
        return;
    }
    Property prop_options("Options",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Options", prop_options);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->prismaticLayerOptionsChanged();
}

void MeshScript::setNumberOfLayers(QString nodeName, int value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setNumberOfLayers: item not recoverable\n");
        return;
    }
    if(value < 1) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setNumberOfLayers: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_NbLayers("Number of layers",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Number of layers", prop_NbLayers);
}

void MeshScript::setFirstLayerHeight(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setFirstLayerHeight: item not recoverable\n");
        return;
    }
    int verify = aNode->getPropertyValue<int>("Options");
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setFirstLayerHeight: invalid value\n");
        return;
    }
    if(verify != 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setFirstLayerHeight: action not possible\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_firstLayerHeight("First layer height",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("First layer height", prop_firstLayerHeight);
}

void MeshScript::setTotalThickness(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setTotalThickness: item not recoverable\n");
        return;
    }
    int verify = aNode->getPropertyValue<int>("Options");
    if(value < 0) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setTotalThickness: invalid value\n");
        return;
    }
    if(verify != 1) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setTotalThickness: action not possible\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_totalThickness("Total thickness",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Total thickness", prop_totalThickness);
}

void MeshScript::setExpansionRatio(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setExpansionRatio: item not recoverable\n");
        return;
    }
    if(value < 1) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setExpansionRatio: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_growthRate("Expansion ratio",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Expansion ratio", prop_growthRate);
}

void MeshScript::setAlgorithm(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setAlgorithm: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("Pre") == 0) data.setValue(0);
    else if (value.compare("Post") == 0) data.setValue(1);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setAlgorithm: invalid value\n");
        return;
    }
    Property prop_generationAlgorithm("Algorithm",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Algorithm", prop_generationAlgorithm);
}

void MeshScript::setBoundaryMeshType(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setBoundaryMeshType: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("Hybrid") == 0) data.setValue(0);
    else if (value.compare("Tetrahedral") == 0) data.setValue(1);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setBoundaryMeshType: invalid value\n");
        return;
    }
    Property prop_boundaryMeshType("Boundary mesh type",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Boundary mesh type", prop_boundaryMeshType);
}

void MeshScript::setGuidingVectorsSmoothingSteps(QString nodeName, int value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setGuidingVectorsSmoothingSteps: item not recoverable\n");
        return;
    }
    if(value < 1) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setGuidingVectorsSmoothingSteps: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_NbGuidingVectorSmoothingSteps("Guiding vectors smoothing steps",data,Property::PropertyGroup_Smoothing);
    aNode->replaceProperty("Guiding vectors smoothing steps", prop_NbGuidingVectorSmoothingSteps);

}

void MeshScript::setThicknessSmoothingSteps(QString nodeName, int value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setThicknessSmoothingSteps: item not recoverable\n");
        return;
    }
    if(value < 1) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setThicknessSmoothingSteps: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_NbLayerThicknessSmoothingSteps("Thickness smoothing steps",data,Property::PropertyGroup_Smoothing);
    aNode->replaceProperty("Thickness smoothing steps", prop_NbLayerThicknessSmoothingSteps);
}

void MeshScript::setGVSCurvatureSensivity(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setGVSCurvatureSensivity: item not recoverable\n");
        return;
    }
    if(value < 1e-5) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setGVSCurvatureSensivity: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_curvatureSensitivityForGuidingVectorsSmoothing("Guiding vector smoothing - curvature sensitivity",data,Property::PropertyGroup_Smoothing);
    aNode->replaceProperty("Guiding vector smoothing - curvature sensitivity", prop_curvatureSensitivityForGuidingVectorsSmoothing);
}

void MeshScript::setTSCurvatureSensivity(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setTSCurvatureSensivity: item not recoverable\n");
        return;
    }
    if(value < 1e-5) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setTSCurvatureSensivity: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_curvatureSensitivityForThickessSmoothing("Thickness smoothing - curvature sensitivity",data,Property::PropertyGroup_Smoothing);
    aNode->replaceProperty("Thickness smoothing - curvature sensitivity", prop_curvatureSensitivityForThickessSmoothing);
}

void MeshScript::setCurvatureSensivity(QString nodeName, double value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setCurvatureSensivity: item not recoverable\n");
        return;
    }
    if(value < 1e-5) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setCurvatureSensivity: invalid value\n");
        return;
    }
    QVariant data;
    data.setValue(value);
    Property prop_curvatureSensitivityForShrink("Curvature sensitivity",data,Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Curvature sensitivity", prop_curvatureSensitivityForShrink);
}

void MeshScript::setLockBoundary(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setLockBoundary: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("Free") == 0) data.setValue(false);
    else if (value.compare("Locked") == 0) data.setValue(true);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setLockBoundary: invalid value\n");
        return;
    }
    Property prop_lockBoundary("Lock boundary",data,Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Lock boundary", prop_lockBoundary);
}

void MeshScript::setCheckSelfIntersections(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setCheckSelfIntersections: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("On") == 0) data.setValue(true);
    else if (value.compare("Off") == 0) data.setValue(false);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setCheckSelfIntersections: invalid value\n");
        return;
    }
    Property prop_checkSelfIntersection("Check self intersections",data,Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Check self intersections", prop_checkSelfIntersection);
}

void MeshScript::setCheckMutualIntersections(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshPrismaticLayer) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setCheckMutualIntersections: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("On") == 0) data.setValue(true);
    else if (value.compare("Off") == 0) data.setValue(false);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setCheckMutualIntersections: invalid value\n");
        return;
    }
    Property prop_checkMutualIntersections("Check mutual intersections",data,Property::PropertyGroup_Advanced);
    aNode->replaceProperty("Check mutual intersections", prop_checkMutualIntersections);
}

//Mesh metric

void MeshScript::setMetricType(QString nodeName, QString value)
{
    SimulationNodeClass *aNode = this->getNode(nodeName);
    if (aNode == nullptr || aNode->getType() != SimulationNodeClass::nodeType_meshMeshMetric) {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMetricType: item not recoverable\n");
        return;
    }
    QVariant data;
    if (value.compare("Modified mean ratio") == 0) data.setValue(0);
    else if (value.compare("Modified condition number") == 0) data.setValue(1);
    else if (value.compare("Modified volume-length") == 0) data.setValue(2);
    else {
        mw->myPythonConsole->pyQtScrCons->stdErr("error in setMetricType: invalid value\n");
        return;
    }
    Property prop_meshMetric("Metric type",data,Property::PropertyGroup_Definition);
    aNode->replaceProperty("Metric type", prop_meshMetric);
    emit mw->getSimulationManager()->myDetailViewer->myGeneralDelegate->meshMetricChanged();
}

//Utility

SimulationNodeClass *MeshScript::getNode(QString nodeName)
{
    QStandardItem *meshRootItem = mw->getSimulationManager()->Mesh_RootItem;
    int n = meshRootItem->rowCount();
    SimulationNodeClass *aNode = nullptr;
    for(int i = 0; i < n; i++){
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(meshRootItem->child(i));
        QString name = item->data(Qt::DisplayRole).value<QString>();
        if (nodeName.compare(name)==0) {
            aNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
            break;
        }
    }
    return aNode;
}

QExtendedStandardItem *MeshScript::getNSItem(QString nodeName)
{
    QStandardItem *namedSelectionRootItem = mw->getSimulationManager()->NamedSelection_RootItem;
    int n = namedSelectionRootItem->rowCount();
    QExtendedStandardItem *nSItem = nullptr;
    for(int i = 0; i < n; i++){
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(namedSelectionRootItem->child(i));
        QString name = item->data(Qt::DisplayRole).value<QString>();
        if (nodeName.compare(name)==0) {
            nSItem = item;
            break;
        }
    }
    return nSItem;
}

QExtendedStandardItem *MeshScript::getMeshItem(QString nodeName)
{
    QStandardItem *meshRootItem = mw->getSimulationManager()->Mesh_RootItem;
    int n = meshRootItem->rowCount();
    QExtendedStandardItem *meshItem = nullptr;
    for(int i = 0; i < n; i++){
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(meshRootItem->child(i));
        QString name = item->data(Qt::DisplayRole).value<QString>();
        if (nodeName.compare(name)==0) {
            meshItem = item;
            break;
        }
    }
    return meshItem;
}
