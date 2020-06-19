#ifndef CONTEXTMENUBUILDER_H
#define CONTEXTMENUBUILDER_H

//! Qt
#include <QModelIndex>
#include <QMenu>
#include <QAction>

//! OCC
#include <AIS_InteractiveContext.hxx>
#include <TopAbs_ShapeEnum.hxx>

class contextMenuBuilder
{

public:

    static void buildMeshContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);                   //!ok
    static void buildContactContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);                //!ok
    static void buildGeometryContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);               //!ok
    static void buildCoordinateSystemMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);              //!ok
    static void buildModelRootContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);              //!ok
    static void buildImportContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled= true);                       //!ok
    static void buildRemotePointContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);            //!ok
    static void buildNamedSelectionContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);         //!ok
    static void buildStructuralAnalysisContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);
    static void buildThermalAnalysisContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);
    static void buildParticlesInFieldsContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);
    static void buildStructuralSolutionContextMenu(QMenu *contextMenu, bool addCommonActions=true, bool isEnabled = true);     //!ok
    static void buildThermalResultsContextMenu(QMenu *contextMenu, bool addCommonActions = true, bool isEnabled = true);
    static void buildCombinedAnalysisResultsContextMenu(QMenu* contextMenu, bool addCommonActions, bool isEnabled);

    static void buildPostObjectContextMenu(QMenu *contextMenu);
    static void buildCommonActions(QMenu* contextMenu, bool isEnabled = true);

    //! -----------------
    //! "common" actions
    //! -----------------
    static void addActionDelete(QMenu *contextMenu);
    static void addActionDuplicate(QMenu *contextMenu);
    static void addActionSuppress(QMenu *contextMenu);
    static void addActionUnsuppress(QMenu *contextMenu);
    static void addActionRename(QMenu *contextMenu);
    static void addActionRenameBasedOnDefinition(QMenu *contextMenu);
    static void addActionSuppressAllOther(QMenu *contextMenu);
    static void addActionInvertSuppressionSet(QMenu *contextMenu);
    static void addActionUnsuppressAllBodies(QMenu* contextMenu);
    static void addActionDeleteAllChildrenItems(QMenu *contextMenu);
    static void addActionCreateNamedSelection(QMenu *contextMenu);
};

#endif // CONTEXTMENUBUILDER_H
