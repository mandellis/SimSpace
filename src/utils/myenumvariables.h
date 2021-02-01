#ifndef MYENUMVARIABLES
#define MYENUMVARIABLES

#include <QMetaType>
#include <QObject>

enum button
{
    button_selectSolid,
    button_selectFace,
    button_selectEdge,
    button_selectVertex
};

enum geometryImportStatus
{
    geometryImportStatus_Loaded,
    geometryImportStatus_NotLoaded
};

#endif // MYENUMVARIABLES
