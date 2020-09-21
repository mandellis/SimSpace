#ifndef SELECTIONMODES
#define SELECTIONMODES

enum CurGlobalSelectionMode
{
    CurGlobalSelectionMode_Single,
    CurGlobalSelectionMode_Multiple
};

enum CurSelectionMode
{
    CurSelection_Nothing,
    CurSelection_Solid,
    CurSelection_Face,
    CurSelection_Edge,
    CurSelection_Vertex,
    CurSelection_PointCoordinatesPicking
};

enum SelectionType
{
    SelectionType_Geometry,
    SelectionType_Mesh
};

#endif // SELECTIONMODES

