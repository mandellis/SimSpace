#ifndef WORKINGMODE
#define WORKINGMODE

enum curWorkingMode
{
    curWorkingMode_onModel = 1,
    curWorkingMode_onCoordinateSystem = 2,
    curWorkingMode_onContact = 3,
    curWorkingMode_onNamedSelections = 4,
    curWorkingMode_onMesh = 5,
    curWorkingMode_onBC = 6,
    curWorkingMode_onSolution = 7,
    curWorkingMode_none = 8
};

#endif // WORKINGMODE
