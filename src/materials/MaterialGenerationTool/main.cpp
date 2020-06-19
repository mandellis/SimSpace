#include "MaterialgenerationTool.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MaterialGenerationTool w;
    w.show();

    return a.exec();
}
