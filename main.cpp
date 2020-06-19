#include "mainwindow.h"
#include <QApplication>
//#include "StdoutRedirector/widget.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWindow w;
    //w.setWindowFlags(Qt::WindowStaysOnTopHint);

    //! PLEASE do not use resizing functions here
    w.show();

    return a.exec();
}
