#include "consolereader.h"
#include <windows.h>
#include <iostream>

#include <QWinEventNotifier>
#include <QDebug>

//! ----------------------------------------------------------
//! function: constructor
//! details:
//! ----------------------------------------------------------
consoleReader::consoleReader(QObject *parent):QObject(parent)
{
    we = new QWinEventNotifier(GetStdHandle(STD_OUTPUT_HANDLE),this);
    connect(we,SIGNAL(activated(HANDLE)),this,SLOT(notify()));
}

//! ----------------------------------------------------------
//! function: notify
//! details:
//! ----------------------------------------------------------
void consoleReader::notify()
{
        std::string line;
        std::getline(std::cin, line);
        const QString &s = QString::fromStdString(line);
        emit textReceived(s);
}
