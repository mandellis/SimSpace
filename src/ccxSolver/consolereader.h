#ifndef CONSOLEREADER_H
#define CONSOLEREADER_H

#include <QObject>
#include <QWinEventNotifier>
#include <windows.h>

class consoleReader: public QObject
{
    Q_OBJECT

public:

    explicit consoleReader(QObject *parent=0);

private:

    QWinEventNotifier *we;

private slots:

    void notify();

signals:

    void textReceived(const QString& text);
};

#endif // CONSOLEREADER_H
