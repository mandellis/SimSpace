#ifndef INPUTFILEGENERATOR_H
#define INPUTFILEGENERATOR_H

#include <QObject>
#include <QThread>
#include <vector>

class QPushButton;
class QProgressIndicator;

class inputFileGenerator: public QThread
{
    Q_OBJECT

private:

    void run() Q_DECL_OVERRIDE;

public:

    inputFileGenerator(QObject *parent = 0);
    void setParameters(std::vector<void *> parameters);
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    std::vector<void*> myParameters;
    QString myInputFileName;
    QProgressIndicator *myProgressIndicator;

public slots:

    //! write a CCX input file
    bool writeCCX();

signals:

    //! emit input file written (true/false)
    void inputFileWritten(bool isDone);
};

#endif // INPUTFILEGENERATOR_H
