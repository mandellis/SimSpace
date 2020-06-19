#ifndef MATERIALGENERATIONTOOL_H
#define MATERIALGENERATIONTOOL_H

#include <QWidget>

class MaterialGenerationTool : public QWidget
{
    Q_OBJECT

public:

    MaterialGenerationTool(QWidget *parent = 0);
    ~MaterialGenerationTool();

private:

    void createContent();

};

#endif // MATERIALGENERATIONTOOL_H
