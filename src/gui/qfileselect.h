#ifndef QFILESELECT_H
#define QFILESELECT_H

#include <QWidget>
#include <QString>
#include <QLineEdit>

class QPushButton;
class QFileDialog;

class QFileSelect: public QWidget
{
    Q_OBJECT

public:

    QFileSelect(QString currentDir="", QWidget *parent=0);

    inline QString getText() { return mySelectedFile; }
    inline void setText(const QString &text) { myTextEdit->setText(text); }

private:

    QFileDialog *myFileDialog;
    QLineEdit *myTextEdit;
    QPushButton *myButton;

    QString myCurDir;
    QString mySelectedFile;

    void createContent();

private:

    QString getDirFromFile(const QString &filePath);

private slots:

    void browse();
    void updateText(QString text);
    void emitEditingFinished();

signals:

    void editingFinished();
};

#endif // QFILESELECT_H
