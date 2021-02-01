#ifndef TEXTVIEWER_H
#define TEXTVIEWER_H

#include <QWidget>

class QPlainTextEdit;

class TextViewer: public QWidget
{
    Q_OBJECT

private:

    QPlainTextEdit *my_TextEdit;

public:

    TextViewer(QWidget *parent=0);

private slots:

    void Scroll();

public slots:

    void insertPlaintext(QString text);
    void clear();

};

#endif // TEXTVIEWER_H
