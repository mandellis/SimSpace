#include "textviewer.h"
#include <QPlainTextEdit>
#include <QVBoxLayout>

// constructor
TextViewer::TextViewer(QWidget *parent): QWidget(parent)
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    my_TextEdit = new QPlainTextEdit();
    my_TextEdit->setReadOnly(true);
    QFont my_font("Arial",12);
    my_TextEdit->setFont(my_font);

    layout->addWidget(my_TextEdit);
    connect(my_TextEdit,SIGNAL(textChanged()),this,SLOT(Scroll()));
}

// slot: move the cursor to the end
void TextViewer::Scroll()
{
    my_TextEdit->moveCursor(QTextCursor::End);
}

// slot: change the text
void TextViewer::insertPlaintext(QString text)
{
    my_TextEdit->insertPlainText(text);
}

// slot: clear
void TextViewer::clear()
{
    my_TextEdit->clear();
}
