#include "qfileselect.h"

#include <QHBoxLayout>
#include <QPushButton>
#include <QDir>
#include <QFileDialog>

#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
QFileSelect::QFileSelect(QString currentDir, QWidget *parent): QWidget(parent)
{
    if(currentDir=="") myCurDir = QDir::current().absolutePath();

    mySelectedFile="";

    this->createContent();
    connect(myButton,SIGNAL(clicked(bool)),this,SLOT(browse()));
    connect(myTextEdit,SIGNAL(textEdited(QString)),this,SLOT(updateText(QString)));
    connect(myTextEdit,SIGNAL(returnPressed()),this,SLOT(emitEditingFinished()));
}

//! ---------------------
//! function: updateText
//! details:
//! ---------------------
QString QFileSelect::getDirFromFile(const QString &filePath)
{
    QString tmp = filePath;
    int n = tmp.split("/").last().length();
    tmp.chop(n);
    return tmp;
}

//! ---------------------
//! function: updateText
//! details:
//! ---------------------
void QFileSelect::updateText(QString text)
{
    mySelectedFile = text;
    myCurDir = getDirFromFile(text);
    //cout<<"____file: "<<mySelectedFile.toStdString()<<"____"<<endl;
    //cout<<"____directory: "<<myCurDir.toStdString()<<"____"<<endl;
}

//! -----------------
//! function: browse
//! details:
//! -----------------
void QFileSelect::browse()
{
    //cout<<"____function called____"<<endl;

    if (myFileDialog->exec() && myFileDialog->result() == QDialog::Accepted)
    {
        QString selectedFile = myFileDialog->selectedFiles()[0];
        myTextEdit->setText(selectedFile);
        myCurDir = getDirFromFile(selectedFile);
        mySelectedFile = selectedFile;
        emitEditingFinished();
    }
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void QFileSelect::createContent()
{
    myFileDialog = new QFileDialog(this);

    QHBoxLayout *l = new QHBoxLayout();

    myTextEdit = new QLineEdit();
    myTextEdit->setContentsMargins(0,0,0,0);
    myTextEdit->setFrame(true);
    myTextEdit->setText(myCurDir);

    myButton = new QPushButton("...");
    myButton->setContentsMargins(0,0,0,0);
    myButton->setFixedHeight(myTextEdit->sizeHint().height()+2);
    myButton->setFixedWidth(24);

    l->addWidget(myTextEdit);
    l->addWidget(myButton);

    l->setContentsMargins(0,0,0,0);
    l->setMargin(0);
    l->setSpacing(0);
    this->setLayout(l);
}

//! ------------------------------
//! function: emitEditingFinished
//! details:
//! ------------------------------
void QFileSelect::emitEditingFinished()
{
    //cout<<"____editing finished____"<<endl;
    emit editingFinished();
}
