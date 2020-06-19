#include "usermessagesmodel.h"
#include <userMessage.h>

using namespace std;

//! -------------------------
//! function: constructor II
//! details:
//! -------------------------
userMessagesModel::userMessagesModel(const std::vector<userMessage> &messages, QAbstractTableModel *parent):
    QAbstractTableModel(parent),
    myMessages(messages)
{
    std::cout<<"userMessagesModel::userMessagesModel()->____constructor I called____"<<std::endl;
}

//! ------------------------------
//! function: constructor I
//! details:  needs setMessages()
//! ------------------------------
userMessagesModel::userMessagesModel(QObject *parent):
    QAbstractTableModel(parent),
    myMessages(std::vector<userMessage>())
{
    std::cout<<"userMessagesModel::userMessagesModel()->____constructor II called____"<<std::endl;
}

//! ----------------------
//! function: setMessages
//! details:
//! ----------------------
void userMessagesModel::setMessages(const std::vector<userMessage> &messages)
{
    myMessages = messages;
}
//! -------------------
//! function: rowCount
//! details:
//! -------------------
int userMessagesModel::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    int NbRows = int(myMessages.size());
    return NbRows;
}

//! ----------------------
//! function: columnCount
//! details:
//! ----------------------
int userMessagesModel::columnCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    return 3;   //! number of fields within the userMessage struct
}

//! ---------------------
//! function: headerData
//! details:
//! ---------------------
QVariant userMessagesModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    QVariant data;

    if (role!= Qt::DisplayRole) return QVariant();

    //! ------------------
    //! horizontal header
    //! ------------------
    if (orientation == Qt::Horizontal)
    {
        if (role == Qt::DisplayRole)
        {
            data.setValue(this->getHorizontalHeaderString(section));
        }
    }
    //! ----------------
    //! vertical header
    //! ----------------
    if(orientation == Qt::Vertical)
    {
        if(role == Qt::DisplayRole)
        {
            return section+1;
        }
    }
    return data;
}

//! ------------------------------------
//! function: getHorizontalHeaderString
//! details:  helper
//! ------------------------------------
QString userMessagesModel::getHorizontalHeaderString(int section) const
{
    QString header;
    switch(section)
    {
    case 0: header = QString("Time stamp"); break;
    case 1: header = QString("Status"); break;
    case 2: header = QString("Message"); break;
    }
    return header;
}

//! ---------------
//! function: data
//! details:
//! ---------------
QVariant userMessagesModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid()) return QVariant();

    if (role == Qt::DisplayRole)
    {
        int row = index.row();
        int col = index.column();
        QVariant data;

        const userMessage &message = myMessages.at(row);

        //! ------------------------------------------
        //! select the field of the message structure
        //! ------------------------------------------
        switch(col)
        {
        case 0: data.setValue(message.timeStamp); break;
        case 1: message.isDone==true? data.setValue(QString("OK")): data.setValue(QString("KO")); break;
        case 2: data.setValue(message.message); break;
        }
        return data;
    }
    if (role == Qt::TextAlignmentRole)
    {
        return Qt::AlignCenter;
    }
    return QVariant();
}

//! -------------------------------------
//! function: flags
//! details:  the items are not editable
//! -------------------------------------
Qt::ItemFlags userMessagesModel::flags(const QModelIndex &index) const
{
    return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
}

//! ------------------------
//! function: appendMessage
//! details:
//! ------------------------
bool userMessagesModel::appendMessage(const userMessage &aMessage)
{
    int pos = this->rowCount(QModelIndex());
    this->beginInsertRows(QModelIndex(),pos,pos);
    myMessages.push_back(aMessage);
    this->endInsertRows();

    //! -----------------------------------------------------
    //! inform the viewer about the position of the new data
    //! -----------------------------------------------------
    QModelIndex topLeftIndex = this->createIndex(pos,0);
    QModelIndex bottomRightIndex = this->createIndex(pos,0);

    emit dataChanged(topLeftIndex,bottomRightIndex);
    return true;
}

//! ------------------------
//! function: removeMessage
//! details:
//! ------------------------
bool userMessagesModel::removeMessage(int row)
{
    cout<<"userMessagesModel::removeMessage()->____function called____"<<endl;
    cout<<"userMessagesModel::removeMessage()->____initial number of rows: "<<this->rowCount(QModelIndex())<<"____"<<endl;

    //! -------------
    //! sanity check
    //! -------------
    if(row>this->rowCount(QModelIndex())-1 || row<0) return false;

    this->beginRemoveRows(QModelIndex(),row,row);
    myMessages.erase(myMessages.begin()+row);
    /*
    cout<<"userMessagesModel::removeMessage()->____final number of rows: "<<this->rowCount(QModelIndex())<<"____"<<endl;
    if(myMessages.size()==0)
    {
        cout<<"userMessagesModel::removeMessage()->____empty data____"<<endl;
    }
    else
    {
        cout<<"userMessagesModel::removeMessage()->____"<<myMessages.at(0).message.toStdString()<<"____"<<endl;
    }
    */
    this->endRemoveRows();

    //! -----------------------------------------------------
    //! inform the viewer about the position of the new data
    //! -----------------------------------------------------
    QModelIndex topLeftIndex = this->createIndex(row,0);
    QModelIndex bottomRightIndex = this->createIndex(row,0);
    emit dataChanged(topLeftIndex,bottomRightIndex);
    return true;
}

//! ----------------
//! function: clear
//! details:
//! ----------------
bool userMessagesModel::clear()
{
    this->beginRemoveRows(QModelIndex(),0,this->rowCount(QModelIndex())-1);
    myMessages.clear();
    this->endRemoveRows();

    //! -----------------------------------------------------
    //! inform the viewer about the position of the new data
    //! -----------------------------------------------------
    QModelIndex topLeftIndex = this->createIndex(0,0);
    QModelIndex bottomRightIndex = this->createIndex(0,0);

    emit dataChanged(topLeftIndex,bottomRightIndex);
    return true;
}

//! ------------------------
//! function: write
//! details:  write on disk
//! ------------------------
#include <fstream>
void userMessagesModel::write(const std::string &fullFilePath)
{
    if(fullFilePath.size()==0) return;
    std::fstream fs(fullFilePath);
    if(!fs.is_open()) return;

    //! -----------------------------
    //! write the number of messages
    //! -----------------------------
    int NbMessages = int(myMessages.size());
    fs<<NbMessages;

    //! -------------------
    //! write the messages
    //! -------------------
    for(int i=0; i<NbMessages; i++)
    {
        const userMessage &amsg = myMessages.at(i);
        fs<<amsg.timeStamp.toStdString()<<endl;
        fs<<amsg.isDone<<endl;
        fs<<amsg.message.toStdString()<<endl;
    }
}

//! ----------------------
//! function: constructor
//! details:  from disk
//! ----------------------
userMessagesModel::userMessagesModel(const std::string &fullFileName, QAbstractTableModel *parent):
    QAbstractTableModel(parent)
{
    myMessages = std::vector<userMessage>();
    std::vector<userMessage> messages;
    bool isDone = this->read(fullFileName,messages);
    if(isDone) myMessages = messages;
}

//! -------------------------------------------------------------
//! function: read
//! details:  read from disk - used in the constructor from disk
//! -------------------------------------------------------------
bool userMessagesModel::read(const std::string &fullFileName, std::vector<userMessage> &messages)
{
    if(fullFileName=="")
    {
        messages = std::vector<userMessage>();
        return false;
    }
    fstream fs(fullFileName);
    if(!fs.is_open())
    {
        messages = std::vector<userMessage>();
        return false;
    }

    //! ----------------------------
    //! read the number of messages
    //! ----------------------------
    int NbMessages;
    fs>>NbMessages;

    //! -------------
    //! sanity check
    //! -------------
    if(NbMessages<1) return false;

    //! ------------------
    //! read the messages
    //! ------------------
    for(int i=0; i<NbMessages; i++)
    {
        userMessage aMsg;
        std::string timeStampS;
        fs>>timeStampS;
        aMsg.timeStamp = QString::fromStdString(timeStampS);
        int isDoneInt;
        fs>>isDoneInt;
        isDoneInt==0? aMsg.isDone = false: aMsg.isDone = true;
        std::string msgS;
        fs>>msgS;
        aMsg.message = QString::fromStdString(msgS);
        messages.push_back(aMsg);
    }
    return true;
}
