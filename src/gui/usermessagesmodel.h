#ifndef USERMESSAGESMODEL_H
#define USERMESSAGESMODEL_H

//! ---
//! Qt
//! ---
#include <QAbstractTableModel>
#include <QObject>
#include <QVariant>

//! ----------------
//! custom includes
//! ----------------
#include <userMessage.h>

//! ----
//! C++
//! ----
#include <vector>
#include <string>
#include <iostream>

class userMessagesModel: public QAbstractTableModel
{

    Q_OBJECT

public:

    //! ------------
    //! constructor
    //! ------------
    userMessagesModel(QObject *parent=0);
    userMessagesModel(const std::vector<userMessage> &messages, QAbstractTableModel *parent=0);

    //! ------------------------
    //! constructor - from disk
    //! ------------------------
    userMessagesModel(const std::string &fullFileName, QAbstractTableModel *parent=0);

    //! -------------
    //! set messages
    //! -------------
    void setMessages(const std::vector<userMessage> &messages);

private:

    std::vector<userMessage> myMessages;

    //! ---------------------------------------------
    //! read - used within the constructor from disk
    //! ---------------------------------------------
    bool read(const std::string &fullFileName, std::vector<userMessage> &messages);

public:

    int rowCount(const QModelIndex &parent) const;
    int columnCount(const QModelIndex &parent) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    QVariant data(const QModelIndex &index, int role) const;
    Qt::ItemFlags flags(const QModelIndex &index) const;

    bool appendMessage(const userMessage &aMessage);
    bool removeMessage(int row);
    bool clear();

    //! --------------
    //! write on disk
    //! --------------
    void write(const std::string &fullFilePath);

    //! -------
    //! helper
    //! -------
    QString getHorizontalHeaderString(int section) const;
};

#endif // USERMESSAGESMODEL_H
