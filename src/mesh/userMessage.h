#ifndef USERMESSAGE_H
#define USERMESSAGE_H

//! ---
//! Qt
//! ---
#include <QMetaType>
#include <QString>

//! ----
//! C++
//! ----
#include <string>
#include <chrono>
#include <ctime>

class userMessage
{
public:

    QString timeStamp;
    bool isDone;
    QString message;

    /*
    //! --------------------
    //! default constructor
    //! --------------------
    userMessage()
    {
        isDone = false;
        message = "";
        timeStamp = QString::fromStdString(getTimeStamp());
    }
    */

    //! ------------
    //! constructor
    //! ------------
    userMessage(bool status = false, const QString &msg = ""):isDone(status),message(msg)
    {
        isDone = status;
        message = msg;
        timeStamp = QString::fromStdString(getTimeStamp());
    }

    //! ----------------
    //! copy constuctor
    //! ----------------
    userMessage(const userMessage &other)
    {
        isDone = other.isDone;
        message = other.message;
        timeStamp = other.timeStamp;
    }

    //! -----------
    //! assignment
    //! -----------
    userMessage operator = (const userMessage &other)
    {
        isDone = other.isDone;
        message = other.message;
        timeStamp = other.timeStamp;
        return *this;
    }

    //! ------------
    //! operator ==
    //! ------------
    inline bool operator == (const userMessage &other)
    {
        if(isDone!=other.isDone) return false;
        if(message!=other.message) return false;
        if(timeStamp!=other.timeStamp) return false;
        return true;
    }

private:

    //! -----------
    //! time stamp
    //! -----------
    std::string getTimeStamp()
    {
        auto theNow = std::chrono::system_clock::now();
        std::time_t theNow_time = std::chrono::system_clock::to_time_t(theNow);
        std::string timeStamp(std::ctime(&theNow_time));
        return timeStamp;
    }
};

Q_DECLARE_METATYPE(userMessage)

#endif // USERMESSAGE_H
