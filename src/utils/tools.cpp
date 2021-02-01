//! ----------------
//! custom includes
//! ----------------
#include "tools.h"
#include "property.h"
#include "qextendedstandarditem.h"
#include "src/main/mydefines.h"

//! ---
//! Qt
//! ---
#include <QDir>
#include <QApplication>

//! ----
//! OCC
//! ----
#include <TopAbs_ShapeEnum.hxx>

//! ----
//! C++
//! ----
#include <set>
#include <vector>

tools::tools() { ; }

//! ----------------------------
//! function: changeIconOpacity
//! details:
//! ----------------------------
void tools::changeIconOpacity(QExtendedStandardItem *item, bool isOpaque)
{
    SimulationNodeClass *theCurNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    QIcon theCurrentIcon =item->getIcon(theCurNode->getType());
    QPixmap p;
    if(isOpaque) p =  theCurrentIcon.pixmap(QSize(32,32), QIcon::Disabled, QIcon::Off);
    else p =  theCurrentIcon.pixmap(QSize(32,32), QIcon::Normal, QIcon::On);
    QIcon theModifiedIcon(p);
    item->setIcon(theModifiedIcon);
}

//! ---------------------------------------------------
//! function: clearDir
//! details:  completely delete the content of a dir
//! ---------------------------------------------------
void tools::clearDir(const QString &path)
{
    cout<<"tools::clearDir()->____function called____"<<endl;
    QDir dir(path);
    if(QDir(dir).exists())
    {
        cout<<"tools::clearDir()->____the directory exists____"<<endl;
        dir.setFilter(QDir::NoDotAndDotDot | QDir::Files);
        foreach(QString dirItem, dir.entryList())
        {
            if(dir.remove(dirItem)==false)
            {
                cerr<<"tools::clearDir()->____cannot remove file: "<<dirItem.toStdString()<<"____"<<endl;
            }
            else
            {
                cout<<"tools::clearDir()->____removing file: "<<dirItem.toStdString()<<"____"<<endl;
            }
        }
        dir.setFilter(QDir::NoDotAndDotDot | QDir::Dirs);
        foreach(QString dirItem, dir.entryList())
        {
            QDir subDir(dir.absoluteFilePath(dirItem));
            if(subDir.removeRecursively()==false)
            {
                cerr<<"tools::clearDir()->____cannot remove directory: "<<subDir.absolutePath().toStdString()<<"____"<<endl;
            }
            else
            {
                cout<<"tools::clearDir()->____removing: "<<subDir.absolutePath().toStdString()<<"____"<<endl;
            }
        }
    }
    dir.cdUp();
}

//! -------------------------------------
//! function: writeVectorOfLocations
//! details:  write std::vector<GeometryTag>
//! -------------------------------------
void tools::writeVectorOfLocations(const std::vector<GeometryTag> &vecLocs, std::ofstream &os)
{
    //!write the number of components
    os<<vecLocs.size()<<endl;

    for(std::vector<GeometryTag>::const_iterator it=vecLocs.begin(); it!=vecLocs.end(); ++it)
    {
        const GeometryTag &loc = *it;
        os<<loc.parentShapeNr<<endl;
        os<<loc.subTopNr<<endl;
        os<<loc.isParent<<endl;
        os<<loc.subShapeType<<endl;
    }
}

//! ------------------------------------
//! function: readVectorOfLocations
//! details:  read std::vector<GeometryTag>
//! ------------------------------------
std::vector<GeometryTag> tools::readVectorOfLocations(std::ifstream &is)
{
    //! read the number of components
    int N;
    is>>N;

    int subShapeTypeInt;
    std::vector<GeometryTag> vecLocs;
    for(int i=0; i<N; i++)
    {
        GeometryTag loc;
        is>>loc.parentShapeNr;
        is>>loc.subTopNr;
        is>>loc.isParent;
        is>>subShapeTypeInt;
        TopAbs_ShapeEnum subShapeType= static_cast<TopAbs_ShapeEnum>(subShapeTypeInt);
        loc.subShapeType = subShapeType;
        vecLocs.push_back(loc);
    }
    return vecLocs;
}

//! ------------------------
//! function: writeQVariant
//! details:
//! ------------------------
void tools::writeQVariant(const QVariant &var, std::ofstream &os)
{
    //! retrieve the id of the QVariant
    int id = QMetaType::type(var.typeName());

    switch(id)
    {
    case 1:     //! bool
        os<<id<<std::endl;
        os<<var.toBool()<<std::endl;
        break;

    case 2:     //! int
        os<<id<<std::endl;
        os<<var.toInt()<<std::endl;
        break;

    case 6:     //! float
        os<<id<<std::endl;
        os<<var.toDouble()<<std::endl;
        break;

    case 10:     //! QString
        os<<id<<'\0';
        os<<var.toString().toStdString()<<std::endl;
        break;

    case 34:    //! char
    {
        os<<id<<std::endl;
        QChar c_qt = var.toChar();
        char c = c_qt.toLatin1();
        os<<c<<std::endl;
    }
        break;

    case 38:    //! double
        os<<id<<std::endl;
        os<<var.toFloat()<<std::endl;
        break;
    }
}

//! -----------------------
//! function: readQVariant
//! details:
//! -----------------------
QVariant tools::readQVariant(std::ifstream &is)
{
    //cout<<"readQVariant->____function called____"<<endl;

    //! retrieve the id of the QVariant datum
    int id;
    is>>id;

    QVariant rVar;

    switch(id)
    {
    case 1:     //! bool
    {
        bool val;
        is>>val;
        rVar.setValue(val);
    }
        break;

    case 2:     //! int
    {
        int val;
        is>>val;
        rVar.setValue(val);
    }
        break;

    case 6:     //! double
    {
        double val;
        is>>val;
        rVar.setValue(val);
    }
        break;

    case 10:     //! QString
    {
        std::string val;
        std::getline(is,val);
        QString qval = QString::fromStdString(val);
        qval.remove(0,1);
        rVar.setValue(qval);
    }
        break;

    case 34:    //! char
    {
        char c;
        is>>c;
        QChar c_qt = c;
        rVar.setValue(c_qt);
    }
        break;

    case 38:    //! float
    {
        float val;
        is>>val;
        rVar.setValue(val);
    }
        break;
    }
    return rVar;
}

//! -------------------
//! function: getScope
//! details:
//! -------------------
std::vector<GeometryTag> tools::getScope(SimulationNodeClass *aNode)
{
    cout<<"tools::getScope->____function called____"<<endl;

    std::vector<GeometryTag> vecLoc;
    QExtendedStandardItem *item = aNode->getPropertyItem("Scoping method");
    if(item==Q_NULLPTR) return vecLoc;
    Property::ScopingMethod theScopingMethod = item->data(Qt::UserRole).value<Property>().getData().value<Property::ScopingMethod>();
    if(theScopingMethod == Property::ScopingMethod_GeometrySelection)
    {
        vecLoc = aNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
        cout<<"tools::getScope->____function called on \"Geometry\" scope: "<<vecLoc.size()<<" shapes____"<<endl;
    }
    else if(theScopingMethod ==Property::ScopingMethod_NamedSelection)
    {
        void *p = aNode->getPropertyValue<void*>("Named selection");
        QExtendedStandardItem* item1 = static_cast<QExtendedStandardItem*>(p);
        SimulationNodeClass *node1 = item1->data(Qt::UserRole).value<SimulationNodeClass*>();
        vecLoc = node1->getPropertyValue<std::vector<GeometryTag>>("Tags");
        cout<<"tools::getScope->____function called on \"Named selection\"_scope: "<<vecLoc.size()<<" shapes____"<<endl;
    }
    return vecLoc;
}

//! -------------------
//! function: getScope
//! details:
//! -------------------
std::vector<GeometryTag> tools::getScope(QExtendedStandardItem *item)
{
    std::vector<GeometryTag> vecLoc;    
    if(item==Q_NULLPTR) return vecLoc;
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    vecLoc = tools::getScope(node);
    return vecLoc;
}

//! ------------------------------
//! function: writeColor
//! details:  write in RGB format
//! ------------------------------
void tools::writeColor(QColor color, std::ofstream &os)
{
    int r,g,b;
    color.getRgb(&r,&g,&b);
    char s[32];
    sprintf(s,"%d\t%d\t%d\n",r,g,b);
    os<<s;
}

//! ------------------------------
//! function: readColor
//! details:  read the RGB format
//! ------------------------------
QColor tools::readColor(std::ifstream &is)
{
    int r,g,b;
    std::string val;
    std::getline(is,val);
    sscanf(val.c_str(),"%d%d%d",&r,&g,&b);
    QColor color(r,g,b);
    cout<<"tools::readColor___(R,G,B) = ("<<r<<", "<<g<<", "<<b<<")"<<endl;
    return color;
}

//! --------------------------
//! function: getWidgetByName
//! details:
//! --------------------------
QWidget* tools::getWidgetByName(const QString &widgetName)
{
    QWidgetList widgets = QApplication::allWidgets();
    QWidget *widget;
    for(int i=0; i<widgets.length(); i++)
    {
        widget = widgets.at(i);
        if(widget->objectName()==widgetName) return widget;
    }
    cout<<"tools::getWidgetByName->____"<<widgetName.toStdString()<<" not found____"<<endl;
    return Q_NULLPTR;
}

//! ------------------------
//! function: getWorkingDir
//! details:
//! ------------------------
QString tools::getWorkingDir()
{
    ifstream is;
    QString settingsFileName = QString(SYSTEM_PROGRAM_DATA).append("\\WB\\settings.txt");
    //cout<<"tools::getWorkingDir()->____Settings file name: "<<settingsFileName.toStdString()<<"____"<<endl;
    is.open(settingsFileName.toStdString());
    if(is.is_open())
    {
        std::string val;
        char tmp[512],wd[512];
        std::getline(is,val);
        is.close();
        sscanf(val.c_str(),"%s%s",tmp,wd);
        return QString::fromLatin1(wd);
    }
    return "";
}


//! --------------------------------------------------------
//! function: clearFromDuplicates
//! details:  eliminate the duplicated values from a vector
//! --------------------------------------------------------
std::vector<int> tools::clearFromDuplicates(const std::vector<int> &aVec)
{
    std::set<int> B(aVec.begin(), aVec.end());
    std::vector<int> vec;
    for(int n=0;n<B.size();n++)
    {
        int x = *std::next(B.begin(), n);
        vec.push_back(x);
    }
    return vec;
}

//! ---------------------------
//! function: timeStamp
//! details:  time stamp label
//! ---------------------------
QString tools::timeStamp()
{
    QDateTime dateTime;
    QString dateFormat = "dd/MM/yyyy";
    QString dateString = dateTime.currentDateTime().toString(dateFormat);
    QString timeFormat = "hh:mm";
    QString timeString = dateTime.currentDateTime().toString(timeFormat);
    QString timeStamp;
    timeStamp.append("Date: ").append(dateString).append("\n").append("Time: ").append(timeString);

    return timeStamp;
}

//! ------------------------------
//! function: getPathOfExecutable
//! details:
//! ------------------------------
#include <Windows.h>
std::string tools::getPathOfExecutable()
{
    LPWSTR buffer;
    GetModuleFileName(NULL, buffer, MAX_PATH);
    std::string aString;
    aString.reserve(wcslen(buffer));
    for (;*buffer; buffer++) aString += (char)*buffer;
    std::string::size_type pos = aString.find_last_of("\\/");
    if (pos == std::string::npos) return "";
    return aString.substr(0, pos);
}
