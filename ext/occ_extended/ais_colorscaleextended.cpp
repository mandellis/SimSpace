#include "ais_colorscaleextended.h"

#include <AIS_InteractiveContext.hxx>
#include <Aspect_TypeOfColorScaleData.hxx>
#include <Aspect_TypeOfColorScalePosition.hxx>
#include <Aspect_Window.hxx>
#include <Geom_Line.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Graphic3d_ArrayOfPolygons.hxx>
#include <Graphic3d_ArrayOfPolylines.hxx>
#include <Graphic3d_AspectFillArea3d.hxx>
#include <Graphic3d_AspectText3d.hxx>
#include <Graphic3d_GraphicDriver.hxx>
#include <Graphic3d_ArrayOfTriangles.hxx>
#include <Prs3d_LineAspect.hxx>
#include <Prs3d_Root.hxx>
#include <Prs3d_ShadingAspect.hxx>
#include <Prs3d_Text.hxx>
#include <Prs3d_TextAspect.hxx>
#include <SelectMgr_EntityOwner.hxx>
#include <SelectMgr_Selection.hxx>
#include <Select3D_SensitiveBox.hxx>
#include <Select3D_SensitiveSegment.hxx>
#include <StdPrs_Curve.hxx>
#include <V3d_Viewer.hxx>
#include <V3d_View.hxx>


IMPLEMENT_STANDARD_RTTIEXT(AIS_ColorScaleExtended,AIS_InteractiveObject)

//=======================================================================
//function : AIS_ColorScaleExtended
//purpose  :
//=======================================================================
AIS_ColorScaleExtended::AIS_ColorScaleExtended() :
    myMin (0.0),
    myMax (1.0),
    myTitle (""),
    myFormat ("%.4g"),
    myInterval (10),
    myColorType (Aspect_TOCSD_AUTO),
    myLabelType (Aspect_TOCSD_AUTO),
    myAtBorder (Standard_True),
    myReversed (Standard_False),
    myIsLogarithmic (Standard_False),
    myLabelPos (Aspect_TOCSP_RIGHT),
    myTitlePos (Aspect_TOCSP_CENTER),
    myXPos (0),
    myYPos (0),
    myBreadth (0),
    myHeight (0),
    myTextHeight(20)
{
}

//=======================================================================
//function : GetRange
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::GetRange (Standard_Real& theMin, Standard_Real& theMax) const
{
    theMin = myMin;
    theMax = myMax;
}

//=======================================================================
//function : GetLabel
//purpose  :
//=======================================================================
TCollection_ExtendedString AIS_ColorScaleExtended::GetLabel (const Standard_Integer theIndex) const
{
    if (GetLabelType() == Aspect_TOCSD_USER)
    {
        if (theIndex <= 0 || theIndex > myLabels.Length())
        {
            return "";
        }
        return myLabels.Value (theIndex);
    }

    // value to be shown depends on label position
    Standard_Real aVal = IsLabelAtBorder() ? GetIntervalValue (theIndex - 1) :
                                             0.5 * (GetIntervalValue (theIndex - 1) + GetIntervalValue (theIndex));

    const TCollection_AsciiString aFormat = Format();
    Standard_Character aBuf[1024];
    sprintf (aBuf, aFormat.ToCString(), aVal);
    return TCollection_ExtendedString (aBuf);
}

//=======================================================================
//function : GetIntervalColor
//purpose  :
//=======================================================================
Quantity_Color AIS_ColorScaleExtended::GetIntervalColor (const Standard_Integer theIndex) const
{
    if (GetColorType() == Aspect_TOCSD_USER)
    {
        if (theIndex <= 0 || theIndex > myColors.Length())
        {
            return Quantity_Color();
        }
        return myColors.Value (theIndex);
    }

    return Quantity_Color (HueFromValue (theIndex - 1, 0, GetNumberOfIntervals() - 1), 1.0, 1.0, Quantity_TOC_HLS);
}

//=======================================================================
//function : GetLabels
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::GetLabels (TColStd_SequenceOfExtendedString& theLabels) const
{
    theLabels.Clear();
    for (Standard_Integer i = 1; i <= myLabels.Length(); i++)
        theLabels.Append (myLabels.Value (i));
}

//=======================================================================
//function : GetColors
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::GetColors (Aspect_SequenceOfColor& theColors) const
{
    theColors.Clear();
    for (Standard_Integer i = 1; i <= myColors.Length(); i++)
        theColors.Append (myColors.Value (i));
}

//=======================================================================
//function : SetMin
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetMin (const Standard_Real theMin)
{
    SetRange (theMin, GetMax());
}

//=======================================================================
//function : SetMax
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetMax (const Standard_Real theMax)
{
    SetRange (GetMin(), theMax);
}

//=======================================================================
//function : SetRange
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetRange (const Standard_Real theMin, const Standard_Real theMax)
{
    if (myMin == theMin && myMax == theMax)
        return;

    myMin = Min (theMin, theMax);
    myMax = Max (theMin, theMax);
}

//=======================================================================
//function : SetLabelType
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetLabelType (const Aspect_TypeOfColorScaleData theType)
{
    if (myLabelType == theType)
        return;

    myLabelType = theType;
}

//=======================================================================
//function : SetColorType
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetColorType (const Aspect_TypeOfColorScaleData theType)
{
    if (myColorType == theType)
        return;

    myColorType = theType;
}

//=======================================================================
//function : SetNumberOfIntervals
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetNumberOfIntervals (const Standard_Integer theNum)
{
    if (myInterval == theNum || theNum < 1)
        return;

    myInterval = theNum;
}

//=======================================================================
//function : SetTitle
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetTitle (const TCollection_ExtendedString& theTitle)
{
    if (myTitle == theTitle)
        return;

    myTitle = theTitle;
}

//=======================================================================
//function : SetFormat
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetFormat (const TCollection_AsciiString& theFormat)
{
    if (myFormat == theFormat)
        return;

    myFormat = theFormat;
}

//=======================================================================
//function : SetLabel
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetLabel (const TCollection_ExtendedString& theLabel, const Standard_Integer theIndex)
{
    Standard_Integer i = (theIndex <= 0 ? myLabels.Length() + 1 : theIndex);
    while (i > myLabels.Length())
        myLabels.Append (TCollection_ExtendedString());
    myLabels.SetValue (i, theLabel);
}

//=======================================================================
//function : SetIntervalColor
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetIntervalColor (const Quantity_Color& theColor, const Standard_Integer theIndex)
{
    Standard_Integer i = (theIndex <= 0 ? myColors.Length() + 1 : theIndex);
    while (i > myColors.Length())
        myColors.Append (Quantity_Color());
    myColors.SetValue (i, theColor);
}

//=======================================================================
//function : SetLabels
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetLabels (const TColStd_SequenceOfExtendedString& theSeq)
{
    myLabels.Clear();
    for (Standard_Integer i = 1; i <= theSeq.Length(); i++)
        myLabels.Append (theSeq.Value (i));
}

//=======================================================================
//function : SetColors
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetColors (const Aspect_SequenceOfColor& theSeq)
{
    myColors.Clear();
    for (Standard_Integer i = 1; i <= theSeq.Length(); i++)
        myColors.Append (theSeq.Value (i));
}

//=======================================================================
//function : SetLabelPosition
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetLabelPosition (const Aspect_TypeOfColorScalePosition thePos)
{
    if (myLabelPos == thePos)
        return;

    myLabelPos = thePos;
}

//=======================================================================
//function : SetTitlePosition
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetTitlePosition (const Aspect_TypeOfColorScalePosition thePos)
{
    if (myTitlePos == thePos)
        return;

    myTitlePos = thePos;
}

//=======================================================================
//function : SetReversed
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetReversed (const Standard_Boolean theReverse)
{
    if (myReversed == theReverse)
        return;

    myReversed = theReverse;
}

//=======================================================================
//function : SetLabelAtBorder
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetLabelAtBorder (const Standard_Boolean theOn)
{
    if (myAtBorder == theOn)
        return;

    myAtBorder = theOn;
}

//=======================================================================
//function : GetPosition
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::GetPosition (Standard_Real& theX, Standard_Real& theY) const
{
    theX = myXPos;
    theY = myYPos;
}

//=======================================================================
//function : SetPosition
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetPosition (const Standard_Integer theX, const Standard_Integer theY)
{
    if (myXPos == theX && myYPos == theY)
        return;

    myXPos = theX;
    myYPos = theY;
}

//=======================================================================
//function : SetXPosition
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetXPosition (const Standard_Integer theX)
{
    SetPosition (theX, GetYPosition());
}

//=======================================================================
//function : SetYPosition
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetYPosition (const Standard_Integer theY)
{
    SetPosition (GetXPosition(), theY);
}

//=======================================================================
//function : GetSize
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::GetSize (Standard_Integer& theBreadth, Standard_Integer& theHeight) const
{
    theBreadth = myBreadth;
    theHeight = myHeight;
}

//=======================================================================
//function : SetSize
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetSize (const Standard_Integer theBreadth, const Standard_Integer theHeight)
{
    if (myBreadth == theBreadth && myHeight == theHeight)
        return;

    myBreadth = theBreadth;
    myHeight = theHeight;
}

//=======================================================================
//function : SetBreadth
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetBreadth (const Standard_Integer theWidth)
{
    SetSize (theWidth, GetHeight());
}

//=======================================================================
//function : SetHeight
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SetHeight (const Standard_Integer theHeight)
{
    SetSize (GetBreadth(), theHeight);
}

//=======================================================================
//function : SizeHint
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::SizeHint (Standard_Integer& theWidth, Standard_Integer& theHeight) const
{
    Standard_Integer aNum = GetNumberOfIntervals();

    Standard_Integer aSpacer = 5;
    Standard_Integer aTextWidth = 0;
    Standard_Integer aTextHeight = TextHeight ("");
    Standard_Integer aColorWidth = 20;

    if (GetLabelPosition() != Aspect_TOCSP_NONE)
    {
        for (Standard_Integer idx = (IsLabelAtBorder() ? 0 : 1); idx <= aNum; idx++)
            aTextWidth = Max (aTextWidth, TextWidth (GetLabel (idx)));
    }

    Standard_Integer aScaleWidth = 0;
    Standard_Integer aScaleHeight = 0;

    aScaleWidth = aColorWidth + aTextWidth + ( aTextWidth ? 3 : 2 ) * aSpacer;
    aScaleHeight = (Standard_Integer)( 1.5 * ( aNum + (IsLabelAtBorder() ? 2 : 1) ) * aTextHeight );

    Standard_Integer aTitleWidth = 0;
    Standard_Integer aTitleHeight = 0;
    if (GetTitle().Length())
    {
        aTitleHeight = TextHeight (GetTitle()) + aSpacer;
        aTitleWidth =  TextWidth (GetTitle()) + 10;
    }

    theWidth = Max (aTitleWidth, aScaleWidth);
    theHeight = aScaleHeight + aTitleHeight;
}

//=======================================================================
//function : Format
//purpose  :
//=======================================================================
TCollection_AsciiString AIS_ColorScaleExtended::Format() const
{
    return GetFormat();
}

//=======================================================================
//function : GetIntervalValue
//purpose  :
//=======================================================================
Standard_Real AIS_ColorScaleExtended::GetIntervalValue (const Standard_Integer theIndex) const
{
    if (GetNumberOfIntervals() <= 0)
        return 0.;

    if (IsLogarithmic())
    {
        Standard_Real aMin = myMin > 0 ? myMin : 1.0;
        Standard_Real aDivisor = std::pow (myMax/aMin, 1.0/myInterval);
        return aMin*std::pow (aDivisor,theIndex);
    }

    Standard_Real aNum = 0;
    if (GetNumberOfIntervals() > 0)
        aNum = GetMin() + theIndex * ( Abs (GetMax() - GetMin()) / GetNumberOfIntervals() );
    return aNum;
}

//=======================================================================
//function : HueFromValue
//purpose  :
//=======================================================================
Standard_Integer AIS_ColorScaleExtended::HueFromValue (const Standard_Integer theValue,
                                                       const Standard_Integer theMin, const Standard_Integer theMax)
{
    Standard_Integer aMinLimit (0), aMaxLimit (230);

    Standard_Integer aHue = aMaxLimit;
    if (theMin != theMax)
        aHue = (Standard_Integer)( aMaxLimit - ( aMaxLimit - aMinLimit ) * ( theValue - theMin ) / ( theMax - theMin ) );

    aHue = Min (Max (aMinLimit, aHue), aMaxLimit);

    return aHue;
}

//=======================================================================
//function : FindColor
//purpose  :
//=======================================================================
Standard_Boolean AIS_ColorScaleExtended::FindColor (const Standard_Real theValue,
                                                    Quantity_Color& theColor) const
{
    if (theValue < myMin || theValue > myMax || myMax < myMin)
    {
        theColor = Quantity_Color();
        return Standard_False;
    }

    if (GetColorType() == Aspect_TOCSD_USER)
    {
        Standard_Integer anIndex = 0;
        if (Abs (myMax - myMin) > Precision::Approximation())
        {
            anIndex = (theValue - myMin < Precision::Confusion())
                    ? 1
                    : Standard_Integer (Ceiling (( theValue - myMin ) / ( (myMax - myMin) / myInterval)));
        }

        if (anIndex <= 0 || anIndex > myColors.Length())
        {
            theColor = Quantity_Color();
            return Standard_False;
        }

        theColor = myColors.Value (anIndex);
        return Standard_True;
    }

    return FindColor (theValue, myMin, myMax, myInterval, theColor);
}

//=======================================================================
//function : FindColor
//purpose  :
//=======================================================================
Standard_Boolean AIS_ColorScaleExtended::FindColor (const Standard_Real theValue,
                                                    const Standard_Real theMin,
                                                    const Standard_Real theMax,
                                                    const Standard_Integer theColorsCount,
                                                    Quantity_Color& theColor)
{
    if (theValue < theMin || theValue > theMax || theMax < theMin)
        return Standard_False;

    else
    {
        Standard_Real anIntervNumber = 0;
        if(Abs (theMax-theMin) > Precision::Approximation())
            anIntervNumber = Floor (Standard_Real (theColorsCount) * ( theValue - theMin ) / ( theMax - theMin ));

        Standard_Integer anInterv = Standard_Integer (anIntervNumber);

        theColor = Quantity_Color (HueFromValue (anInterv, 0, theColorsCount - 1), 1.0, 1.0, Quantity_TOC_HLS);

        return Standard_True;
    }
}

//=======================================================================
//function : Compute
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::Compute(const occHandle(PrsMgr_PresentationManager3d)& /*thePresentationManager*/,
                                     const occHandle(Prs3d_Presentation)& thePresentation,
                                     const Standard_Integer /*theMode*/)
{
    occHandle(V3d_Viewer) aViewer= GetContext()->CurrentViewer();
    aViewer->InitActiveViews();
    Standard_Integer aNum = GetNumberOfIntervals();
    Aspect_TypeOfColorScalePosition aLabPos = GetLabelPosition();

    Standard_Integer aSpacer = 10;
    Standard_Integer aTextWidth = 1;
    Standard_Integer aTextHeight = myTextHeight;
    Standard_Boolean toDrawLabel = GetLabelPosition() != Aspect_TOCSP_NONE;
    TCollection_ExtendedString aTitle = GetTitle();
    Standard_Integer aTitleHeight = aSpacer;
    //Quantity_Color aFgColor (hasOwnColor ? myOwnColor : Quantity_NOC_BLACK);
    Quantity_Color aFgColor = Quantity_NOC_BLACK;

    //! -----------
    //! Draw title
    //! -----------
    if (aTitle.Length())
    {
        aTitleHeight += myTextHeight + aSpacer;
        // cese
        drawText (thePresentation, aTitle, (Standard_Integer)myXPos + aSpacer, myHeight - ((Standard_Integer)myYPos - 1.0 * aSpacer + aTitleHeight), aFgColor);
    }

    Standard_Boolean toReverse = IsReversed();

    Aspect_SequenceOfColor aColors;
    for (Standard_Integer i = 1; i <= aNum; i++)
    {
        if (toReverse)
        {
            aColors.Prepend (GetIntervalColor (i));
        }
        else
        {
            aColors.Append (GetIntervalColor (i));
        }
    }

    TColStd_SequenceOfExtendedString aLabels;
    Standard_Integer aLabCount = IsLabelAtBorder() ? aNum + 1 : aNum;
    for (Standard_Integer i = 1; i <= aLabCount; i++)
    {
        if (toReverse)
        {
            aLabels.Prepend (GetLabel (i));
        }
        else
        {
            aLabels.Append (GetLabel (i));
        }
    }

    if (toDrawLabel)
        for (Standard_Integer i = 1; i <= aLabels.Length(); i++)
            aTextWidth = Max (aTextWidth, TextWidth (aLabels.Value (i)));

    Standard_Integer aSpc = ( myHeight - ( ( Min (aLabCount, 2) + Abs (aLabCount - aNum - 1) ) * aTextHeight ) - aTitleHeight );
    Standard_Real aVal = aSpc != 0 ? 1.0 * ( aLabCount - Min (aLabCount, 0) ) * aTextHeight / aSpc : 0;
    Standard_Real anIPart;
    Standard_Real anFPart = modf (aVal, &anIPart);
    Standard_Integer aFilter = (Standard_Integer)anIPart + ( anFPart != 0 ? 1 : 0 );

    Standard_Real aStep = 1.0 * ( myHeight - (aLabCount - aNum + Abs (aLabCount - aNum - 1)) * aTextHeight - aTitleHeight ) / aNum;

    Standard_Integer anAscent = 0;

    //! GB: increase the maximum width (breadth...) of the colorbox
    Standard_Integer aColorBreadth = 35;

    // OCC
    //Standard_Integer aColorBreadth = Max (5, Min (20, myBreadth - aTextWidth - 3 * aSpacer));

    if (aLabPos == Aspect_TOCSP_CENTER || !toDrawLabel)
        aColorBreadth += aTextWidth;

    //! ------------
    //! Draw colors
    //! ------------
    Standard_Integer aX = (Standard_Integer)myXPos + aSpacer;
    if (aLabPos == Aspect_TOCSP_LEFT) aX += aTextWidth + ( aTextWidth ? 1 : 0 ) * aSpacer;

    Standard_Real anOffset = 1.0 * aTextHeight / 2 * ( aLabCount - aNum + Abs (aLabCount - aNum - 1) );
    anOffset += 2*aSpacer;
    Handle (Graphic3d_Group) aGroup = Prs3d_Root::CurrentGroup (thePresentation);
    Handle (Graphic3d_ArrayOfTriangles) aPrim;
    Standard_Integer anEdgesPerColor = 6;
    Standard_Integer aVerticiesPerColor = 4;
    aPrim = new Graphic3d_ArrayOfTriangles (aColors.Length()*aVerticiesPerColor, aColors.Length()*anEdgesPerColor, 0, 1);
    Standard_Integer aVertIndex = 1;
    for (Standard_Integer i = 1; i <= aColors.Length() && aStep > 0; i++)
    {
        Standard_Integer aY = (Standard_Integer)( myYPos + ( i - 1 )* aStep + anOffset );
        Standard_Integer aColorHeight = (Standard_Integer)( myYPos + ( i ) * aStep + anOffset ) - aY;
        aPrim->AddVertex (gp_Pnt (aX, aY, 0.0), aColors.Value( i ));
        aPrim->AddVertex (gp_Pnt (aX+aColorBreadth, aY, 0.0), aColors.Value( i ));
        aPrim->AddVertex (gp_Pnt (aX, aY+aColorHeight, 0.0), aColors.Value( i ));
        aPrim->AddVertex (gp_Pnt (aX+aColorBreadth, aY+aColorHeight, 0.0), aColors.Value( i ));
        aPrim->AddEdge(aVertIndex);
        aPrim->AddEdge(aVertIndex+1);
        aPrim->AddEdge(aVertIndex+2);
        aPrim->AddEdge(aVertIndex+1);
        aPrim->AddEdge(aVertIndex+2);
        aPrim->AddEdge(aVertIndex+3);
        aVertIndex += 4;

        //! ---------------------------------------------------------------
        //! experimental - draw a black thin frame aroung each colored box
        //! ---------------------------------------------------------------
        drawFrame(thePresentation,aX,aY,aColorBreadth,aColorHeight,Quantity_NOC_BLACK);
    }
    aGroup->AddPrimitiveArray (aPrim);

    if (aStep > 0)
        drawFrame (thePresentation, aX - 1, (Standard_Integer)(myYPos + anOffset - 1), aColorBreadth + 2, (Standard_Integer)(aColors.Length() * aStep + 2), aFgColor);

    //! ------------
    //! Draw Labels
    //! ------------
    anOffset = 1.0 * Abs (aLabCount - aNum - 1) * ( aStep - aTextHeight ) / 2 + 1.0 * Abs (aLabCount - aNum - 1) * aTextHeight / 2;
    anOffset += 2*aSpacer;
    if (toDrawLabel && aLabels.Length() && aFilter > 0)
    {
        Standard_Integer i1 = 0;
        Standard_Integer i2 = aLabCount - 1;
        Standard_Integer aLast1( i1 ), aLast2( i2 );
        aX = (Standard_Integer)myXPos + aSpacer;
        switch (aLabPos)
        {
        case Aspect_TOCSP_NONE:
        case Aspect_TOCSP_LEFT:
            break;
        case Aspect_TOCSP_CENTER:
            aX += ( aColorBreadth - aTextWidth ) / 2;
            break;
        case Aspect_TOCSP_RIGHT:
            aX += aColorBreadth + aSpacer;
            break;
        }
        while (i2 - i1 >= aFilter || ( i2 == 0 && i1 == 0 ))
        {
            Standard_Integer aPos1 = i1;
            Standard_Integer aPos2 = aLabCount - 1 - i2;
            if (aFilter && !( aPos1 % aFilter ))
            {
                drawText (thePresentation, aLabels.Value (i1 + 1), aX, (Standard_Integer)( myYPos + i1 * aStep + anAscent + anOffset ), aFgColor);
                aLast1 = i1;
            }
            if (aFilter && !( aPos2 % aFilter ))
            {
                drawText (thePresentation, aLabels.Value (i2 + 1), aX, (Standard_Integer)( myYPos + i2 * aStep + anAscent + anOffset ), aFgColor);
                aLast2 = i2;
            }
            i1++;
            i2--;
        }
        Standard_Integer aPos = i1;
        Standard_Integer i0 = -1;
        while (aPos <= i2 && i0 == -1)
        {
            if (aFilter && !( aPos % aFilter ) && Abs (aPos - aLast1) >= aFilter && Abs (aPos - aLast2) >= aFilter)
                i0 = aPos;
            aPos++;
        }

        if (i0 != -1)
            drawText (thePresentation, aLabels.Value (i0 + 1), aX, (Standard_Integer)( myYPos + i0 * aStep + anAscent + anOffset ), aFgColor);
    }
}

//=======================================================================
//function : drawFrame
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::drawFrame (const occHandle(Prs3d_Presentation)& thePresentation,
                                        const Standard_Integer theX, const Standard_Integer theY,
                                        const Standard_Integer theWidth, const Standard_Integer theHeight,
                                        const Quantity_Color& theColor)
{
    Handle (Graphic3d_Group) aGroup = Prs3d_Root::CurrentGroup (thePresentation);
    occHandle(Graphic3d_ArrayOfPolylines) aPrim = new Graphic3d_ArrayOfPolylines(5);
    aPrim->AddVertex (theX,theY,0.0);
    aPrim->AddVertex (theX+theWidth,theY,0.0);
    aPrim->AddVertex (theX+theWidth,theY+theHeight,0.0);
    aPrim->AddVertex (theX,theY+theHeight,0.0);
    aPrim->AddVertex (theX,theY,0.0);

    //! thickness changed to 1.5 from 1.0
    occHandle(Prs3d_LineAspect) anAspect =
            new Prs3d_LineAspect (theColor, Aspect_TOL_SOLID, 1.5);
    anAspect->SetColor (theColor);
    aGroup->SetPrimitivesAspect (anAspect->Aspect());
    aGroup->AddPrimitiveArray (aPrim);
}

//=======================================================================
//function : drawText
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::drawText (const occHandle(Prs3d_Presentation)& thePresentation,
                                       const TCollection_ExtendedString& theText,
                                       const Standard_Integer theX, const Standard_Integer theY,
                                       const Quantity_Color& theColor)
{
    if (!myDrawer->HasOwnTextAspect())
    {
        myDrawer->SetTextAspect (new Prs3d_TextAspect());
        *myDrawer->TextAspect()->Aspect() = *myDrawer->Link()->TextAspect()->Aspect();
    }
    occHandle(Prs3d_TextAspect) anAspect = myDrawer->TextAspect();
    anAspect->SetColor (theColor);
    anAspect->SetHeight (myTextHeight);
    anAspect->SetHorizontalJustification (Graphic3d_HTA_LEFT);
    anAspect->SetVerticalJustification (Graphic3d_VTA_BOTTOM);
    anAspect->Aspect()->SetTextZoomable (Standard_True);
    anAspect->Aspect()->SetTextAngle (0.0);
    anAspect->Aspect()->SetTextFontAspect (Font_FA_Regular);
    anAspect->Aspect()->SetFont("Helvetica");
    Prs3d_Text::Draw (Prs3d_Root::CurrentGroup (thePresentation), anAspect, theText, gp_Ax2 (gp_Pnt (theX, theY, 0.0), gp::DZ()), Standard_False);
}

//=======================================================================
//function : TextWidth
//purpose  :
//=======================================================================
Standard_Integer AIS_ColorScaleExtended::TextWidth (const TCollection_ExtendedString& theText) const
{
    Standard_Integer aWidth, anAscent, aDescent;
    TextSize (theText, GetTextHeight(), aWidth, anAscent, aDescent);
    return aWidth;
}

//=======================================================================
//function : TextHeight
//purpose  :
//=======================================================================
Standard_Integer AIS_ColorScaleExtended::TextHeight (const TCollection_ExtendedString& theText) const
{
    Standard_Integer aWidth, anAscent, aDescent;
    TextSize (theText, GetTextHeight(), aWidth, anAscent, aDescent);
    return anAscent+aDescent;
}

//=======================================================================
//function : TextSize
//purpose  :
//=======================================================================
void AIS_ColorScaleExtended::TextSize (const TCollection_ExtendedString& theText, const Standard_Integer theHeight, Standard_Integer& theWidth, Standard_Integer& theAscent, Standard_Integer& theDescent) const
{
    const occHandle(Graphic3d_CView)& aView = GetContext()->CurrentViewer()->ActiveViewIterator().Value()->View();
    Standard_ShortReal aWidth(10.0), anAscent(1.0), aDescent(1.0);
    TCollection_AsciiString aText (theText.ToExtString(), '?');
    GetContext()->CurrentViewer()->Driver()->TextSize (aView, aText.ToCString(), (Standard_ShortReal)theHeight, aWidth, anAscent, aDescent);
    theWidth = (Standard_Integer)aWidth;
    theAscent = (Standard_Integer)anAscent;
    theDescent = (Standard_Integer)aDescent;
}
