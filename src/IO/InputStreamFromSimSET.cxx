/*
    Copyright (C) 2016, UCL
    Copyright (C) 2018, University of Hull
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    See STIR/LICENSE.txt for details
*/
#include "stir/IO/InputStreamFromSimSET.h"
#include "stir/info.h"
#include "stir/warning.h"
#include "stir/error.h"

extern "C" {
#include <print.header.h>
}


START_NAMESPACE_STIR

const char * const
InputStreamFromSimSET::registered_name =
        "SimSET_History_File";

InputStreamFromSimSET::
InputStreamFromSimSET()
{
    set_defaults();
}


InputStreamFromSimSET::
~InputStreamFromSimSET()
{
    // Try to close history file if opened
    // -Keep in mind that this is not were the actual reading will
    // take place, but in the InputStreamFromSimSET.
    // Here we just need the header.
    if (historyFile != 0)
        fclose(historyFile);
}

std::string
InputStreamFromSimSET::
method_info() const
{
    std::ostringstream s;
    s << this->registered_name;
    return s.str();
}

void
InputStreamFromSimSET::set_defaults()
{
    startFileIndex = 0;
}

Succeeded
InputStreamFromSimSET::
set_up(const std::string & _history_params_filename)
{

    history_filename = _history_params_filename;

    if (set_up_standard_hist_file() == Succeeded::no)
    {
        return set_up_custom_hist_file();
    }

    return Succeeded::yes;
}


Succeeded
InputStreamFromSimSET::set_up_standard_hist_file()
{
    LbUsFourByte		numBluePhotons;				/* Number of blue photons for this decay */
    LbUsFourByte		numPinkPhotons;				/* Number of pink photons for this decay */

    LbUsFourByte		numPhotonsProcessed;		/* Number of photons processed */
    LbUsFourByte		numDecaysProcessed;			/* Number of decays processed */
    PHG_Decay		   	decay;						/* The decay */
    PHG_Decay		   	nextDecay;					/* The decay */
    PHG_DetectedPhoton	detectedPhoton;				/* The detected photon */
    PHG_TrackingPhoton	trackingPhoton;				/* The tracking photon */
    PHG_TrackingPhoton	*bluePhotons = 0;			/* Blue photons for current decay */
    PHG_TrackingPhoton	*pinkPhotons = 0;			/* Pink photon for current decay*/
    double 				angle_norm;					/* for normalizing photon direction */
    Boolean				isOldPhotons1;				/* is this a very old history file--must be
                                                     read using oldReadEvent */
    Boolean				isOldPhotons2;				/* is this a moderately old history file--must be
                                                     read using oldReadEvent */
    Boolean				isOldDecays;				/* is this an old history file--must be
                                                     read using oldReadEvent */
    Boolean				isPHGList;					/* this is PHG history file */
    Boolean				isColList;					/* this is collimator history file */
    Boolean				isDetList;					/* this is detector history file */

    return Succeeded::yes;
}

Succeeded
InputStreamFromSimSET::set_up_custom_hist_file()
{

    return Succeeded::yes;
}

END_NAMESPACE_STIR


//phgrdhstHdrParams.reset(new PhoHFileHdrTy);
//headerHk.reset(new LbHdrHkTy);

//// Let's read the header in the proper structure.
//char *history_cstr = new char[history_filename.length() + 1];
//strcpy(history_cstr, history_filename.c_str());

///* Open data file file */
//if ((historyFile = LbFlFileOpen(history_cstr, "rb")) == nullptr)
//{
//    warning("CListModeDataSimSET: Unable to open history file.");
//    return Succeeded::no;
//}

///* Read in the header and verify it is the right type of file */
//if (PhgHdrGtParams(historyFile, phgrdhstHdrParams.get(), headerHk.get()) == false)
//{
//    warning("CListModeDataSimSET: Unable to read the header. Wrong type or version?.");
//    return Succeeded::no;
//}

///* Display the header */
//info("History file Header: ");
//display(phgrdhstHdrParams.get());


//delete [] history_cstr;
