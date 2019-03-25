//
//
/*!
  \file
  \ingroup IO
  \brief Declaration of class stir::InputStreamFromSimSET

  \author Nikos Efthimiou

*/
/*
    Copyright (C) 2019, University of Hull
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

#ifndef __stir_IO_InputStreamFromSimSET_H__
#define __stir_IO_InputStreamFromSimSET_H__

#include "stir/shared_ptr.h"
#include "stir/Succeeded.h"
#include "stir/listmode/CListRecordSimSET.h"
#include "stir/RegisteredParsingObject.h"

#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhgMath.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "ProdTbl.h"
#include "PhoTrk.h"
#include "SubObj.h"
#include "EmisList.h"
#include "Collimator.h"
#include "Detector.h"
#include "phg.h"
#include "PhgBin.h"
}

/* LOCAL TYPES */
typedef enum  {Null, Decay, Photon} EventTy;

START_NAMESPACE_STIR

//! A helper class to read data from a (presumably binary) stream 
/*! 
*/
class InputStreamFromSimSET
{
public:

    static const char * const registered_name;

    typedef std::vector<LbUsFourByte>::size_type SavedPosition;
    //! Constructor taking a stream
    /*! Data will be assumed to start at the current position reported by seekg().
      If reset() is used, it will go back to this starting position.*/
    InputStreamFromSimSET();

    ~InputStreamFromSimSET();

    //! gives method information
    std::string method_info() const;

    //! Must be called before calling for the first event.
    Succeeded set_up(const std::string & _history_params_filename);

    Succeeded set_up_standard_hist_file();

    Succeeded set_up_custom_hist_file();

    inline
    Succeeded get_next_record(CListRecordSimSET& record) const;

    //! go back to starting position
    inline
    Succeeded reset();

    //! save current "get" position in an internal array
    /*! \return an "index" into the array that allows you to go back.
      \see set_get_position
  */
    inline
    SavedPosition save_get_position();

    //! set current "get" position to previously saved value
    inline
    Succeeded set_get_position(const SavedPosition&);

    //! Function that enables the user to store the saved get_positions
    /*! Together with set_saved_get_positions(), this allows
      reinstating the saved get_positions when
      reopening the same stream.
  */
    inline
    std::vector<LbUsFourByte> get_saved_get_positions() const;
    //! Function that sets the saved get_positions
    /*! Normally, the argument results from a call to
      get_saved_get_positions() on the same stream.
      \warning There is no check if the argument actually makes sense
      for the current stream.
  */
    inline
    void set_saved_get_positions(const std::vector<LbUsFourByte>& );

    inline unsigned long get_total_number_of_events() const;

protected:

    void set_defaults();

private:
    //! The history file we are going to process
    FILE *historyFile;

    std::string history_filename;

    //! Input header
    unique_ptr<PhoHFileHdrTy> phgrdhstHdrParams;
    //! Hook to header
    unique_ptr<LbHdrHkTy> headerHk;
    //! Type of current event.
    //! \warning This might vary in custom history files
    EventTy eventType;
    //! Current file index
    LbUsFourByte curFileIndex;

    LbUsFourByte startFileIndex;

    LbUsFourByte numPhotons;

    std::vector<LbUsFourByte> saved_get_positions;

};

END_NAMESPACE_STIR

#include "stir/IO/InputStreamFromSimSET.inl"

#endif
