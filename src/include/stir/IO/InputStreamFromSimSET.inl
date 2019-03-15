/*
 *  Copyright (C) 2015, 2016 University of Leeds
    Copyright (C) 2016, UCL
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
/*!
  \file
  \ingroup IO
  \brief Implementation of class stir::InputStreamFromROOTFile

  \author Nikos Efthimiou
  \author Harry Tsoumpas
*/

#include "stir/IO/InputStreamFromSimSET.h"

START_NAMESPACE_STIR

//unsigned long int
//InputStreamFromSimSET::
//get_total_number_of_events() const
//{
////    return nentries;
//}

Succeeded
InputStreamFromSimSET::
reset()
{
//    current_position = starting_stream_position;
    return Succeeded::yes;
}

InputStreamFromSimSET::SavedPosition
InputStreamFromSimSET::
save_get_position()
{
//    assert(current_position <= nentries);
//    saved_get_positions.push_back(current_position);
//    return saved_get_positions.size()-1;
}

Succeeded
InputStreamFromSimSET::
set_get_position(const InputStreamFromSimSET::SavedPosition& pos)
{

    return Succeeded::yes;
}

std::vector<std::streampos>
InputStreamFromSimSET::
get_saved_get_positions() const
{
//    return saved_get_positions;
}

void
InputStreamFromSimSET::
set_saved_get_positions(const std::vector<std::streampos> &poss)
{
//    saved_get_positions = poss;
}

Succeeded
InputStreamFromSimSET::
get_next_record(CListRecordSimSET& record) const
{


//    return
//            record.init_from_data(ring1, ring2,
//                                  crystal1, crystal2,
//                                  time1, time2,
//                                  eventID1, eventID2);
}

END_NAMESPACE_STIR
