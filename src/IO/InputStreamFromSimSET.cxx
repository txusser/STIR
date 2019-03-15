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

START_NAMESPACE_STIR

const char * const
InputStreamFromSimSET::registered_name =
        "SimSET_History_File";

InputStreamFromSimSET::
InputStreamFromSimSET()
{
    set_defaults();
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
}

void
InputStreamFromSimSET::initialise_keymap()
{

}

bool InputStreamFromSimSET::
post_processing()
{

    return false;
}

Succeeded
InputStreamFromSimSET::
set_up(const std::string & header_path)
{
//    if (base_type::set_up(header_path) == Succeeded::no)
//        return  Succeeded::no;

    return Succeeded::yes;
}

END_NAMESPACE_STIR
