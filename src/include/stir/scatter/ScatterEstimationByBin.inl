//
//
/*
    Copyright (C) 2004 - 2012, Hammersmith Imanet Ltd
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
  \ingroup scatter
  \brief Inline implementations of class stir::ScatterEstimationByBin (cross sections)
  
  \author Nikolaos Dikaios
  \author Charalampos Tsoumpas
  \author Kris Thielemans

*/

#include "ScatterEstimationByBin.h"
#include "stir/IO/read_from_file.h"
#include "stir/is_null_ptr.h"
#include "stir/info.h"
#include "stir/error.h"


START_NAMESPACE_STIR

/****************** functions to set images **********************/

shared_ptr<DiscretisedDensity<3,float> >
ScatterEstimationByBin::
get_image_from_file(const std::string& filename)
{
    shared_ptr<DiscretisedDensity<3,float> > _this_image_sptr(
        read_from_file<DiscretisedDensity<3,float> >(filename));

  if (is_null_ptr(_this_image_sptr))
      error("Error reading image %s\n",
            filename.c_str());

  return _this_image_sptr;
}

END_NAMESPACE_STIR
