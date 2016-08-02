//
//
/*
  Copyright (C) 2004- 2009, Hammersmith Imanet Ltd
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
  \ingroup utilities
  \ingroup scatter
  \brief Simulates a scatter sinogram

  \author Nikos Efthimiou

  \par Usage:
  \code
  simulate_scatter parfile
  \endcode
  See stir::ScatterSimulation documentation for the format
  of the parameter file.
*/

#include "stir/scatter/ScatterSimulation.h"
#include "stir/Succeeded.h"


int main (int argc, const char *argv[])
{
    using namespace stir;

    shared_ptr <ScatterSimulation> scatter_simulation;

    KeyParser parser;
    parser.add_start_key("Scatter Simulation");
    parser.add_stop_key("End Scatter Simulation");
    parser.add_parsing_key("Simulation method",
                           & scatter_simulation);

    parser.parse(argv[1]);

    return scatter_simulation->process_data() == stir::Succeeded::yes ?
                EXIT_SUCCESS : EXIT_FAILURE;
}
