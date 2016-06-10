#include "stir/ExamInfo.h"
#include "stir/utilities.h"

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>
#ifdef BOOST_NO_STRINGSTREAM
#include <strstream.h>
#else
#include <sstream>
#endif


#ifndef STIR_NO_NAMESPACES
using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::equal;
#endif


START_NAMESPACE_STIR

ExamInfo*
ExamInfo::ask_parameters()
{
    ExamInfo* ei_ptr = new ExamInfo();
    ei_ptr->imaging_modality = ImagingModality::PT;

    ei_ptr->low_energy_thres =
      ask_num("Please set the low energy window threshold (in keV)",0.0f, 1000.0f, -1.0f);

    ei_ptr->up_energy_thres =
      ask_num("Please set the upper energy window threshold (in keV)",0.0f, 1000.0f, -1.0f);

    cout << ei_ptr->parameter_info() <<endl;
    return ei_ptr;
}

string
ExamInfo::parameter_info()  const
{

#ifdef BOOST_NO_STRINGSTREAM
  // dangerous for out-of-range, but 'old-style' ostrstream seems to need this
  char str[30000];
  ostrstream s(str, 30000);
#else
  std::ostringstream s;
#endif

  s << "Lower energy window threshold:"
      << low_energy_thres << "\n";
  s << "Upper energy window threshold:"
      << up_energy_thres << "\n" << endl;

  return s.str();

}


END_NAMESPACE_STIR
