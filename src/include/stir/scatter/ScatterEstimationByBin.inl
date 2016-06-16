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
#include "stir/ProjDataInterfile.h"
#include "stir/ProjDataInfo.h"
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/IO/read_from_file.h"
#include "stir/is_null_ptr.h"
#include "stir/info.h"
#include "stir/error.h"
#include <fstream>


START_NAMESPACE_STIR

/****************** functions to set images **********************/
void
ScatterEstimationByBin::
set_activity_image_sptr(const shared_ptr<DiscretisedDensity<3,float> >& new_sptr)
{
  if ( !is_null_ptr(new_sptr) )
    this->activity_image_sptr=new_sptr;

  this->remove_cache_for_integrals_over_activity();
}

void
ScatterEstimationByBin::
set_atten_image_sptr(const shared_ptr<DiscretisedDensity<3,float> >& new_sptr)
{
  if ( !is_null_ptr(new_sptr) )
    this->atten_image_sptr=new_sptr;

  this->remove_cache_for_integrals_over_attenuation();
}

void
ScatterEstimationByBin::
set_image_from_file(const std::string& filename,
                    shared_ptr<DiscretisedDensity<3,float> > & _this_image_sptr)
{
  _this_image_sptr=
    read_from_file<DiscretisedDensity<3,float> >(filename);

  if (is_null_ptr(_this_image_sptr))
      error("Error reading image %s\n",
            filename.c_str());

}

void
ScatterEstimationByBin::
set_sub_atten_image_sptr(const shared_ptr<DiscretisedDensity<3,float> >& new_sptr)
{
  if ( !is_null_ptr(new_sptr) )
    this->sub_atten_image_sptr=new_sptr;

  this->sample_scatter_points();
  this->remove_cache_for_integrals_over_attenuation();
}

/****************** functions to set projection data **********************/


void
ScatterEstimationByBin::
set_template_proj_data_info_from_file(const std::string& filename)
{
  this->template_proj_data_filename = filename;

  shared_ptr<ProjData> template_proj_data_sptr =
    ProjData::read_from_file(this->template_proj_data_filename);

  this->template_exam_info_sptr = template_proj_data_sptr->get_exam_info_sptr();

  this->set_template_proj_data_info_sptr(template_proj_data_sptr->get_proj_data_info_ptr()->create_shared_clone());
}

void
ScatterEstimationByBin::
set_atten_coeffs_from_file(const std::string& filename)
{
  this->atten_coeff_filename = filename;

  this->atten_coeffs_sptr =
    ProjData::read_from_file(this->atten_coeff_filename);
}

void
ScatterEstimationByBin::
set_template_proj_data_info_sptr(const shared_ptr<ProjDataInfo>& new_sptr)
{

  this->proj_data_info_ptr = dynamic_cast<ProjDataInfoCylindricalNoArcCorr const *>(new_sptr->clone());

  if (is_null_ptr(this->proj_data_info_ptr))
    {
      error("ScatterEstimationByBin can only handle non-arccorrected data");
    }

  // find final size of detection_points_vector
  this->total_detectors =
    this->proj_data_info_ptr->get_scanner_ptr()->get_num_rings()*
    this->proj_data_info_ptr->get_scanner_ptr()->get_num_detectors_per_ring ();

  // reserve space to avoid reallocation, but the actual size will grow dynamically
  this->detection_points_vector.reserve(total_detectors);

  // remove any cached values as they'd be incorrect if the sizes changes
  this->remove_cache_for_integrals_over_attenuation();
  this->remove_cache_for_integrals_over_activity();
}


void
ScatterEstimationByBin::
set_proj_data_from_file(const std::string& filename,
                        shared_ptr<ProjData>& _this_projdata)
{
  _this_projdata.reset(new ProjDataInterfile(this->template_exam_info_sptr,
                                                          this->proj_data_info_ptr->create_shared_clone(),
                                                          filename));
}


//
//
// NOT SETS THE REST
//
//


float
ScatterEstimationByBin::
dif_Compton_cross_section(const float cos_theta, float energy)
{ 
  const double Re = 2.818E-13;   // aktina peristrofis electroniou gia to atomo tou H
  const double sin_theta_2= 1-cos_theta*cos_theta ;
  const double P= 1/(1+(energy/511.0)*(1-cos_theta)); 
  return static_cast<float>( (Re*Re/2) * P * (1 - P * sin_theta_2 + P * P));
}

float
ScatterEstimationByBin::
photon_energy_after_Compton_scatter(const float cos_theta, const float energy)
{
  return static_cast<float>(energy/(1+(energy/511.0)*(1-cos_theta)));   // For an arbitrary energy 
}

float
ScatterEstimationByBin::
photon_energy_after_Compton_scatter_511keV(const float cos_theta)
{
  return 511/(2-cos_theta); // for a given energy, energy := 511 keV
}

float
ScatterEstimationByBin::
total_Compton_cross_section(const float energy)
{
  const double a= energy/511.0;
  const double l= log(1.0+2.0*a); 
  const double sigma0= 6.65E-25;   // sigma0=8*pi*a*a/(3*m*m)
  return static_cast<float>( 0.75*sigma0  * ( (1.0+a)/(a*a)*( 2.0*(1.0+a)/(1.0+2.0*a)- l/a ) + l/(2.0*a) - (1.0+3.0*a)/(1.0+2.0*a)/(1.0+2.0*a) ) ); // Klein - Nishina formula = sigma / sigma0
}


float
ScatterEstimationByBin::
total_Compton_cross_section_relative_to_511keV(const float energy)
{
  const double a= energy/511.0;
  static const double prefactor = 9/(-40 + 27*log(3.)); //Klein-Nishina formula for a=1 & devided with 0.75 == (40 - 27*log(3)) / 9

  return //checked this in Mathematica
    static_cast<float>
	(prefactor*
    (((-4 - a*(16 + a*(18 + 2*a)))/square(1 + 2*a) +       
      ((2 + (2 - a)*a)*log(1 + 2*a))/a)/square(a)
     ));
}

void
ScatterEstimationByBin::
set_attenuation_threshold(float _val)
{
    attenuation_threshold = _val;
}

void
ScatterEstimationByBin::
set_random_point(bool _val)
{
    random = _val;
}

void
ScatterEstimationByBin::
set_cache_enabled(bool _val)
{
    use_cache = _val;
}

END_NAMESPACE_STIR
