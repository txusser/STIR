/*
  Copyright (C) 2004 -  2009 Hammersmith Imanet Ltd
  Copyright (C) 2013 University College London
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
  \brief Implementation of most functions in stir::ScatterEstimationByBin

  \author Charalampos Tsoumpas
  \author Kris Thielemans
*/
#include "stir/scatter/ScatterEstimationByBin.h"
#include "stir/ProjDataInterfile.h"
#include "stir/ProjDataInMemory.h"
#include "stir/ExamInfo.h"
#include "stir/ProjDataInfo.h"
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/Bin.h"
#include "stir/ViewSegmentNumbers.h"
#include "stir/CPUTimer.h"
#include "stir/HighResWallClockTimer.h"
#include "stir/Viewgram.h"
#include "stir/is_null_ptr.h"
#include "stir/IO/read_from_file.h"
#include "stir/info.h"
#include "stir/error.h"
#include <fstream>
#include <boost/format.hpp>

#include "stir/zoom.h"
#include "stir/IO/write_to_file.h"
#include "stir/ArrayFunction.h"
#include "stir/stir_math.h"
#include "stir/NumericInfo.h"

// The calculation of the attenuation coefficients
#include "stir/recon_buildblock/ForwardProjectorByBinUsingRayTracing.h"
#include "stir/recon_buildblock/ForwardProjectorByBinUsingProjMatrixByBin.h"
#include "stir/recon_buildblock/ProjMatrixByBinUsingRayTracing.h"
#include "stir/recon_buildblock/BinNormalisationFromAttenuationImage.h"

START_NAMESPACE_STIR

void
ScatterEstimationByBin::
set_defaults()
{
    this->attenuation_threshold =  0.01 ;
    this->random = true;
    this->use_cache = true;
    this->recompute_initial_activity_image = true;
    this->initial_activity_image_filename = "";
    this->atten_image_filename = "";

    this->template_proj_data_filename = "";
    this->output_proj_data_filename = "";

    this->remove_cache_for_integrals_over_activity();
    this->remove_cache_for_integrals_over_attenuation();
}

void
ScatterEstimationByBin::
initialise_keymap()
{
    this->parser.add_start_key("Scatter Estimation Parameters");
    this->parser.add_stop_key("end Scatter Estimation Parameters");
    this->parser.add_key("projdata template filename",
                         &this->template_proj_data_filename);
    // N.E. 13/07/16: I don't like "input file" for the input data.
    // I try to keep consistency with the reconstruction
    // params.
    this->parser.add_key("input file",
                         &this->input_projdata_filename);
    this->parser.add_key("attenuation image filename",
                         &this->atten_image_filename);
    this->parser.add_key("attenuation projdata filename",
                         &this->atten_coeff_filename);
    this->parser.add_key("background projdata filename",
                         &this->back_projdata_filename);
    this->parser.add_key("recompute attenuation coefficients",
                         &this->recompute_atten_coeff);
    this->parser.add_key("normalization_coefficients_filename",
                         &this->normalization_coeffs_filename);
    this->parser.add_key("recompute initial activity image",
                         &this->recompute_initial_activity_image);
    this->parser.add_key("initial activity image filename",
                         &this->initial_activity_image_filename);

    this->parser.add_key("reconstruction template filename",
                         &this->reconstruction_template_par_filename);
    //    this->reconstruction_template_sptr.reset(new Reconstruction<DiscretisedDensity<3, float > > ( rec_filename ));

    // To this point.
    this->parser.add_key("attenuation_threshold", &this->attenuation_threshold);
    this->parser.add_key("random", &this->random);
    this->parser.add_key("use_cache", &this->use_cache);
    this->parser.add_key("output_filename_prefix", &this->output_proj_data_filename);


}

bool
ScatterEstimationByBin::
post_processing()
{
    // Load the scanner template
    this->set_template_proj_data_info_from_file(this->template_proj_data_filename);

    // create output (has to be AFTER set_template_proj_data_info)
    this->set_proj_data_from_file(this->output_proj_data_filename,
                                  this->output_proj_data_sptr);

    // Load the measured input emission data.
    this->input_projdata_sptr =
            ProjData::read_from_file(this->input_projdata_filename);


    // Load the attenuation image.
    // All relevant initialization have been moved outside
    // the set_image_from_file, since it became more generic.
    shared_ptr<DiscretisedDensity<3, float> > tmp =
            this->get_image_from_file(this->atten_image_filename);
    this->set_atten_image_sptr(tmp);

    // N.E: Load the normalization coefficients if a name has been set.
    if (this->normalization_coeffs_filename.size() > 0)
    {
        this->normalization_factors_sptr =
                ProjData::read_from_file(this->normalization_coeffs_filename);
    }
    else
    {
        this->normalization_factors_sptr.reset(new ProjDataInMemory(this->template_exam_info_sptr,
                                                                    this->output_proj_data_sptr->get_proj_data_info_ptr()->create_shared_clone()));
        this->normalization_factors_sptr->fill(1.f);
    }

    // Recompute the attenuation sinogram if needed.
    if (this->recompute_atten_coeff)
    {
        info(boost::format("Attenuation image data are supposed to be in units cm^-1\n"
                           "\tReference: water has mu .096 cm^-1\n"
                           "\tMax in attenuation image: %g\n") %
             this->atten_image_sptr->find_max());
        this->calculate_atten_coeffs(this->output_proj_data_sptr,
                                     this->atten_image_sptr,
                                     this->atten_coeffs_sptr,
                                     this->atten_coeff_filename);
    }
    else
    {
        this->set_proj_data_from_file(this->atten_coeff_filename,
                                      this->atten_coeffs_sptr);
    }

    // If backround sinogram is set : load it.
    if (this->back_projdata_filename.size() > 0)
        this->back_projdata_sptr =
            ProjData::read_from_file(this->back_projdata_filename);

    //
    // TODO: SSRB
    //


    if (is_null_ptr(this->reconstruction_template_sptr))
    {

        KeyParser local_parser;
        local_parser.add_start_key("Reconstruction");
        local_parser.add_stop_key("End Reconstruction");
        local_parser.add_parsing_key("reconstruction method", &this->reconstruction_template_sptr);
        local_parser.parse(this->reconstruction_template_par_filename.c_str());

        // Probably this cannot work ... I'll test it later.
        //this->reconstruction_template_sptr->initialise(this->reconstruction_template_par_filename);
    }

    //
    // If itinial emission is set then load it,
    // otherwise iterate once, over the initial data.
    //

    if (!this->recompute_initial_activity_image && this->initial_activity_image_filename.size() > 0 )
        this->set_activity_image_sptr( get_image_from_file(this->initial_activity_image_filename) );
    else if (this->recompute_initial_activity_image)     // Initial reconstruction
    {
        this->_iterate();
    }

    //
    // TODO: Masks
    //


    //    this->
    //  this->set_density_image_for_scatter_points(this->density_image_for_scatter_points_filename);

    return false;
}

ScatterEstimationByBin::
ScatterEstimationByBin()
{
    this->set_defaults();
}


/****************** New processing functions; Initially they are going to be marked with '_'
 * later it will be removed.
 */

Succeeded
ScatterEstimationByBin::
_process_data()
{

    //LOOP { NUM }

    // Iterate

    // Simulate


}

Succeeded
ScatterEstimationByBin::
_iterate(int _current_iter_num,
         shared_ptr<ExamData>& _input_data,
         shared_ptr<ExamData>& _mult_data,
         shared_ptr<ExamData>& _add_data,
         shared_ptr<VoxelsOnCartesianGrid<float> >& _current_estimate_sptr)
{
    if (is_null_ptr(this->reconstruction_template_sptr))
        error("There was an error in the initialisation of the reconstruction object.");

    this->reconstruction_template_sptr->set_input_data(this->input_projdata_sptr);

    if (!is_null_ptr(this->back_projdata_sptr))
        this->reconstruction_template_sptr->set_additive_proj_data_sptr(this->back_projdata_sptr);

    // TODO: Set the multiplicative factor for the Analytic reconstruction.
    // Currently implemented only in iterative reconstruction.
    if (!is_null_ptr(this->atten_coeffs_sptr))
        this->reconstruction_template_sptr->set_normalisation_proj_data_sptr(this->atten_coeffs_sptr);

    // Should the attenuations/normalization be
    // different from the subsampled attenuation ( which should include the bed attenuation) ???

    this->reconstruction_template_sptr->reconstruct();

    this->activity_image_sptr.reset( dynamic_cast < VoxelsOnCartesianGrid<float> * > (
                                         reconstruction_template_sptr->get_target_image().get()));

    if (this->initial_activity_image_filename.length() > 0)
        OutputFileFormat<DiscretisedDensity < 3, float > >::default_sptr()->
                write_to_file(this->initial_activity_image_filename, *activity_image_sptr.get());
}

/****************** functions to compute scatter **********************/

Succeeded
ScatterEstimationByBin::
process_data()
{
   // Moved to SS
}

//xxx double
double
ScatterEstimationByBin::
process_data_for_view_segment_num(const ViewSegmentNumbers& vs_num)
{
}


void
ScatterEstimationByBin::
write_log(const double simulation_time,
          const float total_scatter)
{
    std::string log_filename =
            this->output_proj_data_filename + ".log";
    std::ofstream mystream(log_filename.c_str());

    if (!mystream)
    {
        warning("Cannot open log file '%s'", log_filename.c_str()) ;
        return;
    }

    int axial_bins = 0 ;

    for (int segment_num = this->output_proj_data_sptr->get_min_segment_num();
         segment_num <= this->output_proj_data_sptr->get_max_segment_num();
         ++segment_num)
        axial_bins += this->output_proj_data_sptr->get_num_axial_poss(segment_num);

    const int total_bins =
            this->output_proj_data_sptr->get_num_views() * axial_bins *
            this->output_proj_data_sptr->get_num_tangential_poss();
    mystream << this->parameter_info()
             << "\nTotal simulation time elapsed: "
             <<   simulation_time / 60 << "min"
               << "\nTotal Scatter Points : " << scatt_points_vector.size()
               << "\nTotal Scatter Counts : " << total_scatter
               << "\nActivity image SIZE: "
               << (*this->activity_image_sptr).size() << " * "
               << (*this->activity_image_sptr)[0].size() << " * "  // TODO relies on 0 index
               << (*this->activity_image_sptr)[0][0].size()
            << "\nAttenuation image SIZE: "
            << (*this->atten_image_sptr).size() << " * "
            << (*this->atten_image_sptr)[0].size() << " * "
            << (*this->atten_image_sptr)[0][0].size()
            << "\nTotal bins : " << total_bins << " = "
            << this->output_proj_data_sptr->get_num_views()
            << " view_bins * "
            << axial_bins << " axial_bins * "
            << this->output_proj_data_sptr->get_num_tangential_poss()
            << " tangential_bins\n";
}


//!
//! \brief ScatterEstimationByBin::calculate_atten_coeffs
//! \param template_proj_data_ptr
//! \param atten_image_sptr
//! \param out_proj_data_ptr
//! \param output_filename
//! \return
//! \details Rationale: If the normalization file is set then use it to import the normalization
//!  factors. If not initialize it to 1s.
Succeeded
ScatterEstimationByBin::
calculate_atten_coeffs(shared_ptr<ProjData> template_proj_data_ptr,
                       shared_ptr<VoxelsOnCartesianGrid<float> > atten_image_sptr,
                       shared_ptr<ProjData> out_proj_data_ptr,
                       std::string output_filename)
{
    shared_ptr<ForwardProjectorByBin> forw_projector_ptr;
    shared_ptr<ProjMatrixByBin> PM(new  ProjMatrixByBinUsingRayTracing());
    forw_projector_ptr.reset(new ForwardProjectorByBinUsingProjMatrixByBin(PM));
    info(boost::format("\n\nForward projector used for the calculation of\n"
                       "attenuation coefficients: %1%\n") % forw_projector_ptr->parameter_info());

    if (output_filename.size() == 0)
        out_proj_data_ptr.reset(new ProjDataInMemory(template_proj_data_ptr->get_exam_info_sptr(),// TODO this should possibly come from the image, or say it's an ACF File
                                                     template_proj_data_ptr->get_proj_data_info_ptr()->create_shared_clone()));

    else
        out_proj_data_ptr.reset(new ProjDataInterfile(template_proj_data_ptr->get_exam_info_sptr(),// TODO this should possibly come from the image, or say it's an ACF File
                                                      template_proj_data_ptr->get_proj_data_info_ptr()->create_shared_clone(),
                                                      output_filename, std::ios::in | std::ios::out | std::ios::trunc));

    // fill with it with the normalization data which have been initialized
    // either from a file or with 1s.
    out_proj_data_ptr->fill(*this->normalization_factors_sptr);
    // construct a normalisation object that does all the work for us.
    shared_ptr < VoxelsOnCartesianGrid < float > > tmp_atten_image_sptr(atten_image_sptr->clone());
    shared_ptr<BinNormalisation> normalisation_ptr
            (new BinNormalisationFromAttenuationImage(tmp_atten_image_sptr,
                                                      forw_projector_ptr));

    if (normalisation_ptr->set_up(template_proj_data_ptr->get_proj_data_info_ptr()->create_shared_clone())
            != Succeeded::yes)
    {
        error("calculate_attenuation_coefficients: set-up of normalisation failed\n");
        return Succeeded::no;
    }

    // dummy values currently necessary for BinNormalisation, but they will be ignored
    const double start_frame = 0;
    const double end_frame = 0;
    shared_ptr<DataSymmetriesForViewSegmentNumbers> symmetries_sptr(forw_projector_ptr->get_symmetries_used()->clone());
    normalisation_ptr->apply(*out_proj_data_ptr, start_frame, end_frame, symmetries_sptr);
    return Succeeded::yes;
}

END_NAMESPACE_STIR
