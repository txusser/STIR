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

#include "stir/VoxelsOnCartesianGrid.h"

#include "stir/SSRB.h"
#include "stir/DataProcessor.h"
#include "stir/PostFiltering.h"
#include "stir/scatter/CreateTailMaskFromACFs.h"
#include "stir/scatter/SingleScatterSimulation.h"

#include "stir/zoom.h"
#include "stir/IO/write_to_file.h"
#include "stir/IO/read_from_file.h"
#include "stir/ArrayFunction.h"
#include "stir/stir_math.h"
#include "stir/NumericInfo.h"

#include "stir/SegmentByView.h"

#include "stir/recon_buildblock/IterativeReconstruction.h"

#include <sstream>

// The calculation of the attenuation coefficients
#include "stir/recon_buildblock/ForwardProjectorByBinUsingRayTracing.h"
#include "stir/recon_buildblock/ForwardProjectorByBinUsingProjMatrixByBin.h"
#include "stir/recon_buildblock/ProjMatrixByBinUsingRayTracing.h"

#include "stir/recon_buildblock/BinNormalisationFromAttenuationImage.h"
#include "stir/recon_buildblock/BinNormalisationFromProjData.h"
#include "stir/recon_buildblock/TrivialBinNormalisation.h"

START_NAMESPACE_STIR

void
ScatterEstimationByBin::
set_defaults()
{

    // All recomputes default true
    this->recompute_initial_activity_image = true;
    this->recompute_atten_coeff = true;
    this->recompute_mask_atten_image = true;
    this->recompute_mask_projdata = true;

    this->initial_activity_image_filename = "";
    this->atten_image_filename = "";

    this->output_projdata_filename = "";

    this->num_scatter_iterations = 5;

    this->iterative_method = true;

}

void
ScatterEstimationByBin::
initialise_keymap()
{
    this->parser.add_start_key("Scatter Estimation Parameters");
    this->parser.add_stop_key("end Scatter Estimation Parameters");
    // N.E. 13/07/16: I don't like "input file" for the input data.
    // I try to keep consistency with the reconstruction
    // params.
    this->parser.add_key("input file",
                         &this->input_projdata_filename);
    this->parser.add_key("attenuation image filename",
                         &this->atten_image_filename);

    // MASK parameters
    this->parser.add_key("mask attenuation image filename",
                         &this->mask_atten_image_filename);
    this->parser.add_key("mask postfilter filename",
                         &this->mask_postfilter_filename);
    //-> for attenuation
    this->parser.add_key("recompute mask attenuation image",
                         &this->recompute_mask_atten_image);
    this->parser.add_key("mask attenuation max threshold ",
                         &this->mask_attenuation_image.max_threshold);
    this->parser.add_key("mask attenuation add scalar",
                         &this->mask_attenuation_image.add_scalar);
    this->parser.add_key("mask attenuation min threshold",
                         &this->mask_attenuation_image.min_threshold);
    this->parser.add_key("mask attenuation times scalar",
                         &this->mask_attenuation_image.times_scalar);
    //-> for activity
    this->parser.add_key("mask activity max threshold ",
                         &this->mask_activity_image.max_threshold);
    this->parser.add_key("mask activity add scalar",
                         &this->mask_activity_image.add_scalar);
    this->parser.add_key("mask activity min threshold",
                         &this->mask_activity_image.min_threshold);
    this->parser.add_key("mask activity times scalar",
                         &this->mask_activity_image.times_scalar);
    // END MASK IMAGE
    // MASK PROJDATA
    this->parser.add_key("recompute mask projdata",
                         &this->recompute_mask_projdata);
    this->parser.add_key("mask projdata filename",
                         &this->mask_projdata_filename);
    this->parser.add_key("tail fitting par filename",
                         &this->tail_mask_par_filename);
    // END MASK PROJDATA

    this->parser.add_key("attenuation projdata filename",
                         &this->atten_coeff_filename);
    this->parser.add_key("recompute attenuation coefficients",
                         &this->recompute_atten_coeff);

    this->parser.add_key("background projdata filename",
                         &this->back_projdata_filename);

    this->parser.add_key("normalisation projdata filename",
                         &this->normalisation_projdata_filename);

    this->parser.add_key("recompute initial activity image",
                         &this->recompute_initial_activity_image);
    this->parser.add_key("initial activity image filename",
                         &this->initial_activity_image_filename);

    // ITERATIONS RELATED
    this->parser.add_key("reconstruction template filename",
                         &this->reconstruction_template_par_filename);
    this->parser.add_key("number of scatter iterations",
                         &this->num_scatter_iterations);
    //END ITERATIONS RELATED

    //Scatter simulation
    this->parser.add_key("scatter simulation parameters",
                         &this->scatter_sim_par_filename);

    //    this->reconstruction_template_sptr.reset(new Reconstruction<DiscretisedDensity<3, float > > ( rec_filename ));

    this->parser.add_key("output filename prefix", &this->output_projdata_filename);
}

bool
ScatterEstimationByBin::
post_processing()
{
    shared_ptr<ProjData> input_projdata_sptr =
            ProjData::read_from_file(this->input_projdata_filename);

    //Create the ProjDataInfo to be used in SSRB
    if (input_projdata_sptr->get_max_segment_num() > 0 )
    {
        this->proj_data_info_2d_sptr.reset(
                    dynamic_cast<ProjDataInfoCylindricalNoArcCorr* >
                    (SSRB(*input_projdata_sptr->get_proj_data_info_ptr(),
                          input_projdata_sptr->get_max_segment_num())));

        size_t lastindex = this->input_projdata_filename.find_last_of(".");
        std::string rawname = this->input_projdata_filename.substr(0, lastindex);
        std::string out_filename = rawname + "_2d.hs";
        this->input_projdata_2d_sptr.reset(new ProjDataInterfile(input_projdata_sptr->get_exam_info_sptr(),
                                                                 this->proj_data_info_2d_sptr,
                                                                 out_filename,
                                                                 std::ios::in | std::ios::out | std::ios::trunc));

        SSRB(*this->input_projdata_2d_sptr.get(),
             *input_projdata_sptr.get(),false);
    }
    else
    {
        this->input_projdata_2d_sptr = input_projdata_sptr;
        this->proj_data_info_2d_sptr.reset(
                    input_projdata_sptr->get_proj_data_info_sptr().get());
    }

    if (!this->recompute_initial_activity_image && this->initial_activity_image_filename.size() > 0 )
    {
        this->activity_image_sptr =
                read_from_file<DiscretisedDensity<3,float> >(this->initial_activity_image_filename);
    }

    //
    // Initialise the reconstruction method
    //

    if (this->reconstruction_template_par_filename.size() > 0 )
    {
        KeyParser local_parser;
        local_parser.add_start_key("Reconstruction");
        local_parser.add_stop_key("End Reconstruction");
        local_parser.add_parsing_key("reconstruction method", &this->reconstruction_template_sptr);
        local_parser.parse(this->reconstruction_template_par_filename.c_str());
    }
    else
    {
        error("Please define a reconstruction method.");
    }

    // Load the attenuation image.
    {
        this->atten_image_sptr =
                read_from_file<DiscretisedDensity<3,float> >(this->atten_image_filename);

        info(boost::format("Attenuation image data are supposed to be in units cm^-1\n"
                           "\tReference: water has mu .096 cm^-1\n"
                           "\tMax in attenuation image: %g\n") %
             this->atten_image_sptr->find_max());

        int min_z = this->atten_image_sptr->get_min_index();
        int min_y = this->atten_image_sptr.get()[0][min_z].get_min_index();
        int len_y = this->atten_image_sptr.get()[0][min_z].get_length();
        int len_x = this->atten_image_sptr.get()[0][min_z][min_y].get_length();

        if (len_y != len_x)
            error("The voxels in the x and y dimensions are different. Cannot zoom...  ");

    }

    //
    // ScatterSimulation
    //

    if (this->scatter_sim_par_filename.size() > 0 )
    {
        KeyParser local_parser;
        local_parser.add_start_key("Scatter Simulation");
        local_parser.add_stop_key("End Scatter Simulation");
        local_parser.add_parsing_key("Simulation method", &this->scatter_simulation_sptr);
        local_parser.parse(this->scatter_sim_par_filename.c_str());

        // The image is provided to the simulation.
        // and it will override anything that the ScatterSimulation.par file has done.
        this->scatter_simulation_sptr->set_density_image_and_subsample(this->atten_image_sptr);
        this->scatter_simulation_sptr->set_projdata_and_subsample(this->proj_data_info_2d_sptr);
        this->scatter_simulation_sptr->set_exam_info_sptr(this->input_projdata_2d_sptr->get_exam_info_sptr());
    }
    else
    {
        error("Please define a scatter simulation method.");
    }

    if (this->reconstruction_template_sptr->get_registered_name() == "OSMAPOSL" ||
            this->reconstruction_template_sptr->get_registered_name() == "OSSPS")
    {
        if (set_up_iterative() == false)
            error("Initialisation Failed!");

        this->iterative_method = true;
    }
    else if (this->reconstruction_template_sptr->get_registered_name() == "FBP3D" ||
             this->reconstruction_template_sptr->get_registered_name() == "FBP2D")
    {
        if (set_up_analytic() == false)
            error("Initialisation Failed!");

        this->iterative_method = false;
    }

    //     create output (has to be AFTER set_template_proj_data_info)

    this->output_projdata_sptr.reset(new ProjDataInMemory(this->input_projdata_2d_sptr->get_exam_info_sptr(),
                                                          this->input_projdata_2d_sptr->get_proj_data_info_sptr()->create_shared_clone()));



    return false;
}

ScatterEstimationByBin::
ScatterEstimationByBin()
{
    this->set_defaults();
}

bool
ScatterEstimationByBin::
set_up_iterative()
{
    IterativeReconstruction <DiscretisedDensity<3, float> > * iterative_object =
            dynamic_cast<IterativeReconstruction<DiscretisedDensity<3, float> > *> (this->reconstruction_template_sptr.get());



    this->reconstruction_template_sptr->set_input_data(this->input_projdata_2d_sptr);



    //
    // Initialise the mask image.
    //

    if (this->recompute_mask_atten_image)
    {
        // Applying mask
        // 1. Clone from the original image.
        // 2. Apply to the new clone.
        this->mask_atten_image_sptr.reset(this->atten_image_sptr->clone());
        this->apply_mask_in_place(this->mask_atten_image_sptr,
                                  this->mask_attenuation_image);

        if (this->mask_atten_image_filename.size() > 0 )
            OutputFileFormat<DiscretisedDensity<3,float> >::default_sptr()->
                    write_to_file(this->mask_atten_image_filename, *this->mask_atten_image_sptr.get());
    }
    else if (!this->recompute_mask_atten_image && this->mask_atten_image_filename.size() > 0)
    {
        this->mask_atten_image_sptr =
                read_from_file<DiscretisedDensity<3, float> >(this->mask_atten_image_filename);
    }
    else
        error ("Please set the postfilter parameter filename or set to recompute it.");

    // Forward Project the mask
    if (this->recompute_mask_projdata)
    {

        if (is_null_ptr(this->mask_atten_image_sptr))
            error("You cannot forward project if you have not set the mask image");

        shared_ptr<ForwardProjectorByBin> forw_projector_sptr;
        shared_ptr<ProjMatrixByBin> PM(new  ProjMatrixByBinUsingRayTracing());
        forw_projector_sptr.reset(new ForwardProjectorByBinUsingProjMatrixByBin(PM));
        info(boost::format("\n\nForward projector used for the calculation of\n"
                           "attenuation coefficients: %1%\n") % forw_projector_sptr->parameter_info());

        forw_projector_sptr->set_up(this->input_projdata_2d_sptr->get_proj_data_info_ptr()->create_shared_clone(),
                                    this->mask_atten_image_sptr );

        shared_ptr<ProjData> mask_projdata(new ProjDataInMemory(this->input_projdata_2d_sptr->get_exam_info_sptr(),
                                                                this->input_projdata_2d_sptr->get_proj_data_info_ptr()->create_shared_clone()));

        forw_projector_sptr->forward_project(*mask_projdata, *this->mask_atten_image_sptr);

        //add 1 to be able to use create_tail_mask_from_ACFs (which expects ACFs,
        //so complains if the threshold is too low)

        pow_times_add pow_times_add_object(1.0f,
                                           1.0f,
                                           1.0f,
                                           0.0f,
                                           10000.00f);

        for (int segment_num = mask_projdata->get_min_segment_num();
             segment_num <= mask_projdata->get_max_segment_num();
             ++segment_num)
        {
            SegmentByView<float> segment_by_view =
                    mask_projdata->get_segment_by_view(segment_num);

            in_place_apply_function(segment_by_view,
                                    pow_times_add_object);

            if (!(mask_projdata->set_segment(segment_by_view) == Succeeded::yes))
                warning("Error set_segment %d\n", segment_num);
        }

        if (this->mask_projdata_filename.size() > 0)
            this->mask_projdata_sptr.reset(new ProjDataInterfile(mask_projdata->get_exam_info_sptr(),
                                                                 mask_projdata->get_proj_data_info_ptr()->create_shared_clone(),
                                                                 this->mask_projdata_filename));
        else
            this->mask_projdata_sptr.reset(new ProjDataInMemory(mask_projdata->get_exam_info_sptr(),
                                                                mask_projdata->get_proj_data_info_ptr()->create_shared_clone()));

        CreateTailMaskFromACFs create_tail_mask_from_acfs;
        create_tail_mask_from_acfs.parse(this->tail_mask_par_filename.c_str());

        create_tail_mask_from_acfs.set_input_projdata(mask_projdata);
        create_tail_mask_from_acfs.set_output_projdata(this->mask_projdata_sptr);
        create_tail_mask_from_acfs.process_data();
    }
    else // Load from file
    {
        this->mask_projdata_sptr =
                ProjData::read_from_file(this->mask_projdata_filename);
    }


    //
    // Multiplicative projdata
    //

    shared_ptr<BinNormalisation> _attenuation_correction(new TrivialBinNormalisation());
    shared_ptr<BinNormalisation> _normalisation_coeffs(new TrivialBinNormalisation());

    // Attenuation projdata
    if (this->atten_coeff_filename.size() > 0)
    {
        // Read ProjData and make them 2D

        shared_ptr < ProjData > atten_projdata_sptr;

        atten_projdata_sptr =
                ProjData::read_from_file(this->atten_coeff_filename);

        if( atten_projdata_sptr->get_max_segment_num() > 0)
        {
            size_t lastindex = this->atten_coeff_filename.find_last_of(".");
            std::string rawname = this->atten_coeff_filename.substr(0, lastindex);
            std::string out_filename = rawname + "_2d.hs";

            this->atten_projdata_2d_sptr.reset(new ProjDataInterfile(this->input_projdata_2d_sptr->get_exam_info_sptr(),
                                                                     this->proj_data_info_2d_sptr,
                                                                     out_filename,
                                                                     std::ios::in | std::ios::out | std::ios::trunc));

            SSRB(*this->atten_projdata_2d_sptr.get(),
                 *atten_projdata_sptr.get(), true);
        }
        else
        {
            this->atten_projdata_2d_sptr = atten_projdata_sptr;
        }

        _attenuation_correction.reset(new BinNormalisationFromProjData(this->atten_projdata_2d_sptr));
    }
    else if (this->atten_coeff_filename.size() == 0)
    {
        _attenuation_correction.reset(new BinNormalisationFromAttenuationImage(this->atten_image_sptr));
    }
    //<- End of Attenuation projdata

    // Normalisation ProjData
    if (this->normalisation_projdata_filename.size() > 0 )
    {
        // Read ProjData and make them 2D

        shared_ptr < ProjData > norm_projdata_sptr(ProjData::read_from_file(this->normalisation_projdata_filename));

        if( norm_projdata_sptr->get_max_segment_num() > 0) //If the sinogram is 3D then process it.
        {
            shared_ptr < ProjData> inv_norm_projdata_sptr(new ProjDataInMemory(norm_projdata_sptr->get_exam_info_sptr(),
                                                                               norm_projdata_sptr->get_proj_data_info_ptr()->create_shared_clone()));

            inv_norm_projdata_sptr->fill(*norm_projdata_sptr.get());

            // We need to get norm2d=1/SSRB(1/norm3d))

            pow_times_add pow_times_add_object(0.0f,
                                               1.0f,
                                               -1.0f,
                                               0.0f,
                                               10000.00f);

            for (int segment_num = inv_norm_projdata_sptr->get_min_segment_num();
                 segment_num <= inv_norm_projdata_sptr->get_max_segment_num();
                 ++segment_num)
            {
                SegmentByView<float> segment_by_view =
                        inv_norm_projdata_sptr->get_segment_by_view(segment_num);

                in_place_apply_function(segment_by_view,
                                        pow_times_add_object);

                if (!(inv_norm_projdata_sptr->set_segment(segment_by_view) == Succeeded::yes))
                    warning("Error set_segment %d\n", segment_num);
            }

            size_t lastindex = this->normalisation_projdata_filename.find_last_of(".");
            std::string rawname = this->normalisation_projdata_filename.substr(0, lastindex);
            std::string out_filename = rawname + "_2d.hs";

            this->norm_projdata_2d_sptr.reset(new ProjDataInterfile(this->input_projdata_2d_sptr->get_exam_info_sptr(),
                                                                    this->proj_data_info_2d_sptr,
                                                                    out_filename,
                                                                    std::ios::in | std::ios::out | std::ios::trunc));

            SSRB(*this->norm_projdata_2d_sptr.get(),
                 *inv_norm_projdata_sptr.get(),false);

            for (int segment_num = this->norm_projdata_2d_sptr->get_min_segment_num();
                 segment_num <= this->norm_projdata_2d_sptr->get_max_segment_num();
                 ++segment_num)
            {
                SegmentByView<float> segment_by_view =
                        this->norm_projdata_2d_sptr->get_segment_by_view(segment_num);

                in_place_apply_function(segment_by_view,
                                        pow_times_add_object);

                if (!(this->norm_projdata_2d_sptr->set_segment(segment_by_view) == Succeeded::yes))
                    warning("Error set_segment %d\n", segment_num);
            }


        }
        else // If it is 2D, I assume that it has been processed properly.
        {
            this->norm_projdata_2d_sptr = norm_projdata_sptr;
        }

        _normalisation_coeffs.reset(new BinNormalisationFromProjData(this->norm_projdata_2d_sptr));
    }
    //<- End Normalisation ProjData


    this->_multiplicative_data.reset(new ChainedBinNormalisation(_attenuation_correction,
                                                                 _normalisation_coeffs));

    iterative_object->set_normalisation_sptr(this->_multiplicative_data);

    //
    // Set additive (background) projdata
    //

    if (this->back_projdata_filename.size() > 0)
    {
        shared_ptr<ProjData> back_projdata_sptr =
                ProjData::read_from_file(this->back_projdata_filename);

        if( back_projdata_sptr->get_max_segment_num() > 0)
        {
            size_t lastindex = this->back_projdata_filename.find_last_of(".");
            std::string rawname = this->back_projdata_filename.substr(0, lastindex);
            std::string out_filename = rawname + "_2d.hs";

            this->back_projdata_2d_sptr.reset(new ProjDataInterfile(this->input_projdata_2d_sptr->get_exam_info_sptr(),
                                                                    this->proj_data_info_2d_sptr,
                                                                    out_filename,
                                                                    std::ios::in | std::ios::out | std::ios::trunc));

            SSRB(*this->back_projdata_2d_sptr.get(),
                 *back_projdata_sptr.get(), false);
        }
        else
        {
            this->back_projdata_2d_sptr = back_projdata_sptr;
        }

        iterative_object->set_additive_proj_data_sptr(this->back_projdata_2d_sptr);
    }


    // For testing ... delete later
    //    this->reconstruction_template_sptr->reconstruct();

    //    this->activity_image_sptr = reconstruction_template_sptr->get_target_image();

    //    if (this->initial_activity_image_filename.length() > 0)
    //        OutputFileFormat<DiscretisedDensity < 3, float > >::default_sptr()->
    //                write_to_file(this->initial_activity_image_filename, *activity_image_sptr.get());
    //<- ...

    //    this->
    //  this->set_activity_image(this->activity_image_filename);
    //  this->set_density_image_for_scatter_points(this->density_image_for_scatter_points_filename);
    return true;
}

bool
ScatterEstimationByBin::
set_up_analytic()
{
    //TODO : I have most stuff in tmp.
    return true;
}

void
ScatterEstimationByBin::
apply_mask_in_place(shared_ptr<DiscretisedDensity<3, float> >& arg,
                    const mask_parameters& _this_mask)
{
    // Re-reading every time should not be a problem, as it is
    // just a small txt file.
    PostFiltering filter;
    filter.parser.parse(this->mask_postfilter_filename.c_str());

    // How to use:
    //  pow_times_add(const float add_scalar,
    //  const float mult_scalar, const float power,
    //  const float min_threshold, const float max_threshold)

    // The power is hardcoded because in this context has no
    // use.
    pow_times_add pow_times_add_object(_this_mask.add_scalar,
                                       _this_mask.times_scalar,
                                       1.0,
                                       _this_mask.min_threshold,
                                       _this_mask.max_threshold);

    // 1. filter the image
    filter.filter_ptr->apply(*arg.get());

    // 2. max threshold
    // 3. add scalar
    // 4. min threshold
    // 5. times scalar

    in_place_apply_function(*arg.get(),
                            pow_times_add_object);

}
Succeeded
ScatterEstimationByBin::
process_data()
{
    // Reconstruct an image ... if needed.
    if( is_null_ptr(this->activity_image_sptr) )
    {
        if (iterative_method)
            reconstruct_iterative(0,
                                  this->activity_image_sptr);
        else
            reconstruct_analytic();

        //        // Create the mask image
        //        this->mask_act_image_sptr.reset(this->activity_image_sptr->clone());

        //        // zoom the attenuation to the activity image

        //        apply_mask_in_place(this->mask_act_image_sptr,
        //                            this->mask_activity_image);

        //            std::string output_filename = "./activity_mask";
        if (this->initial_activity_image_filename.size() > 0 )
            OutputFileFormat<DiscretisedDensity < 3, float > >::default_sptr()->
                    write_to_file(this->initial_activity_image_filename, *this->activity_image_sptr.get());
    }

    for (int i_scat_iter = 1;
         i_scat_iter < this->num_scatter_iterations;
         i_scat_iter++)
    {
        // TODO: Mask iteratively...

        // Make a mask out of it


        if (i_scat_iter == 2 ) // Average at 2
        {

        }

        //filter ... optional

        // simulate
        scatter_simulation_sptr->set_output_proj_data(this->output_projdata_sptr);
        scatter_simulation_sptr->set_activity_image_sptr(this->activity_image_sptr);
        scatter_simulation_sptr->process_data();

        // end simulate

        if (iterative_method)
            reconstruct_iterative(i_scat_iter,
                                  this->activity_image_sptr);
        else
        {

            reconstruct_analytic();
            //TODO: threshold ... to cut the negative values
        }
    }

    return Succeeded::yes;
}

Succeeded
ScatterEstimationByBin::
reconstruct_iterative(int _current_iter_num,
                      shared_ptr<DiscretisedDensity<3, float> > & _current_estimate_sptr)
{

    IterativeReconstruction <DiscretisedDensity<3, float> > * iterative_object =
            dynamic_cast<IterativeReconstruction<DiscretisedDensity<3, float> > *> (this->reconstruction_template_sptr.get());

    if (!is_null_ptr(_current_estimate_sptr))
        iterative_object->set_initial_data_ptr(_current_estimate_sptr);

    // Set the previous scatter estimate in the reconstruction process.
    if(_current_iter_num > 0)
        iterative_object->set_additive_proj_data_sptr(this->output_projdata_sptr);

    iterative_object->reconstruct();

    if (is_null_ptr(_current_estimate_sptr))
        _current_estimate_sptr = iterative_object->get_target_image();

    //std::stringstream convert;   // stream used for the conversion

    //    convert << "/Users/nikos/Desktop/Workspace/data/scatters_real/produced_data/initial_image_" <<
    //               _current_iter_num;
    //    std::string output_filename =  convert.str();
    //    OutputFileFormat<DiscretisedDensity < 3, float > >::default_sptr()->
    //            write_to_file(output_filename, *_current_estimate_sptr.get());
}

Succeeded
ScatterEstimationByBin::
reconstruct_analytic()
{
    //TODO
}

/****************** functions to compute scatter **********************/

void
ScatterEstimationByBin::
write_log(const double simulation_time,
          const float total_scatter)
{
    //    std::string log_filename =
    //            this->output_proj_data_filename + ".log";
    //    std::ofstream mystream(log_filename.c_str());

    //    if (!mystream)
    //    {
    //        warning("Cannot open log file '%s'", log_filename.c_str()) ;
    //        return;
    //    }

    //    int axial_bins = 0 ;

    //    for (int segment_num = this->output_proj_data_sptr->get_min_segment_num();
    //         segment_num <= this->output_proj_data_sptr->get_max_segment_num();
    //         ++segment_num)
    //        axial_bins += this->output_proj_data_sptr->get_num_axial_poss(segment_num);

    //    const int total_bins =
    //            this->output_proj_data_sptr->get_num_views() * axial_bins *
    //            this->output_proj_data_sptr->get_num_tangential_poss();
    //    mystream << this->parameter_info()
    //             << "\nTotal simulation time elapsed: "
    //             <<   simulation_time / 60 << "min"
    //               << "\nTotal Scatter Points : " << scatt_points_vector.size()
    //               << "\nTotal Scatter Counts : " << total_scatter
    //               << "\nActivity image SIZE: "
    //               << (*this->activity_image_sptr).size() << " * "
    //               << (*this->activity_image_sptr)[0].size() << " * "  // TODO relies on 0 index
    //               << (*this->activity_image_sptr)[0][0].size()
    //            << "\nAttenuation image SIZE: "
    //            << (*this->atten_image_sptr).size() << " * "
    //            << (*this->atten_image_sptr)[0].size() << " * "
    //            << (*this->atten_image_sptr)[0][0].size()
    //            << "\nTotal bins : " << total_bins << " = "
    //            << this->output_proj_data_sptr->get_num_views()
    //            << " view_bins * "
    //            << axial_bins << " axial_bins * "
    //            << this->output_proj_data_sptr->get_num_tangential_poss()
    //            << " tangential_bins\n";
}

END_NAMESPACE_STIR
