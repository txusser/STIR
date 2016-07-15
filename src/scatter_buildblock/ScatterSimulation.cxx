
#include "stir/scatter/ScatterSimulation.h"


START_NAMESPACE_STIR

ScatterSimulation::
ScatterSimulation()
{
    this->set_defaults();
}

ScatterSimulation::
~ScatterSimulation()
{}

Succeeded
ScatterSimulation::
process_data()
{

}

void
ScatterSimulation::set_defaults()
{

    this->sub_atten_image_filename = "";
    this->sub_vox_xy = -1.f;
    this->sub_vox_z = -1.f;
    this->sub_num_dets_per_ring = -1.f;
    this->sub_num_rings = -1.f;

}

void
ScatterSimulation::
ask_parameters()
{

}

void
ScatterSimulation::initialise_keymap()
{
    this->parser.add_key("recompute attenuation coefficients",
                         &this->recompute_atten_coeff);
    this->parser.add_key("normalization_coefficients_filename",
                         &this->normalization_coeffs_filename);
    this->parser.add_key("recompute subsampled attenuation image",
                         &this->recompute_sub_atten_image);
    this->parser.add_key("subsampled attenuation image filename",
                         &this->sub_atten_image_filename);
    this->parser.add_key("subsampled voxel size plane xy (mm/pixel)",
                         &this->sub_vox_xy);
    this->parser.add_key("subsampled voxel size on z axis (mm/pixel)",
                         &this->sub_vox_z);
    this->parser.add_key("recompute subsampled projdata",
                         &this->recompute_sub_projdata);
    this->parser.add_key("subsampled projdata template filename",
                         &this->sub_proj_data_filename);
    this->parser.add_key("number of subsampled detectors",
                         &this->sub_num_dets_per_ring);
    this->parser.add_key("number of subsampled rings",
                         &this->sub_num_rings);
}


bool
ScatterSimulation::
post_processing()
{

}

END_NAMESPACE_STIR
