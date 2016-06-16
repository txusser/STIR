#include "stir/scatter/ScatterEstimationWithFit.h"
#include "stir/ProjDataInterfile.h"
#include "stir/ProjDataInfo.h"
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/is_null_ptr.h"
#include "stir/IO/read_from_file.h"
#include "stir/info.h"
#include "stir/error.h"
#include <fstream>
#include <boost/format.hpp>

START_NAMESPACE_STIR

ScatterEstimationWithFit::
ScatterEstimationWithFit()
{
    this->set_defaults();
}

void
ScatterEstimationWithFit::set_defaults()
{
    base_type::set_defaults();
    //    this->fitted_output_proj_data_filename = "";

    this->transverse_scatter_point_sampling_in_mm = 20;
    this->axial_scatter_point_sampling_in_mm = 20;

    this->num_detectors_per_ring_for_scatter_simulation=48;
    this->num_detectors_rings_for_scatter_simulation=6;

}

void
ScatterEstimationWithFit::initialise_keymap()
{
    base_type::initialise_keymap();

    this->parser.add_key("transverse_scatter_point_sampling_in_mm", &this->transverse_scatter_point_sampling_in_mm);
    this->parser.add_key("axial_scatter_point_sampling_in_mm", &this->axial_scatter_point_sampling_in_mm);

    this->parser.add_key("num_detectors_per_ring_for_scatter_simulation",
                         &this->num_detectors_per_ring_for_scatter_simulation);
    this->parser.add_key("num_detectors_rings_for_scatter_simulation",
                         &this->num_detectors_rings_for_scatter_simulation);

    //    this->parser.add_key("fitted_output_filename_prefix", &this->fitted_output_proj_data_filename);
}


bool
ScatterEstimationWithFit::
pre_processing()
{

    //Apply corrections (randoms etc) to emission data.

    //  : Reconstruct attenuation image (described in the users guide)

    // Use create_tail_mask_from_ACFs to estimate the region that will be used for scaling

    // Construct a "sub-sampled sinogram" template.

    // Construct a "sub-sampled" attenuation image.

    // In the mMR example at this point there is an SSRB step, which can be optional, I think

    // Create sinograms from mask

}


bool
ScatterEstimationWithFit::
post_processing()
{
    //    //base_type::post_processing();

    //    this->set_activity_image(this->activity_image_filename);
    //    this->set_density_image(this->density_image_filename);

    //    info(boost::format("Attenuation image data are supposed to be in units cm^-1\n"
    //                       "\tReference: water has mu .096 cm^-1\n"
    //                       "\tMax in attenuation image: %g\n") %
    //         this->density_image_sptr->find_max());

    //    const VoxelsOnCartesianGrid<float>& activity_image =
    //            dynamic_cast<const VoxelsOnCartesianGrid<float>&>(*activity_image_sptr);

    //    const CartesianCoordinate3D<float> org_voxel_size = activity_image.get_voxel_size();

    //    const CartesianCoordinate3D<float> zooms =
    //            make_coord(this->axial_scatter_point_sampling_in_mm,
    //                       this->transverse_scatter_point_sampling_in_mm,
    //                       this->transverse_scatter_point_sampling_in_mm) /
    //            org_voxel_size;

    //    const CartesianCoordinate3D<float> offsets_in_mm(0.F,0.F,0.F);

    //    BasicCoordinate<3,int> new_sizes

    //            shared_ptr<DiscretisedDensity<3,float> >
    //            scatter_points_image_sptr(new VoxelsOnCartesianGrid<float>(zoom_image(activity_image,zooms, offsets, new_sizes));
    //            this->set_density_image_for_scatter_points_sptr(scatter_points_image_sptr);

    //    // for template_proj_data_filename;
    //    this->num_detectors_per_ring_for_scatter_simulation;
    //    this->num_detectors_rings_for_scatter_simulation;


    //    this->set_template_proj_data_info(this->template_proj_data_filename);

    //    // create output (has to be AFTER set_template_proj_data_info)
    //    this->set_output_proj_data(this->output_proj_data_filename);

    //    return false;

}

Succeeded
ScatterEstimationWithFit::
process_data()
{
//    this->initialise_cache_for_scattpoint_det_integrals_over_attenuation();
//    this->initialise_cache_for_scattpoint_det_integrals_over_activity();

//    ViewSegmentNumbers vs_num;

//    /* ////////////////// SCATTER ESTIMATION TIME ////////////////
//   */
//    CPUTimer bin_timer;
//    int bin_counter = 0;
//    bin_timer.start();
//    int axial_bins = 0 ;
//    for (vs_num.segment_num()=this->proj_data_info_ptr->get_min_segment_num();
//         vs_num.segment_num()<=this->proj_data_info_ptr->get_max_segment_num();
//         ++vs_num.segment_num())
//        axial_bins += this->proj_data_info_ptr->get_num_axial_poss(vs_num.segment_num());
//    const int total_bins =
//            this->proj_data_info_ptr->get_num_views() * axial_bins *
//            this->proj_data_info_ptr->get_num_tangential_poss();

//    /* ////////////////// end SCATTER ESTIMATION TIME ////////////////
//   */

//    /* Currently, proj_data_info.find_cartesian_coordinates_of_detection() returns
//     coordinate in a coordinate system where z=0 in the first ring of the scanner.
//     We want to shift this to a coordinate system where z=0 in the middle
//     of the scanner.
//     We can use get_m() as that uses the 'middle of the scanner' system.
//     (sorry)
//  */
//#ifndef NDEBUG
//    {
//        CartesianCoordinate3D<float> detector_coord_A, detector_coord_B;
//        // check above statement
//        this->proj_data_info_ptr->find_cartesian_coordinates_of_detection(
//                    detector_coord_A,detector_coord_B,Bin(0,0,0,0));
//        assert(detector_coord_A.z()==0);
//        assert(detector_coord_B.z()==0);
//        // check that get_m refers to the middle of the scanner
//        const float m_first =
//                this->proj_data_info_ptr->get_m(Bin(0,0,this->proj_data_info_ptr->get_min_axial_pos_num(0),0));
//        const float m_last =
//                this->proj_data_info_ptr->get_m(Bin(0,0,this->proj_data_info_ptr->get_max_axial_pos_num(0),0));
//        assert(fabs(m_last + m_first)<m_last*10E-4);
//    }
//#endif
//    this->shift_detector_coordinates_to_origin =
//            CartesianCoordinate3D<float>(this->proj_data_info_ptr->get_m(Bin(0,0,0,0)),0, 0);

//    float total_scatter = 0 ;

//    for (vs_num.segment_num()=this->proj_data_info_ptr->get_min_segment_num();
//         vs_num.segment_num()<=this->proj_data_info_ptr->get_max_segment_num();
//         ++vs_num.segment_num())
//    {
//        for (vs_num.view_num()=this->proj_data_info_ptr->get_min_view_num();
//             vs_num.view_num()<=this->proj_data_info_ptr->get_max_view_num();
//             ++vs_num.view_num())
//        {
//            total_scatter += this->process_data_for_view_segment_num(vs_num);
//            bin_counter +=
//                    this->proj_data_info_ptr->get_num_axial_poss(vs_num.segment_num()) *
//                    this->proj_data_info_ptr->get_num_tangential_poss();

//            /* ////////////////// SCATTER ESTIMATION TIME ////////////////
//           */
//            {
//                // TODO remove statics
//                static double previous_timer = 0 ;
//                static int previous_bin_count = 0 ;
//                std::cerr << bin_counter << " bins  Total time elapsed "
//                          << bin_timer.value() << " sec \tTime remaining about "
//                          << (bin_timer.value()-previous_timer)*
//                             (total_bins-bin_counter)/
//                             (bin_counter-previous_bin_count)/60
//                          << " minutes"
//                          << std::endl;
//                previous_timer = bin_timer.value() ;
//                previous_bin_count = bin_counter ;
//            }
//            /* ////////////////// end SCATTER ESTIMATION TIME ////////////////
//           */
//        }
//    }
//    bin_timer.stop();
//    this->write_log(bin_timer.value(), total_scatter);

//    if (detection_points_vector.size() != static_cast<unsigned int>(total_detectors))
//    {
//        warning("Expected num detectors: %d, but found %d\n",
//                total_detectors, detection_points_vector.size());
//        return Succeeded::no;
//    }

//    return Succeeded::yes;
}


END_NAMESPACE_STIR
