
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
    //
    // Subsampling stuff.
    //
    int sub_num_dets_per_ring = -1;
    int sub_num_rings = -1;
    float sub_vox_xy = -1.f;
    float sub_vox_z  = -1.f;

    if (this->sub_num_dets_per_ring < 0 && this->sub_vox_xy < 0)
        error ("Please set either number of detectors in subsampled scanner and/or "
               "voxel size in subsampled image");
    if ( this->sub_num_rings < 0 && this->sub_vox_z < 0 )
        error ("Please set either subsampled number of rings and/or subsampled voxel z size");

    // First Load if thats the case and then check
    if ( ( !this->recompute_sub_atten_image && (sub_atten_image_filename.size() > 0)))
    {
        this->set_sub_atten_image_sptr(this->get_image_from_file(this->sub_atten_image_filename));
        sub_vox_xy = this->sub_atten_image_sptr->get_voxel_size()[2];
        sub_vox_z = this->sub_atten_image_sptr->get_voxel_size()[1];
    }
    else if (this->recompute_sub_atten_image)
    {
        sub_vox_xy = (this->sub_vox_xy < 0) ? num_dets_to_vox_size(this->sub_num_dets_per_ring, true)
                                            : this->sub_vox_xy;
        sub_vox_z = (this->sub_vox_z < 0) ? num_dets_to_vox_size(this->sub_num_dets_per_ring, false)
                                          : this->sub_vox_z;
        subsample_image(this->atten_image_sptr,
                        this->sub_atten_image_sptr,
                        sub_vox_xy, sub_vox_z,
                        this->sub_atten_image_filename);
    }

    if ( !this->recompute_sub_projdata && (sub_proj_data_filename.size() > 0))
    {
        this->set_sub_proj_data_info_from_file(this->sub_proj_data_filename);

        sub_num_dets_per_ring = this->sub_proj_data_info_ptr->get_scanner_ptr()->get_num_detectors_per_ring();
        sub_num_rings = this->sub_proj_data_info_ptr->get_scanner_ptr()->get_num_rings();
    }
    else if (this->recompute_sub_projdata)
    {
        sub_num_dets_per_ring = (this->sub_num_dets_per_ring < 0 ) ? vox_size_to_num_dets(this->sub_vox_xy, true)
                                                                   : this->sub_num_dets_per_ring;
        sub_num_rings = (this->sub_num_rings < 0) ? vox_size_to_num_dets(this->sub_vox_z, false)
                                                  : this->sub_num_rings;

        if (sub_num_rings % 2 == 0)
            error("this z voxel size would leed to even number of rings, "
                  "which are supported currently.");

        subsample_projdata_info(this->proj_data_info_ptr,
                                sub_num_dets_per_ring, sub_num_rings,
                                this->sub_proj_data_filename);
    }

    // Check voxel size in comparison to detector size xy

    if (this->recompute_sub_atten_image && this->recompute_sub_projdata)
    {
        if (fabs(this->sub_atten_image_sptr->get_voxel_size()[1] -
                 this->sub_proj_data_info_ptr->get_scanner_ptr()->get_ring_spacing()) > 1.E-2)
            error("DataSymmetriesForBins_PET_CartesianGrid can currently only support z-grid spacing "
                  "equal to the ring spacing of the scanner divided by an integer. Sorry\n");
    }
}


void
ScatterSimulation::
subsample_image(shared_ptr<VoxelsOnCartesianGrid<float> >& _this_image_sptr,
                shared_ptr<VoxelsOnCartesianGrid<float> >& _new_image_sptr,
                float& _sub_vox_xy, float& _sub_vox_z,
                std::__cxx11::string output_filename)
{
    CartesianCoordinate3D <float> _cur_vox = _this_image_sptr->get_voxel_size();
    float zoom_xy = _cur_vox.x() / _sub_vox_xy;
    float zoom_z = _cur_vox.z() / _sub_vox_z;
    int size_xy = _this_image_sptr->get_x_size() * zoom_xy + 0.999;

    if (size_xy % 2 == 0)
        size_xy += 1;

    int size_z = _this_image_sptr->get_z_size() * zoom_z + 0.999;
    float scale_att = zoom_xy * zoom_xy * zoom_z;
    // Just multiply with the scale factor.
    pow_times_add pow_times_add_object(0.0, scale_att, 1.0,
                                       NumericInfo<float>().min_value(),
                                       NumericInfo<float>().max_value());
    // zoom image
    const CartesianCoordinate3D<float>
            zooms(zoom_z, zoom_xy, zoom_xy);
    const CartesianCoordinate3D<float>
            offsets_in_mm = _this_image_sptr->get_origin();
    const CartesianCoordinate3D<int>
            new_sizes(size_z, size_xy, size_xy);
    VoxelsOnCartesianGrid<float> new_image =
            zoom_image(*_this_image_sptr.get(), zooms, offsets_in_mm, new_sizes);
    in_place_apply_function(new_image, pow_times_add_object);
    float dd = new_image.find_max();
    _new_image_sptr.reset(new_image.clone());

    // write it to file
    if (output_filename.size() > 0)
        write_to_file(output_filename, new_image);
}


void
ScatterSimulation::
subsample_projdata_info(ProjDataInfoCylindricalNoArcCorr* _original_projdata_info,
                        int new_num_dets_per_ring, int new_num_rings,
                        std::string output_filename)
{

    Scanner tmpl_scanner = *_original_projdata_info->get_scanner_ptr();
    /*
     * To be used later when scatter correction will be 3D
    const bool is_Advance =
            tmpl_scanner.get_type() == Scanner::Advance ||
            tmpl_scanner.get_type() == Scanner::DiscoveryLS;
    const bool is_DiscoveryST =
            tmpl_scanner.get_type() == Scanner::DiscoveryST;
    const bool is_GE =
            is_Advance || is_DiscoveryST;
    */
    std::string new_name = "subsmpl_" + tmpl_scanner.get_name();
    float rings_ratio = static_cast<float>(tmpl_scanner.get_num_rings()) / new_num_rings;
    float distance_between_rings = tmpl_scanner.get_ring_spacing() * rings_ratio;

    // Effectively this is the subsampled voxel size.
    // We can get this information because the image is calculated first.
    float default_bin_size = this->sub_atten_image_sptr->get_voxel_size()[2];

    shared_ptr<Scanner> new_scanner_ptr(
                new Scanner(Scanner::User_defined_scanner,
                            new_name,
                            new_num_dets_per_ring,
                            new_num_rings,
                            tmpl_scanner.get_max_num_non_arccorrected_bins(),
                            tmpl_scanner.get_default_num_arccorrected_bins(),
                            tmpl_scanner.get_inner_ring_radius(),
                            tmpl_scanner.get_average_depth_of_interaction(),
                            distance_between_rings,
                            default_bin_size,
                            static_cast<float>(0.0), /* should not be hardcoded */
                            tmpl_scanner.get_num_axial_blocks_per_bucket(),
                            tmpl_scanner.get_num_transaxial_blocks_per_bucket(),
                            tmpl_scanner.get_num_axial_crystals_per_block(),
                            tmpl_scanner.get_num_transaxial_crystals_per_block(),
                            tmpl_scanner.get_num_axial_crystals_per_singles_unit(),
                            tmpl_scanner.get_num_transaxial_crystals_per_singles_unit(),
                            tmpl_scanner.get_num_detector_layers(),
                            tmpl_scanner.get_energy_resolution(),
                            tmpl_scanner.get_reference_energy()));

    //    if (new_scanner_ptr->check_consistency() != Succeeded::yes)
    //        error("Initialization of the subsampled scanner failed, sheck the number "
    //              "of detectors per ring and the number of rings\n");

    // get FOV in mm

    float max_tang_pos = this->proj_data_info_ptr->get_max_tangential_pos_num();
    Bin max_tang_bin(0,0,0,max_tang_pos);

    float fov_size_mm = this->proj_data_info_ptr->get_s(max_tang_bin) * 2.f;

    float new_angular_increment = _PI / new_num_dets_per_ring;
    float new_s = tmpl_scanner.get_effective_ring_radius() * sin(new_angular_increment);


    int new_num_view = static_cast<int>((new_num_dets_per_ring / 2) + 0.5);
    int new_num_tang_pos = fov_size_mm / new_s;

    int max_delta = 0;
    int span = 1;

    ProjDataInfo* new_projdata_info =
            span == 0 ? ProjDataInfo::ProjDataInfoGE(new_scanner_ptr, max_delta, new_num_view, new_num_tang_pos, false)
                      : ProjDataInfo::ProjDataInfoCTI(new_scanner_ptr, span, max_delta, new_num_view, new_num_tang_pos, false);

    shared_ptr < ProjDataInfo > pr_sptr(new_projdata_info);
    this->sub_proj_data_info_ptr = dynamic_cast<ProjDataInfoCylindricalNoArcCorr*>(pr_sptr->clone());

    if (output_filename.size() > 0)
        shared_ptr<ProjData> proj_data_ptr(new ProjDataInterfile(this->template_exam_info_sptr,
                                                                 pr_sptr, output_filename));
}

Succeeded
ScatterSimulation::
process_data()
{
    this->initialise_cache_for_scattpoint_det_integrals_over_attenuation();
    this->initialise_cache_for_scattpoint_det_integrals_over_activity();
    ViewSegmentNumbers vs_num;
    /* ////////////////// SCATTER ESTIMATION TIME ////////////////
   */
    CPUTimer bin_timer;
    bin_timer.start();
    // variables to report (remaining) time
    HighResWallClockTimer wall_clock_timer;
    double previous_timer = 0 ;
    int previous_bin_count = 0 ;
    int bin_counter = 0;
    int axial_bins = 0 ;
    wall_clock_timer.start();

    for (vs_num.segment_num() = this->proj_data_info_ptr->get_min_segment_num();
         vs_num.segment_num() <= this->proj_data_info_ptr->get_max_segment_num();
         ++vs_num.segment_num())
        axial_bins += this->proj_data_info_ptr->get_num_axial_poss(vs_num.segment_num());

    const int total_bins =
            this->proj_data_info_ptr->get_num_views() * axial_bins *
            this->proj_data_info_ptr->get_num_tangential_poss();
    /* ////////////////// end SCATTER ESTIMATION TIME ////////////////
   */
    /* Currently, proj_data_info.find_cartesian_coordinates_of_detection() returns
     coordinate in a coordinate system where z=0 in the first ring of the scanner.
     We want to shift this to a coordinate system where z=0 in the middle
     of the scanner.
     We can use get_m() as that uses the 'middle of the scanner' system.
     (sorry)
  */
#ifndef NDEBUG
    {
        CartesianCoordinate3D<float> detector_coord_A, detector_coord_B;
        // check above statement
        this->proj_data_info_ptr->find_cartesian_coordinates_of_detection(
                    detector_coord_A, detector_coord_B, Bin(0, 0, 0, 0));
        assert(detector_coord_A.z() == 0);
        assert(detector_coord_B.z() == 0);
        // check that get_m refers to the middle of the scanner
        const float m_first =
                this->proj_data_info_ptr->get_m(Bin(0, 0, this->proj_data_info_ptr->get_min_axial_pos_num(0), 0));
        const float m_last =
                this->proj_data_info_ptr->get_m(Bin(0, 0, this->proj_data_info_ptr->get_max_axial_pos_num(0), 0));
        assert(fabs(m_last + m_first) < m_last * 10E-4);
    }
#endif
    this->shift_detector_coordinates_to_origin =
            CartesianCoordinate3D<float>(this->proj_data_info_ptr->get_m(Bin(0, 0, 0, 0)), 0, 0);
    float total_scatter = 0 ;

    for (vs_num.segment_num() = this->proj_data_info_ptr->get_min_segment_num();
         vs_num.segment_num() <= this->proj_data_info_ptr->get_max_segment_num();
         ++vs_num.segment_num())
    {
        for (vs_num.view_num() = this->proj_data_info_ptr->get_min_view_num();
             vs_num.view_num() <= this->proj_data_info_ptr->get_max_view_num();
             ++vs_num.view_num())
        {
            total_scatter += this->process_data_for_view_segment_num(vs_num);
            bin_counter +=
                    this->proj_data_info_ptr->get_num_axial_poss(vs_num.segment_num()) *
                    this->proj_data_info_ptr->get_num_tangential_poss();
            /* ////////////////// SCATTER ESTIMATION TIME ////////////////
       */
            {
                wall_clock_timer.stop(); // must be stopped before getting the value
                info(boost::format("%1% bins  Total time elapsed %2% sec "
                                   "\tTime remaining about %3% minutes")
                     % bin_counter
                     % wall_clock_timer.value()
                     % ((wall_clock_timer.value() - previous_timer)
                        * (total_bins - bin_counter) / (bin_counter - previous_bin_count) / 60));
                previous_timer = wall_clock_timer.value() ;
                previous_bin_count = bin_counter ;
                wall_clock_timer.start();
            }
            /* ////////////////// end SCATTER ESTIMATION TIME ////////////////
       */
        }
    }

    bin_timer.stop();
    wall_clock_timer.stop();
    this->write_log(wall_clock_timer.value(), total_scatter);

    if (detection_points_vector.size() != static_cast<unsigned int>(total_detectors))
    {
        warning("Expected num detectors: %d, but found %d\n",
                total_detectors, detection_points_vector.size());
        return Succeeded::no;
    }

    return Succeeded::yes;
}

double
ScatterSimulation::
process_data_for_view_segment_num(const ViewSegmentNumbers& vs_num)
{
    // First construct a vector of all bins that we'll process.
    // The reason for making this list before the actual calculation is that we can then parallelise over all bins
    // without having to think about double loops.
    std::vector<Bin> all_bins;
    {
        Bin bin(vs_num.segment_num(), vs_num.view_num(), 0, 0);

        for (bin.axial_pos_num() = this->proj_data_info_ptr->get_min_axial_pos_num(bin.segment_num());
             bin.axial_pos_num() <= this->proj_data_info_ptr->get_max_axial_pos_num(bin.segment_num());
             ++bin.axial_pos_num())
        {
            for (bin.tangential_pos_num() = this->proj_data_info_ptr->get_min_tangential_pos_num();
                 bin.tangential_pos_num() <= this->proj_data_info_ptr->get_max_tangential_pos_num();
                 ++bin.tangential_pos_num())
            {
                all_bins.push_back(bin);
            }
        }
    }
    // now compute scatter for all bins
    double total_scatter = 0;
    Viewgram<float> viewgram =
            this->output_proj_data_sptr->get_empty_viewgram(vs_num.view_num(), vs_num.segment_num());
#ifdef STIR_OPENMP
#pragma omp parallel for reduction(+:total_scatter) schedule(dynamic)
#endif

    for (int i = 0; i < static_cast<int>(all_bins.size()); ++i)
    {
        const Bin bin = all_bins[i];
        unsigned det_num_A = 0; // initialise to avoid compiler warnings
        unsigned det_num_B = 0;
        this->find_detectors(det_num_A, det_num_B, bin);
        const double scatter_ratio =
                scatter_estimate(det_num_A, det_num_B);
        viewgram[bin.axial_pos_num()][bin.tangential_pos_num()] =
                static_cast<float>(scatter_ratio);
        total_scatter += scatter_ratio;
    } // end loop over bins

    if (this->output_proj_data_sptr->set_viewgram(viewgram) == Succeeded::no)
        error("ScatterEstimationByBin: error writing viewgram");

    return static_cast<double>(viewgram.sum());
}


END_NAMESPACE_STIR
