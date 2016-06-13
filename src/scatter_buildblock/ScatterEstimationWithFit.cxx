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
    this->fitted_output_proj_data_filename = "";

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

    this->parser.add_key("num_detectors_per_ring_for_scatter_simulation", &this->num_detectors_per_ring_for_scatter_simulation);
    this->parser.add_key("num_detectors_rings_for_scatter_simulation", &this->num_detectors_rings_for_scatter_simulation);

    this->parser.add_key("fitted_output_filename_prefix", &this->fitted_output_proj_data_filename);
}

bool
ScatterEstimationWithFit::
post_processing()
{
  //base_type::post_processing();

  this->set_activity_image(this->activity_image_filename);
  this->set_density_image(this->density_image_filename);
  info(boost::format("Attenuation image data are supposed to be in units cm^-1\n"
                     "\tReference: water has mu .096 cm^-1\n"
                     "\tMax in attenuation image: %g\n") %
       this->density_image_sptr->find_max());

  const VoxelsOnCartesianGrid<float>& activity_image =
    dynamic_cast<const VoxelsOnCartesianGrid<float>&>(*activity_image_sptr);

  const CartesianCoordinate3D<float> org_voxel_size = activity_image.get_voxel_size();

  const CartesianCoordinate3D<float> zooms =
    make_coord(this->axial_scatter_point_sampling_in_mm,
               this->transverse_scatter_point_sampling_in_mm,
               this->transverse_scatter_point_sampling_in_mm) /
    org_voxel_size;

  const CartesianCoordinate3D<float> offsets_in_mm(0.F,0.F,0.F);

BasicCoordinate<3,int> new_sizes

  shared_ptr<DiscretisedDensity<3,float> >
    scatter_points_image_sptr(new VoxelsOnCartesianGrid<float>(zoom_image(activity_image,zooms, offsets, new_sizes));
  this->set_density_image_for_scatter_points_sptr(scatter_points_image_sptr);

  // for template_proj_data_filename;
  this->num_detectors_per_ring_for_scatter_simulation;
  this->num_detectors_rings_for_scatter_simulation;


  this->set_template_proj_data_info(this->template_proj_data_filename);

  // create output (has to be AFTER set_template_proj_data_info)
  this->set_output_proj_data(this->output_proj_data_filename);

  return false;

}

/****************** functions to set images **********************/
void
ScatterEstimationWithFit::
set_activity_image_sptr(const shared_ptr<DiscretisedDensity<3,float> >& new_sptr)
{
  this->activity_image_sptr=new_sptr;
  this->remove_cache_for_integrals_over_activity();
}

void
ScatterEstimationWithFit::
set_activity_image(const std::string& filename)
{
  this->activity_image_filename=filename;
  this->activity_image_sptr=
    read_from_file<DiscretisedDensity<3,float> >(filename);
  if (is_null_ptr(this->activity_image_sptr))
    {
      error(boost::format("Error reading activity image %s") %
            this->activity_image_filename);
    }
  this->set_activity_image_sptr(this->activity_image_sptr);
}


void
ScatterEstimationWithFit::
set_density_image_sptr(const shared_ptr<DiscretisedDensity<3,float> >& new_sptr)
{
  this->density_image_sptr=new_sptr;
  this->remove_cache_for_integrals_over_attenuation();
}

void
ScatterEstimationWithFit::
set_density_image(const std::string& filename)
{
  this->density_image_filename=filename;
  this->density_image_sptr=
    read_from_file<DiscretisedDensity<3,float> >(filename);
  if (is_null_ptr(this->density_image_sptr))
    {
      error(boost::format("Error reading density image %s") %
            this->density_image_filename);
    }
  this->set_density_image_sptr(this->density_image_sptr);
}

void
ScatterEstimationWithFit::
set_density_image_for_scatter_points_sptr(const shared_ptr<DiscretisedDensity<3,float> >& new_sptr)
{
  this->density_image_for_scatter_points_sptr=new_sptr;
  this->sample_scatter_points();
  this->remove_cache_for_integrals_over_attenuation();
}

void
ScatterEstimationWithFit::
set_density_image_for_scatter_points(const std::string& filename)
{
  this->density_image_for_scatter_points_filename=filename;
  this->density_image_for_scatter_points_sptr=
    read_from_file<DiscretisedDensity<3,float> >(filename);
  if (is_null_ptr(this->density_image_for_scatter_points_sptr))
    {
      error(boost::format("Error reading density_for_scatter_points image %s") %
            this->density_image_for_scatter_points_filename);
    }
  this->set_density_image_for_scatter_points_sptr(this->density_image_for_scatter_points_sptr);
}


/****************** functions to set projection data **********************/

void
ScatterEstimationWithFit::
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
ScatterEstimationWithFit::
set_template_proj_data_info(const std::string& filename)
{
  this->template_proj_data_filename = filename;
  shared_ptr<ProjData> template_proj_data_sptr =
    ProjData::read_from_file(this->template_proj_data_filename);

  this->set_template_proj_data_info_sptr(template_proj_data_sptr->get_proj_data_info_ptr()->create_shared_clone());
}

/*
void
ScatterEstimationByBin::
set_output_proj_data_sptr(const shared_ptr<ProjData>& new_sptr)
{
  this->output_proj_data_sptr = new_sptr;
}
*/

void
ScatterEstimationWithFit::
set_output_proj_data(const std::string& filename)
{
  this->output_proj_data_filename = filename;
  // TODO get ExamInfo from image
  shared_ptr<ExamInfo> exam_info_sptr(new ExamInfo);
  this->output_proj_data_sptr.reset(new ProjDataInterfile(exam_info_sptr,
                                                          this->proj_data_info_ptr->create_shared_clone(),
                                                          this->output_proj_data_filename));
}

END_NAMESPACE_STIR
