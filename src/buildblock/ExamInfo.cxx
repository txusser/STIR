#include "stir/ExamInfo.h"

// TODO, Currently all stir::Scanner types are PET.
//

ExamInfo* ExamInfo::ask_parameters()
{
    imaging_modality = ImagingModality::PT;

    ProjDataInfo * pdi_ptr =
      span==0
      ? ProjDataInfoGE(scanner_ptr,max_delta,num_views,num_tangential_poss,arc_corrected)
      : ProjDataInfoCTI(scanner_ptr,span,max_delta,num_views,num_tangential_poss,arc_corrected);

    cout << pdi_ptr->parameter_info() <<endl;

    return pdi_ptr;

}
