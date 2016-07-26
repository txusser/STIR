/*!
  \file
  \ingroup utilities
  \ingroup recon_buildblock
  \brief Demo for Realtime reconstruction initialization

  \author Nikos Efthimiou
  \par main() for any reconstruction which will be created in realtime.
  \code
  recontest parfile
  \endcode

*/


#include "stir/DiscretisedDensity.h"
#include "stir/IO/read_from_file.h"
#include "stir/recon_buildblock/Reconstruction.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include "stir/Succeeded.h"
#include "stir/CPUTimer.h"
#include "stir/HighResWallClockTimer.h"
#include "stir/IO/write_to_file.h"
#include "stir/IO/read_from_file.h"

#include "stir/ProjData.h"
#include "stir/ProjDataInMemory.h"
#include "stir/ProjDataInterfile.h"
#include "stir/recon_buildblock/BinNormalisationFromProjData.h"
#include "stir/recon_buildblock/BinNormalisationFromAttenuationImage.h"
#include "stir/TextWriter.h"
#include "stir/is_null_ptr.h"
#include "stir/recon_buildblock/ChainedBinNormalisation.h"
#include "stir/DiscretisedDensity.h"

#include "stir/recon_buildblock/IterativeReconstruction.h"

using std::cerr;
using std::cout;
using std::endl;

static void print_usage_and_exit()
{
    std::cerr<<"This executable is able to reconstruct some data without calling a specific reconstruction method, from the code.\n";
    std::cerr<<"but specifing the method in the par file with the \"econstruction method\". \n";
    std::cerr<<"\nUsage:\nrecontest reconstuction_test.par\n";
    std::cerr<<"Example parameter file:\n\n"
            <<"Reconstruction :=\n"
           <<"input file := my_prompts.hs\n"
          <<"attenuation input file := my_acfs.hs\n"
         <<"reconstruction method := OSMAPOSL\n"
        <<"OSMAPOSLParameters := \n"
       <<"objective function type:= PoissonLogLikelihoodWithLinearModelForMeanAndProjData\n"
      <<"PoissonLogLikelihoodWithLinearModelForMeanAndProjData Parameters:=\n"
     <<"input file := <input>.hs\n"
    <<"maximum absolute segment number to process := -1\n"
    <<"projector pair type := Matrix\n"
    <<"Projector Pair Using Matrix Parameters :=\n"
    <<"Matrix type := Ray Tracing\n"
    <<"Ray tracing matrix parameters :=\n"
    <<"number of rays in tangential direction to trace for each bin:= 10\n"
    <<"End Ray tracing matrix parameters :=\n"
    <<"End Projector Pair Using Matrix Parameters :=\n"
    <<"recompute sensitivity := 1\n"
    <<"zoom := 1\n"
    <<"end PoissonLogLikelihoodWithLinearModelForMeanAndProjData Parameters:=\n"
    <<"enforce initial positivity condition:= 1 \n"
    <<"number of subsets:= 1\n"
    <<"number of subiterations:= 1 \n"
    <<"save estimates at subiteration intervals:= 1\n"
    <<"output filename prefix := output\n"
    <<"end OSMAPOSLParameters := \n"
    <<"output filename prefix := output_image \n"
    <<"end reconstruction := \n";
    exit(EXIT_FAILURE);
}

/***********************************************************/

int main(int argc, const char *argv[])
{

    using namespace stir;

    TextPrinter tp("std::cout");
    TextWriterHandle twh;
    twh.set_information_channel(&tp);
    twh.set_warning_channel(&tp);
    twh.set_error_channel(&tp);

    if (argc!=2)
        print_usage_and_exit();

    shared_ptr < Reconstruction < DiscretisedDensity < 3, float > > >
            reconstruction_method_sptr;

    std::string data_filename;
    std::string output_filename;
    std::string attenuation_filename;
    std::string attenuation_image_filename;
    std::string normalisation_filename;
    std::string background_filename;
    std::string mult_output_filename = "mult_sino2";
    bool iterative_method = true;

    IterativeReconstruction <DiscretisedDensity<3, float> > * iterative_object;

    shared_ptr<ProjData> recon_data_sptr;
    shared_ptr<ProjData> atten_data_sptr;
    shared_ptr<ProjData> back_data_sptr;
    shared_ptr<ProjData> norm_data_sptr;
    shared_ptr<ProjData> mult_data_sptr;

    shared_ptr<BinNormalisation> bin_norm_sptr;

    shared_ptr<DiscretisedDensity<3,float> > atten_image_sptr;

    twh.print_information("Starting Genererilised reconstruction ... ");

    KeyParser parser;
    parser.add_start_key("Reconstruction");
    parser.add_stop_key("End Reconstruction");
    parser.add_key("input file", &data_filename);
    parser.add_key("attenuation input file", &attenuation_filename);
    parser.add_key("attenuation image file", &attenuation_image_filename);
    parser.add_key("normalisation input file", &normalisation_filename);
    parser.add_key("background input file", &background_filename);
    parser.add_key("output filename prefix", &output_filename);
    parser.add_parsing_key("reconstruction method", &reconstruction_method_sptr);

    twh.print_information("Parsing par file ... ");
    parser.parse(argv[1]);

    HighResWallClockTimer t;
    t.reset();
    t.start();

    if (data_filename.size() > 0)
    {
        recon_data_sptr =
                ProjData::read_from_file(data_filename);

        reconstruction_method_sptr->set_input_data(recon_data_sptr);

        twh.print_information("Input data loaded ... ");
    }
    else
        return false;

    if (reconstruction_method_sptr->get_registered_name() == "OSMAPOSL" ||
            reconstruction_method_sptr->get_registered_name() == "OSSPS")
    {

        iterative_object =
                dynamic_cast<IterativeReconstruction<DiscretisedDensity<3, float> > *> (reconstruction_method_sptr.get());
        iterative_method = true;
    }
    else if (reconstruction_method_sptr->get_registered_name() == "FBP3D" ||
             reconstruction_method_sptr->get_registered_name() == "FBP2D")
    {
       iterative_method = false;
    }

    if (attenuation_filename.size() > 0 )
    {
        atten_data_sptr =
                ProjData::read_from_file(attenuation_filename);
    }

    if (attenuation_image_filename.size() > 0 )
    {
        atten_image_sptr = read_from_file<DiscretisedDensity<3,float> >(attenuation_image_filename);
    }

    if ( normalisation_filename.size() > 0 )
    {
        norm_data_sptr =
                ProjData::read_from_file(normalisation_filename);
    }

    // Example 1: bin normalisation from image file
    //    if (attenuation_image_filename.size() > 0 && iterative_method)
    //    {
    //        bin_norm_sptr.reset(new BinNormalisationFromAttenuationImage(attenuation_image_filename));
    //        reconstruction_method_sptr->set_normalisation_sptr(bin_norm_sptr);
    //        twh.print_information("Bin normalisation from Attenuation image file loaded ...\n ");
    //    }

    // Example 2: bin normalisation from projdata file
    //    if (attenuation_filename.size() > 0 && iterative_method)
    //    {
    //        bin_norm_sptr.reset(new BinNormalisationFromProjData(attenuation_filename) );
    //        reconstruction_method_sptr->set_normalisation_sptr(bin_norm_sptr);
    //        twh.print_information("Bin normalisation from Attenuation ProjData file loaded ...\n ");
    //    }

    // Example 3: bin normalisation from DiscretisedDensity
    //        if (!is_null_ptr(atten_image_sptr) && iterative_method)
    //        {
    //            bin_norm_sptr.reset(new BinNormalisationFromAttenuationImage(atten_image_sptr));
    //            reconstruction_method_sptr->set_normalisation_sptr(bin_norm_sptr);
    //            twh.print_information("Bin normalisation from DiscretisedDensity data loaded ...\n ");
    //        }

    //Example 4: bin normalisation from ProjData.
    //    if (!is_null_ptr(atten_data_sptr) && iterative_method)
    //    {
    //        bin_norm_sptr.reset(new BinNormalisationFromProjData(atten_data_sptr));
    //        reconstruction_method_sptr->set_normalisation_sptr(bin_norm_sptr);
    //        twh.print_information("Bin normalisation from ProjData loaded ...\n ");
    //    }

    // Example 5: Normalisation from ProjData
    //    if (!is_null_ptr(atten_data_sptr) && iterative_method)
    //    {
    //        reconstruction_method_sptr->set_normalisation_proj_data_sptr(atten_data_sptr);
    //        twh.print_information("Bin normalisation from ProjData loaded ...\n ");
    //    }

    //Example 6: Chained BinNormalisation from ProjData
    if (!is_null_ptr(atten_data_sptr) &&
            !is_null_ptr(norm_data_sptr) && iterative_method)
    {
        shared_ptr<BinNormalisationFromProjData> _atten(new BinNormalisationFromProjData(atten_data_sptr));
        shared_ptr<BinNormalisationFromProjData> _norm(new BinNormalisationFromProjData(norm_data_sptr));
        bin_norm_sptr.reset(new ChainedBinNormalisation(_atten, _norm));

        iterative_object->set_normalisation_sptr(bin_norm_sptr);
        twh.print_information("Chained Bin normalisation from ProjData loaded ...\n ");
    }

    if (background_filename.size() > 0 && iterative_method)
    {
        back_data_sptr =
                ProjData::read_from_file(background_filename);

        iterative_object->set_additive_proj_data_sptr(back_data_sptr);

        twh.print_information("Background (scattered / random) data loaded ... ");
    }


    if (reconstruction_method_sptr->reconstruct() == Succeeded::yes)
    {
        t.stop();
        twh.print_information("Reconstruction finished.");

        std::cout << "Total Wall clock time: " << t.value() << " seconds" << std::endl;
    }
    else
    {
        t.stop();
        return Succeeded::no;
    }

    //
    // Save the reconstruction output from this location.
    //

    if (output_filename.length() > 0 )
    {
        twh.print_information("Saving output image.");
        shared_ptr  < DiscretisedDensity < 3, float > > reconstructed_image =
                reconstruction_method_sptr->get_target_image();

        OutputFileFormat<DiscretisedDensity < 3, float > >::default_sptr()->
                write_to_file(output_filename, *reconstructed_image.get());
        twh.print_information("Output image saved.");
    }

    return Succeeded::yes;

    return EXIT_SUCCESS;

}
