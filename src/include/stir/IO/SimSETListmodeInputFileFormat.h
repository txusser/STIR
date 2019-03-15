/*
    Copyright (C) 2019 University of Hull
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
#ifndef __stir_IO_SimSETListmodeInputFileFormat_h__
#define __stir_IO_SimSETListmodeInputFileFormat_h__


#include "stir/IO/InputFileFormat.h"
#include "stir/listmode/CListModeDataSimSET.h"
#include "stir/error.h"

#include <SystemDependent.h>
#include <LbTypes.h>
#include <LbError.h>
//#include <LbDebug.h>
//#include <LbEnvironment.h>
//#include <LbFile.h>
#include <LbMemory.h>
//#include <LbInterface.h>
#include <LbParamFile.h>
#include <LbHeader.h>
#include <PhgParams.h>
#include <ColTypes.h>
#include <ColParams.h>
#include <DetTypes.h>
#include <DetParams.h>
#include <CylPos.h>
#include <PhoHFile.h>
#include <PhgHdr.h>


START_NAMESPACE_STIR

//!
//! \brief The SimSETListmodeInputFileFormat class
//! \details Class for being able to read list mode data from the SimSET via the listmode-data registry.
//! At the stage we are not able to use SimSET functions to check the header. That's why we check for two
//! parameters the header version and the header size. The first should be 2.0 (as standing) the second
//! 32768.
//! \author Nikos Efthimiou
//!
class SimSETListmodeInputFileFormat :
        public InputFileFormat<CListModeData >
{
public:
    virtual const std::string
    get_name() const
    {  return "SimSET Histrory file"; }

protected:

    virtual
    bool actual_can_read(const FileSignature& signature,
                         std::istream& input) const
    {
        return this->is_SimSET_signature(signature.get_signature());
    }

    // Large portions of this function come from display_header.c
    bool is_SimSET_signature(const char* const signature) const
    {
        LbUsFourByte cursor = 0;
        LbUsFourByte dist = 0;
        PhoHFileHdrKindTy headerKind;
        float headerVer= 0.0;

        /* Get the header size */
        LbUsFourByte headerSize = 0;
        dist = sizeof(LbUsFourByte);

        for (LbUsFourByte i = 0; i < dist; ++i)
        {
            LbUsFourByte mult = sizeof(LbUsFourByte) - i - 1;
            headerSize |= static_cast<unsigned char>(signature[i+cursor]) << 8*mult;
        }
        cursor += dist;

        if (headerSize != HDR_PHG_HEADER_SIZE_ID)
        {
            warning("SimSETListmodeInputFileFormat: Header is missing 'size' parameter");
            return false;
        }

        /* Get the header kind */
        {
            LbFourByte curID = 0;
            dist = sizeof (PhoHFileHdrKindTy);
            for (LbUsFourByte i = 0; i < dist; ++i)
            {
                LbUsFourByte mult = sizeof(LbUsFourByte) - i - 1;
                curID |= static_cast<LbFourByte>(signature[i+cursor]) << 8*mult;
            }
            cursor += dist;
            memcpy(&headerKind, &curID, dist);

            if (headerKind < PhoHFileEn_PHGOLD || headerKind > PhoHFileEn_DET)
            {
                warning("SimSETLjustistmodeInputFileFormat: Header is missing 'kind' parameter or invalid value");
                return false;
            }
        }

        /* Get the header version */
        {
            LbUsFourByte curID = 0;
            dist = sizeof (float);
            cursor = 32;
            for (LbUsFourByte i = 0; i < dist; ++i)
            {
                LbUsFourByte mult =  i ;

                LbUsFourByte curID2 = static_cast<LbUsFourByte>(signature[i+cursor]);
                curID2 = curID2 << 8*mult;
                curID |= curID2;
            }
            cursor += dist;
            memcpy(&headerVer, &curID, sizeof (float));

            headerVer = static_cast<float>(headerVer);
            if(headerVer != PHG_HDR_HEADER_VERSION)
            {
                warning("SimSETLjustistmodeInputFileFormat: Header is not the same version as with "
                        "the linked SimSET.");
                return false;
            }
        }

        int nikos = 0;
    }



public:    virtual unique_ptr<data_type>
    read_from_file(std::istream& input) const
    {
        error("read_from_file for SimSET listmode data with istream not implemented %s:%s. Sorry",
                __FILE__, __LINE__);
        return unique_ptr<data_type>();
    }

    virtual unique_ptr<data_type>
    read_from_file(const std::string& filename) const
    {
        return unique_ptr<data_type>(new CListModeDataSimSET(filename));
    }

};

END_NAMESPACE_STIR

#endif

//        // Create a new Header hook
//        LbHdrHkTy headerHk;
//        headerHk.fileRef  = nullptr;
//        headerHk.headerSize = signature.size();

//        /* Allocate memory for empty header */
//        if ((headerHk.headerData = LbMmAlloc(static_cast<LbUsFourByte>(headerHk.headerSize))) == nullptr) {
//            return false;
//        }

//        /* Initialize the header to "empty" field value */
//        memset(headerHk.headerData, -1, signature.size());

//        const char* const signature_data = signature.get_signature();
//        for ( int i = 0; i < signature.size(); ++i )
//        {
//           size_t n = std::strlen( &signature_data[i] );

//           std::memcpy( headerHk.headerData, &signature_data[i], n );
//        }
