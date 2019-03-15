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

//#include <SystemDependent.h>
//#include <LbTypes.h>
//#include <LbError.h>
//#include <LbDebug.h>
//#include <LbEnvironment.h>
//#include <LbFile.h>
//#include <LbMemory.h>
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


START_NAMESPACE_STIR

//!
//! \brief The SimSETListmodeInputFileFormat class
//! \details Class for being able to read list mode data from the SimSET via the listmode-data registry.
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

    float disHdrVersions = 1.00;
    PhoHFileHdrTy	header;

    virtual
    bool
    actual_can_read(const FileSignature& signature,
                    std::istream& input) const
    {
        return this->is_SimSET_signature(signature.get_signature());
    }

    // Large portions of this function come from display_header.c
    bool is_SimSET_signature(const char* const signature) const
    {

        FILE *historyFile;
        LbUsFourByte		headerSize;	/* Size of header in the file. */
        PhoHFileHdrKindTy	headerKind;	/* The type of this header */
        float headerVersion; /* The version of  this header */


        /* First, read the header size; it is in the first four bytes */
        memcpy(&headerSize, signature, sizeof(LbUsFourByte));
        if (headerSize != sizeof(header)) {
            return 0;
        }
//        sizeof(PhoHFileHdrKindTy);

//        return (
//            standardise_interfile_keyword(keyword) ==
//            standardise_interfile_keyword("SimSET header"));

        int nikos = 0;
    }

public:
    virtual unique_ptr<data_type>
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
