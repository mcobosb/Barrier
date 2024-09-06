/*!
*
 * \class CDataWriter
 * \brief
 * \details TODO 001 This is a more detailed description of the CDataWriter class
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file data_writer.h
 * \brief Contains DataWWriter definitions
 *
 */
#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include <simulation.h>
#include <string>
using std::string;
#include <map>
using std::map;

#include <netcdf.h>


class CSimulation;

class CDataWriter {

    //! The CSimulation class is a friend of the CDataWriter class
    friend class CSimulation;

private:

    //! Ids for x and t
    int n_XId, n_TId;

    //! Id for the netcdf file
    int m_ncId;

    //! Ids for every variable
    int n_VariableId;

    //! Number of dimensions
    int n_DimensionsIds[2];  // Store x and t dimensions

    //! Length of X
    size_t x_len;

    //! Length of T unlimited
    size_t time_len = NC_UNLIMITED;

    //! Dictionary with the information about the variables
    map<string, map<string, string>> m_mVariableDefinitions;

    //! Dictionary with the Ids of every variable
    map<string, int> m_mVariableIds;


public:


    CDataWriter();
    ~CDataWriter();

    const char * put_time(struct tm * tm, const char * str);

    int nDefineNetCDFFile(const CSimulation* m_pSimulation);
    string strGetVariableMetadata(const string& strVariable, const string& strField);
    int nSetOutputData(CSimulation* m_pSimulation) const;
    int nCloseNetCDFFile(CSimulation* m_pSimulation);
};

#endif // DATA_WRITER_H