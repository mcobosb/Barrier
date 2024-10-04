/*!
*
 * \class CDataWriter
 * \brief
 * \details Description of CDataWriter class which contains the method for writing output on NetCDF file
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

    //! Ids for x and t and NetCDF file
    int n_XId, n_TId, m_ncId;

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

    //! Create the file with coordinates, time and variable which is definitions
    void nDefineNetCDFFile(const CSimulation* m_pSimulation);

    //! Return the metadata of every variable
    string strGetVariableMetadata(const string& strVariable, const string& strField);

    //! Save data of every variable into NetCDF file
    void nSetOutputData(CSimulation* m_pSimulation) const;

    //! Close the NetCDF file
    void nCloseNetCDFFile(CSimulation* m_pSimulation) const;
};

#endif // DATA_WRITER_H