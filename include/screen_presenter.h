/*!
*
 * \class CScreenPresenter
 * \brief
 * \details TODO 001 This is a more detailed description of the ScreenPresenter class
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file screen_presenter.h
 * \brief Contains CScreenPresenter definitions
 *
 */

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <sstream>
using std::ostream;
using std::ostringstream;

#include <fstream>
using std::ofstream;


class CScreenPresenter
    {
private:

    //! The SV folder
    string m_strSVDir;

    //! The main output file stream
    ofstream OutStream;

    //! System start-simulation time
    time_t m_tSysStartTime{};

    //! System finish-simulation time
    time_t m_tSysEndTime{};

    //! Last value returned by clock()
    double m_dClkLast{};
    //! Total elapsed CPU time
    double m_dCPUClock{};

    // Utility routines
    void StartClock();
    bool bFindExeDir(char const*);
    // bool bTimeToQuit(void);
    // static int nDoTimeUnits(string const*);
    // int nDoSimulationTimeMultiplier(string const*);
    // static double dGetTimeMultiplier(string const*);
    // static bool bParseDate(string const*, int&, int&, int&);
    // static bool bParseTime(string const*, int&, int&, int&);
    // void DoTimestepTotals(void);
    static string strGetBuild();
    static string strGetComputerName();


public:
    ofstream LogStream;

    // Intel x86, byte order is little-endian, 32-bit
    string const PLATFORM = "Intel x86/GNU C++";

    string const VERSION = "0.0.1";
    string const PROGRAM_NAME = "Saint Venant Equations Solver version 0.0.1 (15 Aug 2024)";
    string const PROGRAM_NAME_SHORT = "SV";


    string const NOTE = "      Note ";
    string const COPYRIGHT = "(C) 2024 Manuel Cobos";
    string const LINE = "-------------------------------------------------------------------------------";
    string const DISCLAIMER1 = "This program is distributed in the hope that it will be useful. but WITHOUT ANY";
    string const DISCLAIMER2 = "WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A";
    string const DISCLAIMER3 = "PARTICULAR PURPOSE. See the GNU General Public License for more details. You";
    string const DISCLAIMER4 = "should have received a copy of the GNU General Public License along with this";
    string const DISCLAIMER5 = "program; if not. contact the Free Software Foundation. Inc.. 675 Mass Ave.";
    string const DISCLAIMER6 = "Cambridge. MA 02139. USA.";

    string const START_NOTICE = "- Started on ";
    string const INITIALIZING_NOTICE = "- Initializing";
    string const RUN_END_NOTICE = "- Run ended at ";

    // Not likely that user will need to change these
    static constexpr int BUF_SIZE = 2048;                                     // Max length (inc. terminating NULL) of any C-type string

    // clock_t is a signed long: see <time.h>
    long const CLOCK_T_MIN = LONG_MIN;
    double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);

    CScreenPresenter();
    ~CScreenPresenter();

    //! Carries out init-of-simulation tidying (error messages etc.)
    void StartingRun(int, char const* []);
    void EndingRun();

    void AnnounceStart();
    void AnnounceLicence();
    // void AnnounceReadBasementDEM() const;
    // static void AnnounceAddLayers();
    // static void AnnounceReadRasterFiles();
    // static void AnnounceReadVectorFiles();
    // void AnnounceReadLGIS() const;
    // void AnnounceReadICGIS() const;
    // void AnnounceReadIHGIS() const;
    static void AnnounceInitializing();
    // void AnnounceReadInitialSuspSedGIS() const;
    // void AnnounceReadInitialFineUnconsSedGIS(int const) const;
    // void AnnounceReadInitialSandUnconsSedGIS(int const) const;
    // void AnnounceReadInitialCoarseUnconsSedGIS(int const) const;
    // void AnnounceReadInitialFineConsSedGIS(int const) const;
    // void AnnounceReadInitialSandConsSedGIS(int const) const;
    // void AnnounceReadInitialCoarseConsSedGIS(int const) const;
    // void AnnounceReadDeepWaterWaveValuesGIS() const;
    // void AnnounceReadSedimentEventInputValuesGIS() const;
    // void AnnounceReadFloodLocationGIS(void) const;
    // void AnnounceReadTideData(void) const;
    // static void AnnounceReadSCAPEShapeFunctionFile(void);
    // static void AnnounceAllocateMemory(void);
    // static void AnnounceIsRunning(void);
    // static void AnnounceSimEnd(void);
};
