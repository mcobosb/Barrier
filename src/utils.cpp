

#include <error_handling.h>
#include <string>
using std::string;

#include <vector>
using std::vector;

#include <sstream>
using std::ostream;
using std::stringstream;

#include <fstream>
using std::ofstream;

#include <algorithm>
using std::transform;

#include "utils.h"



/**
 * @brief Convert string to lowercase (creates a copy)
 * 
 * @param strIn Pointer to input string (left unchanged)
 * @return Lowercase version of input
 * 
 * @note Uses std::transform with tolower
 */
string strToLower(string const *strIn)
{
    string strOut = *strIn;
    transform(strIn->begin(), strIn->end(), strOut.begin(), tolower);
    return strOut;
}

/**
 * @brief Split string by delimiter (Python-like split() function)
 * 
 * Helper function that populates a pre-constructed vector.
 * Empty items are ignored (e.g., "a,,b" with ',' delimiter gives ["a", "b"]).
 * 
 * @param s Pointer to input string
 * @param delim Delimiter character (e.g., ',', ' ', ':')
 * @param elems Pointer to output vector (modified in-place)
 * @return Reference to elems vector
 * 
 * @note Source: stackoverflow.com/questions/236129/split-a-string-in-c
 * @see VstrSplit() for version that returns new vector
 */
vector<string> VstrSplit2(string const *s, char const delim, vector<string> *elems)
{
    stringstream ss(*s);
    string item;
    while (getline(ss, item, delim))
    {
        if (!item.empty())
            elems->push_back(item);
    }
    return *elems;
}

/**
 * @brief Split string by delimiter (convenience wrapper)
 * 
 * Returns a new vector of split strings. Calls VstrSplit2() internally.
 * 
 * Example: VstrSplit("a,b,c", ',') returns ["a", "b", "c"]
 * 
 * @param s Pointer to input string
 * @param delim Delimiter character
 * @return Vector of split strings (empty items ignored)
 */
vector<string> VstrSplit(string const *s, char const delim)
{
    vector<string> elems;
    VstrSplit2(s, delim, &elems);
    return elems;
}


/**
 * @brief Validate if string can be parsed as double
 * 
 * Uses std::istringstream to attempt parsing.
 * 
 * @param str String to validate (e.g., "3.14", "-1.5e-3")
 * @return true if valid double, false otherwise
 * 
 * @warning Does NOT detect trailing garbage (e.g., "3.14abc" returns true)
 * @note Source: stackoverflow.com/questions/392981
 */
bool bIsStringValidDouble(string& str)
{
    std::istringstream iStr(str);
    double dDummy;

    if (!(iStr >> dDummy))
        return false;

    return true;
}

/**
 * @brief Validate if string can be parsed as integer
 * 
 * Trims leading whitespace, removes optional sign (+/-), then checks for digits only.
 * 
 * @param str String to validate (modified: whitespace/sign removed)
 * @return true if valid integer, false otherwise
 * 
 * @note Accepts: "42", "-17", "+3", "  123"
 * @note Rejects: "3.14", "1e5", "12abc"
 * @note Source: stackoverflow.com/questions/2844817
 */
bool bIsStringValidInt(string& str)
{
    // Trim leading whitespace
    size_t nPos = str.find_first_not_of(" \t");
    if (nPos != string::npos)
        str = str.substr(nPos);

    // If the first character is the sign, remove it
    if ((str[0] == '-') || (str[0] == '+'))
        str.erase(0, 1);

    // Now check that the string contains only numbers
    return (str.find_first_not_of("0123456789") == string::npos);
}

/**
 * @brief Find minimum value in vector
 * 
 * Simple O(n) linear search.
 * 
 * @param vec Input vector (must not be empty)
 * @return Minimum value
 * 
 * @warning Undefined behavior if vector is empty (no check)
 * @note Used in TVD flux limiters (Superbee, MinMod)
 */
double dMinVectorValue(const vector<double>& vec) {
    // if (vec.empty()) {
    //     std::cerr << "Vector is empty." << std::endl;
    //     return std::numeric_limits<int>::max(); // Return maximum possible int value
    // }

    double minValue = vec[0]; // Starting value

    for (size_t i = 1; i < vec.size(); ++i) {
        if (vec[i] < minValue) {
            minValue = vec[i]; // Update the min value
        }
    }

    return minValue;
}

/**
 * @brief Find maximum value in vector
 * 
 * Simple O(n) linear search.
 * 
 * @param vec Input vector (must not be empty)
 * @return Maximum value
 * 
 * @warning Undefined behavior if vector is empty (no check)
 * @note Used in TVD flux limiters (Superbee)
 */
double dMaxVectorValue(const vector<double>& vec) {
    // if (vec.empty()) {
    //     std::cerr << "Vector is empty." << std::endl;
    //     return std::numeric_limits<int>::max(); // Return maximum possible int value
    // }

    double maxValue = vec[0]; // Starting value

    for (size_t i = 1; i < vec.size(); ++i) {
        if (vec[i] > maxValue) {
            maxValue = vec[i]; // Update the max value
        }
    }

    return maxValue;
}

