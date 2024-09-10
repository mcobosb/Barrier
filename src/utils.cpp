

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



//===============================================================================================================================
//! Returns the lower case version of an string, leaving the original unchanged
//===============================================================================================================================
string strToLower(string const *strIn)
{
    string strOut = *strIn;
    transform(strIn->begin(), strIn->end(), strOut.begin(), tolower);
    return strOut;
}

//===============================================================================================================================
//! From http://stackoverflow.com/questions/236129/split-a-string-in-c They implement (approximately) Python's split() function. This first version puts the results into a pre-constructed string vector. It ignores empty items
//===============================================================================================================================
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

//===============================================================================================================================
//! From http://stackoverflow.com/questions/236129/split-a-string-in-c They implement (approximately) Python's split() function. This second version returns a new string vector (it calls the first version)
//===============================================================================================================================
vector<string> VstrSplit(string const *s, char const delim)
{
    vector<string> elems;
    VstrSplit2(s, delim, &elems);
    return elems;
}


//===============================================================================================================================
//! Checks to see if a string can be read as a valid double number. Does not find trailing (i.e.post-number) rubbish, but then neither does strtod(). From https://stackoverflow.com/questions/392981/how-can-i-convert-string-to-double-in-c
//===============================================================================================================================
bool bIsStringValidDouble(string& str)
{
    std::istringstream iStr(str);
    double dDummy;

    if (!(iStr >> dDummy))
        return false;

    return true;
}

//===============================================================================================================================
//! Checks to see if a string can be read as a valid integer, from https://stackoverflow.com/questions/2844817/how-do-i-check-if-a-c-string-is-an-int
//===============================================================================================================================
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

//===============================================================================================================================
//! Compute the min value of a vector
//===============================================================================================================================
double dMinVectorValue(const vector<double>& vec) {
    // if (vec.empty()) {
    //     std::cerr << "El vector está vacío." << std::endl;
    //     return std::numeric_limits<int>::max(); // Retorna el mayor valor posible de int
    // }

    double minValue = vec[0]; // Starting value

    for (int i = 1; i < vec.size(); ++i) {
        if (vec[i] < minValue) {
            minValue = vec[i]; // Update the min value
        }
    }

    return minValue;
}

//===============================================================================================================================
//! Compute the max value of a vector
//===============================================================================================================================
double dMaxVectorValue(const vector<double>& vec) {
    // if (vec.empty()) {
    //     std::cerr << "El vector está vacío." << std::endl;
    //     return std::numeric_limits<int>::max(); // Retorna el mayor valor posible de int
    // }

    double maxValue = vec[0]; // Starting value

    for (int i = 1; i < vec.size(); ++i) {
        if (vec[i] > maxValue) {
            maxValue = vec[i]; // Update the max value
        }
    }

    return maxValue;
}