#ifndef UTILS_H
#define UTILS_H

#include <string>
using std::string;

#include <vector>
using std::vector;


//===================================================== hard-wired constants ====================================================
constexpr char COLON = ':';
constexpr char COMMA = ',';
constexpr char DASH = '-';
constexpr char PATH_SEPARATOR = '/'; // Works for Windows too!
constexpr char QUOTE1 = ';';
constexpr char QUOTE2 = '#';
constexpr char SLASH = '/';
constexpr char SPACE = ' ';
constexpr char TILDE = '~';

// Time unit codes
constexpr int TIME_UNKNOWN = -1;
constexpr int TIME_SECONDS = 0;
constexpr int TIME_HOURS = 1;
constexpr int TIME_DAYS = 2;
constexpr int TIME_MONTHS = 3;
constexpr int TIME_YEARS = 4;


bool bIsStringValidDouble(string &);
bool bIsStringValidInt(string &);
double dMinVectorValue(const vector<double> &);
double dMaxVectorValue(const vector<double> &);

string strToLower(string const*);
vector<string> VstrSplit2(string const*, char, vector<string>*);
vector<string> VstrSplit(string const*, char);
#endif // UTILS_H
