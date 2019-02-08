#ifndef UTIL_HPP_3C4E662C
#define UTIL_HPP_3C4E662C

#include <vector>
#include <sstream>

/// Create a std::vector initialized with a c-string (const char*)
template <class T>
std::vector<T> make_vector (const char * STR)
{
	std::vector<T> V;
	std::stringstream ss(STR);
	T value;
	ss >> value;
	while (ss) {
		V.push_back(value);
		ss >> value;
	}
	return V;
}

/// Show the elapsed time ts in secondes in the hour-minutes-seconds format
void seconds_hms(std::ostream & os, int ts)
{
	int ts_hour = (int)floor(ts / 3600.0);
	int ts_min = (int)floor((ts - ts_hour * 3600.0) / 60.0);
	int ts_sec = ts - 3600 * ts_hour - 60 * ts_min;

	os << ts_hour << "h" << ts_min  << "m" << ts_sec  << "s";
}

#endif /* end of include guard: UTIL_HPP_3C4E662C */

