#ifndef EXECCHRONO_HPP_22940E83
#define EXECCHRONO_HPP_22940E83

#include <ctime>
#include <cmath>
#include <string>
#include <iostream>
#include "message.hpp"

/**
@file ExecChrono.hpp
Usage:
@code{.cpp}
ExecChrono MM; // the chrono will start automatically (but we can use 'start' if we want)
// Doing something...
MM.stop();
@endcode
*/
class ExecChrono
{
private:
	std::string m_label;
	clock_t m_start, m_end;
	
public:
	ExecChrono():m_label("Measured time") { start(); }
	ExecChrono(const char * label):m_label(label) { start(); }
	void start() { m_start = clock(); }
	void stop() { 
		m_end = clock();
		
		// Now print it
		double ts = (double)(m_end - m_start) / CLOCKS_PER_SEC;
		//double ts_hour = floor(ts / 3600.0);
		//double ts_min  = floor((ts - ts_hour * 3600.0) / 60.0);
		//double ts_sec  = ts - 3600.0 * ts_hour - 60.0 * ts_min;
		
		std::cout << m_label << ": " << msg::HumanReadableSeconds(ts) << std::endl; 
		//<< ts_hour << "h " << ts_min << "m " << ts_sec << "s" << std::endl;
	}
};

#endif /* end of include guard: EXECCHRONO_HPP_22940E83 */
