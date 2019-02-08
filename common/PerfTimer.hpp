#ifndef PERFTIMER_HPP_88457644
#define PERFTIMER_HPP_88457644







#include <ctime>
#include <cmath>
#include <chrono>

/**
@file PerfTimer.hpp
Example of usage:
@code{.cpp}
PerfTimer ptimer; // the chrono will start automatically
                  // (but we can use 'start' or 'reset' if we want)
	
// pointless loop to cause a delay
for (int x = 1; x < 10000000; ++x) {
	x * x * x;
}

std::cout << msg::HumanReadableSeconds(ptimer.getElapsedTimeSeconds()) << std::endl;
@endcode
*/
class PerfTimer
{
private:
	std::chrono::high_resolution_clock::time_point epoch;
	
public:
	PerfTimer() { reset(); }
	void reset() { epoch = std::chrono::high_resolution_clock::now(); }
	void start() { reset(); }
	double getElapsedTimeSeconds() { 
		auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - epoch).count();
		return (duration_ms * 0.001);
	}
	double getIntermediateElapsedTimeSeconds() {
		double ret = getElapsedTimeSeconds();
		reset();
		return ret;
	}
};

#endif /* end of include guard: PERFTIMER_HPP_88457644 */
