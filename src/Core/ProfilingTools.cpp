#ifdef ROCKABLE_ENABLE_PROFILING
#include <MATools.hxx>
#endif /* ROCKABLE_ENABLE_PROFILING */

namespace RockableProfiler
{
	/* brief This function initializes timers if the ROCKABLE_ENABLE_PROFILING is set to on otherwise do nothing */ 
	void Initialize()
	{
#ifdef ROCKABLE_ENABLE_PROFILING
		MATools::MATimer::initialize();
#endif /* ROCKABLE_ENABLE_PROFILING */
	}

	/* brief This function finalizes timers if the ROCKABLE_ENABLE_PROFILING is set to on otherwise do nothing */ 
	void Finalize()
	{
#ifdef ROCKABLE_ENABLE_PROFILING
		MATools::MATimer::finalize();
#endif /* ROCKABLE_ENABLE_PROFILING */
	}
};
