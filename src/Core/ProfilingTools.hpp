#pragma once

#ifdef ROCKABLE_ENABLE_PROFILING
#include <MATools.hxx>
#endif /* ROCKABLE_ENABLE_PROFILING */

namespace RockableProfiler
{
	/* brief This function initializes timers if the ROCKABLE_ENABLE_PROFILING is set to on otherwise do nothing */ 
	void Initialize();

	/* brief This function finalizes timers if the ROCKABLE_ENABLE_PROFILING is set to on otherwise do nothing */ 
	void Finalize();

#ifndef ROCKABLE_ENABLE_PROFILING
	#ifdef START_TIMER /* remove these 3 following lines when profiler will be removed in toofus */
		#undef START_TIMER
	#endif
	#define START_TIMER(XNAME) 
	#define Catch_Time_Section(X) 
	#define Catch_Nested_Time_Section(X) 

	/**
	 * @brief This function captures the runtime of a given section.
	 * @param [in] lambda section that the user wants to measure.
	 * @return The runtime of lambda 
	 */
	template<typename Lambda>
		double chrono_section(Lambda&& lambda_function) {return double(-1);}

	/**
	 * @brief This function captures the runtime of a given section and add it in the current MATimerNode named a_name.
	 * @param [in] a_name of the chrono section measured.
	 * @param [in] a_lambda is the section captured.
	 */
	template<typename Lambda>
		void add_capture_chrono_section(std::string a_name, Lambda&& a_lambda_function) {return;}
#endif /* ndef ROCKABLE_ENABLE_PROFILING */

	/* Create an object ProfilerManager at the begining of Rockable */ 
	class ProfilerManager
	{
		public:
			ProfilerManager() {RockableProfiler::Initialize();}
			void DisableTimeTable()
			{ 
#ifdef ROCKABLE_ENABLE_PROFILING
				MATools::MATimer::Optional::disable_print_timetable();
#endif /* ROCKABLE_ENABLE_PROFILING */
			}
			void DisableWriteFile()
			{
#ifdef ROCKABLE_ENABLE_PROFILING
				MATools::MATimer::Optional::disable_write_file();
#endif /* ROCKABLE_ENABLE_PROFILING */
			}
			~ProfilerManager() {RockableProfiler::Finalize();}
	};
};
