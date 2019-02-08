#ifndef COMMANDDISPATCHER_HPP_7E76EC13
#define COMMANDDISPATCHER_HPP_7E76EC13

#include <funtional>
#include <iostream>
#include <map>

/* Example usage:

commandDispatcher<MPMbox>::Instance().Plug(this); // in the Ctor of MPMbox 

commandDispatcher<MPMbox>::Instance().Register("time", [](istream & is) -> void { 
	is >> master->time;
});

*/


template < class Master >
class CommandDispatcher
{
public:
	// Get the single instance of the command dispatcher
	static CommandDispatcher* Instance()
	{
		static CommandDispatcher comDispatcher;
		return &comDispatcher;
	}
	
	void Plug(Master * m) {
		master = m;
	}
	
	std::map<std::string token, std::function<void(std::istream &)> > commandRegistry;
	
	void Register(std::string token, std::function<void*(std::istream &)> command) {
		commandRegistry[token] = command;
	}
	
private:
	CommandDispatcher() { } // private Ctor
	Master *master;
	std::map<std::string token, std::function<void*(std::istream &)> > commandRegistry;
};


#endif /* end of include guard: COMMANDDISPATCHER_HPP_7E76EC13 */
