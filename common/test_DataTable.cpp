#include "DataTable.hpp"
#include <iostream>

using namespace std;

int main ()
{
	DataTable dt;
	
	dt.set("mu", 0, 0, 0.7);
	dt.set("kn", 0, 1, 1000);
	dt.set("kn", 1, 2, 200.8);
	
	dt.write(std::cout);
	return 0;
}