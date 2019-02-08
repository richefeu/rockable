#ifndef DATATABLE_HPP_60C51008
#define DATATABLE_HPP_60C51008

/// @file
/// @brief Classes to store and manage parameters of groups and between groups
/// @author Vincent Richefeu <Vincent.Richefeu@3sr-grenoble.fr>,
/// Lab 3SR, Grenoble University

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <set>

/**
Usage example:
@code
DataTable dt;

dt.set("mu", 0, 0, 0.7);
dt.set("kn", 0, 1, 1000);
dt.set("kn", 1, 2, 200.8);

dt.write(std::cout);
@endcode
*/


/// @brief A class for definition and storage of data for interaction between each group.
/// The class is designed for a rapid access to the data with the method get
class DataTable
{
public:
	size_t ngroup;
	std::vector <std::vector <std::vector <double> > > tables;
	std::map <std::string, size_t> data_id;
	std::set < std::tuple<size_t,size_t,size_t> > defined;

public:
	
	DataTable() { set_ngroup(1); }

	void clear() {
		for (size_t t = 0 ; t < tables.size() ; ++t) {
			for (size_t i = 0 ; i < tables[t].size() ; ++i) tables[t][i].clear();
			tables[t].clear();
		}
		tables.clear();
		data_id.clear();
		set_ngroup(1);
	}

	void set_ngroup(size_t n) {
		ngroup = n;
		for (size_t t = 0 ; t < tables.size() ; ++t) {
			tables[t].resize(ngroup);
			for (size_t i = 0 ; i < tables[t].size() ; ++i) tables[t][i].resize(ngroup, 0.0);
		}
	}

	size_t get_ngroup() const {
		return ngroup;
	}

	bool exists (const std::string & name) const {
		std::map<std::string, size_t >::const_iterator ip = data_id.find(name);
		if (ip == data_id.end()) {
			return false;
		}
		return true;
	}

	size_t add (const std::string & name) {
		std::map<std::string, size_t >::const_iterator ip = data_id.find(name);
		if (ip != data_id.end()) {
			return (size_t) ip->second;
		}

		data_id[name] = tables.size();
		std::vector<std::vector<double> > table;
		table.resize(ngroup);
		for (size_t i = 0; i < table.size(); i++) table[i].resize(ngroup, 0.0);
		tables.push_back(table);
		return (tables.size() - 1);
	}

	size_t get_id (const std::string & name) const {
		std::map<std::string, size_t >::const_iterator ip = data_id.find(name);
		if (ip != data_id.end()) {
			return (size_t) ip->second;
		}
		else {
			return 0;
		}
	}

	double get (size_t id, size_t g1, size_t g2) const {
		return tables[id][g1][g2];
	}
	
	bool isDefined(size_t id, size_t g1, size_t g2) const {
		return ( defined.find( std::tuple<size_t,size_t,size_t>(id, g1, g2) ) != defined.end() );
	}
	
	void set (size_t id, size_t g1, size_t g2, double val) {
		if (g1 >= ngroup) set_ngroup(g1 + 1);
		if (g2 >= ngroup) set_ngroup(g2 + 1);
		if (id < tables.size()) {
			tables[id][g1][g2] = val;
			tables[id][g2][g1] = val;
			defined.insert( std::tuple<size_t,size_t,size_t>(id, g1, g2) );
		}
	}
	
	// A self-add method to set a parameter
	size_t set (const std::string & name, size_t g1, size_t g2, double val) {
		size_t id = add(name);
		set (id, g1, g2, val);
		return id;
	}
};

#endif /* end of include guard: DATATABLE_HPP_60C51008 */
