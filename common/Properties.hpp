#ifndef PROPERTIES_HPP_60C51008
#define PROPERTIES_HPP_60C51008

/// @file
/// @brief Store and manage parameters of groups and between groups
/// @author Vincent Richefeu <Vincent.Richefeu@3sr-grenoble.fr>,
/// Lab 3SR, Grenoble University

#include <string>
#include <iostream>
#include <vector>
#include <map>

/// @brief Definition and storage of global properties (e.g. density) for each group
class Properties
{
public:
	size_t ngroup;
	std::vector<std::vector<double> > prop;
	std::map <std::string, size_t> prop_id;

public:
	
	Properties() {
		set_ngroup(1);
	}

	/// Clear all properties
	void clear() {
		for (size_t t = 0 ; t < prop.size() ; ++t) {
			prop[t].clear();
		}
		prop.clear();
		prop_id.clear();
		set_ngroup(1);
	}

	/// set the number of groups
	void set_ngroup(size_t n) {
		ngroup = n;
		for (size_t t = 0 ; t < prop.size() ; ++t) {
			prop[t].resize(ngroup, 0.0);
		}
	}

	bool exists(const std::string & name) const {
		std::map<std::string, size_t >::const_iterator ip = prop_id.find(name);
		if (ip == prop_id.end()) {
			return false;
		}
		return true;
	}

	size_t add(const std::string & name) {		
		std::map<std::string, size_t >::const_iterator ip = prop_id.find(name);
		if (ip != prop_id.end()) return (size_t) ip->second;
		
		prop_id[name] = prop.size();
		std::vector<double> g;
		g.resize(ngroup, 0.0);
		prop.push_back(g);
		return (prop.size() - 1);
	}

	size_t get_id (const std::string & name) const {
		std::map<std::string, size_t >::const_iterator ip = prop_id.find(name);
		if (ip != prop_id.end()) {
			return (size_t) ip->second;
		}
		else {
			return 0;
		}
	}

	double get(size_t id, size_t g) const {
		return prop[id][g];
	}

	void set(size_t id, size_t g, double val) {
		if (g >= ngroup) set_ngroup(g + 1);
		if (id < prop.size()) prop[id][g] = val;
	}
	
	size_t set(const std::string & name, size_t g, double val) {
		size_t id = add(name);
		set (id, g, val);
		return id;
	}
};

#endif /* end of include guard: PROPERTIES_HPP_60C51008 */
