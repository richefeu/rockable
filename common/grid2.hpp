#ifndef GRID2_HPP_5BC47A5E
#define GRID2_HPP_5BC47A5E

/// @file
/// @brief A two-dimentional grid (Simplified 2D matrix)
/// @author Vincent Richefeu <Vincent.Richefeu@3sr-grenoble.fr>,
/// Lab 3SR, Grenoble University

/// Grid nrows x ncols
template <class T>
class grid2
{
public:
	vector<T> Data;
	uint nrows, ncols;

	grid2 () : nrows(0), ncols(0) { }
	grid2 (uint r, uint c = 1) : Data(r * c, 0), nrows(r), ncols(c) { }

	void resize(uint r, uint c = 1) {
		Data.resize(r * c);
		nrows = r;
		ncols = c;
	}
	
	void fill(const T& val = 0) {
		std::fill(Data.begin(), Data.end(), val);
	}

	T& operator()(uint r, uint c) {
		return Data[c * nrows + r];
	}
	
	const T& operator()(uint r, uint c) const {
		return Data[c * nrows + r];
	}

	// input/output
	friend ostream& operator << (ostream& pStr, const grid2 & Grd) {
		for (uint i = 0 ; i < Grd.Data.size() - 1 ; ++i) pStr << Grd.Data[i] << ' ';
		pStr << Grd.Data[Grd.Data.size() - 1];
		return pStr;
	}

	friend istream& operator >> (istream& pStr, grid2 & Grd) {
		T val;
		Grd.Data.clear();
		for (size_t i = 0 ; i < Grd.Data.size() - 1 ; ++i) {
			pStr >> val;
			Grd.Data.push_back(val);
		}
		return pStr;
	}
};

typedef grid2<double>       grid2r;
typedef grid2<int>          grid2i;
typedef grid2<unsigned int> grid2ui;
typedef grid2<bool>         grid2b;

#endif /* end of include guard: GRID2_HPP_5BC47A5E */
