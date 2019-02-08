#ifndef MAT4SYM_HPP_0DEED3DC
#define MAT4SYM_HPP_0DEED3DC

class mat4sym {
	
public:
	
	double xx, xy;
	double     yy;

	mat4sym(): xx(0), xy(0), yy(0) { }
	mat4sym(double XX, double XY, double YY): xx(XX), xy(XY), yy(YY) { }
	
	// Constants
	static mat4sym zero() { return mat4sym(); }
	static mat4sym unit() { return mat4sym(1,0,1); }
	static mat4sym one()  { return mat4sym(1,1,1); }
	
	void reset() { xx = xy = yy = 0; }

	void eigenvalues(double & v1, double & v2, bool & swapped) const {
		v1 = 0.5 * (xx + yy) + sqrt((0.5 * (xx - yy)) * (0.5 * (xx - yy)) + xy * xy);
		v2 = 0.5 * (xx + yy) - sqrt((0.5 * (xx - yy)) * (0.5 * (xx - yy)) + xy * xy);
		if (v2 > v1) {
			double swap = v1;
			v1 = v2;
			v2 = swap;
			swapped = true;
		}
	}

	mat4sym& operator += (const mat4sym &a) {
		xx += a.xx;
		xy += a.xy;
		yy += a.yy;
		return *this;
	}
	
	mat4sym& operator -= (const mat4sym &a) {
		xx -= a.xx;
		xy -= a.xy;
		yy -= a.yy;
		return *this;
	}
	
	mat4sym& operator *= (double k) {
		xx *= k;
		xy *= k;
		yy *= k;
		return *this;
	}
	
	mat4sym& operator /= (double k) {
		xx /= k;
		xy /= k;
		yy /= k;
		return *this;
	}

	friend mat4sym operator + (const mat4sym &a, const mat4sym &b) {
		return mat4sym (a.xx + b.xx, a.xy + b.xy, a.yy + b.yy);
	}
	
	friend mat4sym operator - (const mat4sym &a, const mat4sym &b) {
		return mat4sym (a.xx - b.xx, a.xy - b.xy, a.yy - b.yy);
	}
	
	friend mat4sym operator - (const mat4sym &a) {
		return mat4sym (-a.xx, -a.xy, -a.yy);
	}
	
	friend mat4sym operator * (const mat4sym &a, double k) {
		return mat4sym (k * a.xx, k * a.xy, k * a.yy);
	}

	friend mat4sym operator * (double k, const mat4sym &a) {
		return mat4sym (k * a.xx, k * a.xy, k * a.yy);
	}

	friend mat4sym operator / (const mat4sym &a, double k) {
		return mat4sym (a.xx / k, a.xy / k, a.yy / k);
	}

	friend vec2r operator * (const mat4sym &a, const vec2r &b) {
		return vec2r(a.xx * b.x + a.xy * b.y, a.xy * b.x + a.yy * b.y);
	}

	friend vec2r operator * (const vec2r &b, const mat4sym &a) {
		return vec2r(a.xx * b.x + a.xy * b.y, a.xy * b.x + a.yy * b.y);
	}
	
	friend std::ostream & operator << (std::ostream& pStr, const mat4sym& pV) {
		return (pStr <<  pV.xx << CommBox().sep << pV.xy << CommBox().sep << pV.yy << std::endl);
	}
};

#endif /* end of include guard: MATSYM4_HPP_0DEED3DC */
