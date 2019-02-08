#ifndef MAT4_HPP_A9AFC5DF
#define MAT4_HPP_A9AFC5DF

#include "vec2.hpp"

class mat4 {

public:

	double xx, xy;
	double yx, yy;

	mat4() : xx(0), xy(0), yx(0), yy(0) { }
	//mat4(mat4sym & m): xx(m.xx), xy(m.xy), yx(m.xy), yy(m.yy) { }
	mat4(double XX, double XY, double YX, double YY) : xx(XX), xy(XY), yx(YX), yy(YY) { }
	mat4 (const double M[]): xx(M[0]), xy(M[1]), yx(M[2]), yy(M[3]) { }

	// Constants
	static mat4 unit() { return mat4(1, 0, 0, 1); }
	static mat4 zero() { return mat4(1, 1, 1, 1); }
	static mat4 one() { return mat4(1, 1, 1, 1); }
	
	void reset() { xx = xy = yx = yy = 0; }
	
	void reset (const double val) {
		xx = xy = val;
		yx = yy = val;
	}

	void set_diag (const double XX, const double YY) {
		xx = XX;
		yy = YY;
	}
	
	mat4 transpose() { return mat4(xx, yx, xy, yy); }

	void eigenvalues(double & v1, double & v2, bool & swapped) const {
		/// @fixme seems to ok only for symmetric matrix
		v1 = 0.5 * (xx + yy) + sqrt((0.5 * (xx - yy)) * (0.5 * (xx - yy)) + xy * xy);
		v2 = 0.5 * (xx + yy) - sqrt((0.5 * (xx - yy)) * (0.5 * (xx - yy)) + xy * xy);
		if (v2 > v1) {
			double swap = v1;
			v1 = v2;
			v2 = swap;
			swapped = true;
		}
	}

	// works for non-symmetric matrices as well. Has problems. Maybe a tolerance should be included
	void eigen(mat4 & V, mat4 & D) {
	    double T = xx + yy;
	    double det = xx * yy - xy * yx;
	    double L1 = T / 2 + sqrt(T * T / 4 - det); // eigenval
	    double L2 = T / 2 - sqrt(T * T / 4 - det); // eigenval
	    D.xx = L1;
	    D.xy = D.yx = 0.0;
	    D.yy = L2;
	    // Eigenvectors organized vertically
	    if (yx != 0) {
	        V.xx = L1 - yy;
	        V.yx = yx;
	        V.xy = L2 - yy;
	        V.yy = yx;
	    }
	    else if (xy != 0) {
	        V.xx = xy;
	        V.yx = L1 - xx;
	        V.xy = xy;
	        V.yy = L2 - xx;
	    }
	    else if (xy == 0 and yx == 0) {
	        V.xx = 1;
	        V.yx = 0;
	        V.xy = 0;
	        V.yy = 1;
	    }
    
	    // Normalizing vectors
	    vec2r v1(V.xx, V.yx);
	    vec2r v2(V.xy, V.yy);
	    v1 = v1 / sqrt(v1.x * v1.x + v1.y * v1.y);  // use normalized()...
	    v2 = v2 / sqrt(v2.x * v2.x + v2.y * v2.y);

	    // Putting them back in the V matrix
	    V.reset();
	    V.xx = v1.x;
	    V.yx = v1.y;
	    V.xy = v2.x;
	    V.yy = v2.y;
	}

	int sym_eigen (mat4 & V, mat4 & D) const {
		int rot = 0;
		vec2r B;
		vec2r Z;

		// Save the input matrix in orig, use new matrix inp
		mat4 A = *this;
		// Set vectors to the identity matrix
		V.xx = 1;
		V.xy = 0;
		V.yx = 0;
		V.yy = 1;
		// Set B and D values to the diagonal of the input matrix
		B.x = D.xx = A.xx;
		B.y = D.yy = A.yy;

		// Rotate until off-diagonal elements of input matrix are zero
		for (int sweep = 0; sweep++ < 50;) {
			double sum = fabs(A.xy);
			double thresh;

			if (fabs(sum) < 1.0e-15) return rot;

			thresh = (sweep < 4) ? sum * 0.2 / 4.0 : 0.0;  // First three sweeps?
			double g = 100.0 * fabs(A.xy);  //TBC!!

			// After 4 sweeps, skip the rotation if the
			// off-diagonal element is small.
			if ((sweep > 4) && (g < 1.0e-15)) A.xy = 0.0;
			else if (fabs(A.xy) > thresh) {
					
				double h = D.yy - D.xx;
				double c, s, t;  // cosine, sine, tangent of rotation angle
				double tau;

				if (g < 1.0e-20) t = A.xy / h;
				else {
					double theta = 0.5 * h / A.xy;
					t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
					if (theta < 0.0) t = -t;
				}

				c = 1.0 / sqrt(1.0 + t * t); // cosine of rotation angle
				s = t * c;                   // sine of rotation angle
				tau = s / (1.0 + c);

				h = t * A.xy;
				Z.x -= h;
				Z.y += h;
				D.xx -= h;
				D.yy += h;
				A.xy = 0.0;
					
				g = V.xx;
				h = V.xy;
				V.xx = g - s * (h + g * tau);
				V.xy = h + s * (g - h * tau);
					
				g = V.yx;
				h = V.yy;
				V.yx = g - s * (h + g * tau);
				V.yy = h + s * (g - h * tau);
								
				rot++;
			}
		

			// Set the eigen values
			B += Z;
			D.xx = B.x;
			D.yy = B.y;
			//D = B;
			Z.x = 0.0;
			Z.y = 0.0;
		}
		return -1;  // Non-normal return - too many rotations
	}
	
	bool inverse() {
		double det = xx * yy - xy * yx;
		if (fabs(det) < 1.0e-20) return false; // inverse cannot be calculated

		double swap = xx;
		xx = yy;
		yy = swap;

		double inv_det = 1.0 / det;
		xx *=  inv_det;
		xy *= -inv_det;
		yy *=  inv_det;
		yx *= -inv_det;

		return true;
	}

	mat4 inverse2() {
    
	    double det = xx * yy - xy * yx;
		//if (fabs(det) < 1.0e-20) return false; // inverse cannot be calculated
	    double xx1(xx), xy1(xy), yx1(yx), yy1(yy);
		double swap = xx1;
		xx1 = yy1;
		yy1 = swap;

		double inv_det = 1.0 / det;
		xx1 *=  inv_det;
		xy1 *= -inv_det;
		yy1 *=  inv_det;
		yx1 *= -inv_det;

		return mat4(xx1, xy1, yx1, yy1);   
	}

	double det() const { return (xx * yy - xy * yx); }
	
	void svd ( mat4 & U, mat4 & S, mat4 & V) const {
	    // taken from http://www.lucidarme.me/?p=4802
	    // U matrix
	    double val1 = xx * yx + xy * yy;
	    double val2 = xx * xx + xy * xy - yx * yx - yy * yy;
	    double val3 = xx * xy + yx * yy;
	    double val4 = xx * xx - xy * xy + yx * yx - yy * yy;

	    double theta = 0.5 * atan2(2*val1, val2);
	    U.xx = cos(theta);
	    U.xy = -sin(theta);
	    U.yx = sin(theta);
	    U.yy = cos(theta);

	    //Singular value matrix (S)
	    double S1 = xx * xx + xy * xy + yx * yx + yy * yy;
	    double S2 = sqrt(val2 * val2 + 4.0 * val1 * val1);
	    //singular values
	    double sv1 = sqrt((S1 + S2) / 2.0);
	    double sv2 = sqrt((S1 - S2) / 2.0);

	    S.xx = sv1;
	    S.yy = sv2;
	    S.xy = 0;
	    S.yx = 0;

	    // V matrix
	    double phi = 0.5 * atan2(2.0 * val3, val4);
	    double s11 = (xx * cos(theta) + yx * sin(theta)) * cos(phi) + (xy * cos(theta) + yy * sin(theta)) * sin(phi);
	    double s22 = (xx * sin(theta) - yx * cos(theta)) * sin(phi) + (-xy * sin(theta) + yy * cos(theta)) * cos(phi);
	    V.xx = s11 / fabs(s11) * cos(phi);
	    V.xy = -s22 / fabs(s22) * sin(phi);
	    V.yx = s11 / fabs(s11) * sin(phi);
	    V.yy = s22 / fabs(s22) * cos(phi);
	}

	bool square_root(mat4 & SqR) const {
		double tau = xx + yy;
		double delta = xx * yy - xy * yx;
	
		if (delta == 0.0) return false;
	
		double s = sqrt(delta);
		double t = sqrt(tau + 2 * s);
	
		SqR.xx = (xx + s) / t;
		SqR.xy = (xy) / t;
		SqR.yx = (yx) / t;
		SqR.yy = (yy + s) / t;
	
		return true;
	}
    
	// =======================
	//  Arithmetic operations
	// =======================
	
	mat4& operator += (const mat4 &a) {
		xx += a.xx;
		xy += a.xy;
		yx += a.yx;
		yy += a.yy;
		return *this;
	}
	
	mat4& operator -= (const mat4 &a) {
		xx -= a.xx;
		xy -= a.xy;
		yx -= a.yx;
		yy -= a.yy;
		return *this;
	}
	
	mat4& operator *= (double k) {
		xx *= k;
		xy *= k;
		yx *= k;
		yy *= k;
		return *this;
	}

	mat4& operator /= (double k) {
		xx /= k;
		xy /= k;
		yx /= k;
		yy /= k;
		return *this;
	}

	// Comparisons
    bool operator == (const mat4 & other) const {
		return (
			   this->xx == other.xx && this->xy == other.xy 
			&& this->yx == other.yx && this->yy == other.yy 
		);
	}
	
	bool operator != (const mat4 & other) const {
		return !(*this == other);
	}

	// =========
	//  FRIENDS
	// =========
	
	friend mat4 operator + (const mat4 &a, const mat4 &b) {
		return mat4(
			a.xx + b.xx, a.xy + b.xy,
			a.yx + b.yx, a.yy + b.yy
		);
	}

	friend mat4 operator - (const mat4 &a, const mat4 &b) {
		return mat4(
			a.xx - b.xx, a.xy - b.xy,
			a.yx - b.yx, a.yy - b.yy
		);
	}

	friend mat4 operator - (const mat4 &a) {
		return mat4(
			-a.xx, -a.xy,
			-a.yx, -a.yy
		);
	}

	friend mat4 operator * (const mat4 &a, double k) {
		return mat4(
			k * a.xx, k * a.xy,
			k * a.yx, k * a.yy
		);
	}

	friend mat4 operator * (double k, const mat4 &a) {
		return mat4(
			k * a.xx, k * a.xy,
			k * a.yx, k * a.yy
		);
	}

	friend mat4 operator / (const mat4 &a, double k) {
		return mat4(
			a.xx / k, a.xy / k,
			a.yx / k, a.yy / k
		);
	}

	friend vec2r operator * (const mat4 &a, const vec2r &b) {
		return vec2r(a.xx * b.x + a.xy * b.y, a.yx * b.x + a.yy * b.y);
	}

	friend vec2r operator * (const vec2r &b, const mat4 &a) {
		return vec2r(a.xx * b.x + a.xy * b.y, a.yx * b.x + a.yy * b.y);
	}

	friend mat4 operator * (const mat4 &a, const mat4 &b) {
		return mat4(
			a.xx * b.xx + a.xy * b.yx, a.xx * b.xy + a.xy * b.yy,
			a.yx * b.xx + a.yy * b.yx, a.yx * b.xy + a.yy * b.yy
		);
	}

	friend std::ostream & operator << (std::ostream& pStr, const mat4& pV) {
		return (pStr <<  pV.xx << CommBox().sep << pV.xy << CommBox().sep << pV.yx << CommBox().sep << pV.yy << '\n');
	}
};





#endif /* end of include guard: MAT4_HPP_A9AFC5DF */
