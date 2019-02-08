#ifndef MAT9SYM_HPP_929FBB93
#define MAT9SYM_HPP_929FBB93

/// Symetric matrix 3x3
template <class T>
class mat9sym
{
public:
	T xx, xy, xz;
	T     yy, yz;
	T         zz;

	mat9sym () : xx(0), xy(0), xz(0), yy(0), yz(0), zz(0) { }
	mat9sym (const T XX, const T XY, const T XZ, const T YY, const T YZ, const T ZZ )
		: xx(XX), xy(XY), xz(XZ), yy(YY), yz(YZ), zz(ZZ) { }


	mat9sym (const mat9<T> &M) : xx(M.xx), xy(0.5 * (M.xy + M.yx)), xz(0.5 * (M.xz + M.zx)), yy(M.yy), yz(0.5 * (M.yz + M.zy)), zz(M.zz) { }
	mat9sym (const mat9sym &M) : xx(M.xx), xy(M.xy), xz(M.xz), yy(M.yy), yz(M.yz), zz(M.zz) { }

	static mat9sym zero() { return mat9sym(); }
	static mat9sym one()  { return mat9sym(1, 1, 1, 1, 1, 1); }
	static mat9sym unit() { return mat9sym(1, 0, 0, 1, 0, 1); }

	void set_diag (const T XX, const T YY, const T ZZ ) {
		xx = XX;
		yy = YY;
		zz = ZZ;
	}

	// Arithmetic operations
	friend mat9sym operator + (const mat9sym &a, const mat9sym &b) {
		return mat9sym(a.xx + b.xx, a.xy + b.xy, a.xz + b.xz,
		                            a.yy + b.yy, a.yz + b.yz,
		                                         a.zz + b.zz
		              );
	}

	friend mat9sym operator - (const mat9sym &a, const mat9sym &b) {
		return mat9sym(a.xx - b.xx, a.xy - b.xy, a.xz - b.xz,
		                            a.yy - b.yy, a.yz - b.yz,
		                                         a.zz - b.zz
		              );
	}

	friend mat9sym operator - (const mat9sym &a) {
		return mat9sym(-a.xx, -a.xy, -a.xz,
		                      -a.yy, -a.yz,
		                             -a.zz
		              );
	}

	friend mat9sym operator * (const mat9sym &a, T k) {
		return mat9sym(k * a.xx, k * a.xy, k * a.xz,
		                         k * a.yy, k * a.yz,
		                                   k * a.zz
		              );
	}

	friend mat9sym operator * (const mat9sym &a, const mat9sym &b) {
		return mat9sym(a.xx * b.xx + a.xy * b.xy + a.xz * b.xz, a.xy * b.xx + a.yy * b.xy + a.yz * b.xz, a.xz * b.xx + a.yz * b.xy + a.zz * b.xz,
		                                                        a.xy * b.xy + a.yy * b.yy + a.yz * b.yz, a.xz * b.xy + a.yz * b.yy + a.zz * b.yz,
		                                                                                                 a.xz * b.xz + a.yz * b.yz + a.zz * b.zz
		              );
	}

	friend mat9sym operator * (T k, const mat9sym &a) {
		return mat9sym(
			k * a.xx, k * a.xy, k * a.xz,
			          k * a.yy, k * a.yz,
			                    k * a.zz
		);
	}

	friend vec3<T> operator * (const mat9sym &a, const vec3<T> &v) {
		return vec3<T>(
			a.xx * v.x + a.xy * v.y + a.xz * v.z,
			a.xy * v.x + a.yy * v.y + a.yz * v.z,
			a.xz * v.x + a.yz * v.y + a.zz * v.z
		);
	}

	friend mat9sym operator / (const mat9sym &a, T K) {
		T k = 0.0;
		if (K != 0.0) k = 1.0 / K;
		return mat9sym(
			k * a.xx, k * a.xy, k * a.xz,
			          k * a.yy, k * a.yz,
			                    k * a.zz
		);
	}
	
	/// Dyadic product (tensorial product or otimes)
	friend mat9sym dyadic_product(const vec3<T> &a, const vec3<T> &b) {
		return mat9sym(
			a.x * b.x, a.x * b.y, a.x * b.z,
			           a.y * b.y, a.y * b.z,
			                      a.z * b.z
		);
	}

	void operator += (const mat9sym &a) {
		xx += a.xx;
		xy += a.xy;
		xz += a.xz;
		yy += a.yy;
		yz += a.yz;
		zz += a.zz;
	}

	void operator -= (const mat9sym &a) {
		xx -= a.xx;
		xy -= a.xy;
		xz -= a.xz;
		yy -= a.yy;
		yz -= a.yz;
		zz -= a.zz;
	}

	void operator *= (T k) {
		xx *= k;
		xy *= k;
		xz *= k;
		yy *= k;
		yz *= k;
		zz *= k;
	}

	void operator /= (T K) {
		T k = 0.0;
		if (K != 0.0) k = 1.0 / K;
		xx *= k;
		xy *= k;
		xz *= k;
		yy *= k;
		yz *= k;
		zz *= k;
	}

	void setZero() {
		xx = xy = xz = yy = yz = zz = 0.0;
	}

	void setIdentity() {
		xx = yy = zz = 1.0;
		xy = xz = yz = 0.0;
	}

	T det() const {
		return (xx * (yy * zz - yz * yz) - xy * (xy * zz - yz * xz) + xz * (xy * yz - yy * xz));
	}

	mat9sym<T> get_inverse() const
	{
		double det = det();
		double invdet;
		if (fabs(det) < 1e-20) invdet = 0.0; // this is a choice. Why not!
		else                   invdet = 1.0 / det;
		return mat9sym<T>(
			(yy * zz - yz * yz) * invdet,
			-(xy * zz - yz * xz) * invdet,
			(xy * yz - xz * yy) * invdet,

			-(xy * zz - xz * yz) * invdet,
			(xx * zz - xz * xz) * invdet,
			-(xx * yz - xz * xy) * invdet,

			(xy * yz - xz * yy) * invdet,
			-(xx * yz - xy * xz) * invdet,
			(xx * yy - xy * xy) * invdet
		);
	}
		
/*
		/// Compute eigenvectors (stored as columns in V) and corresponding eigenvalues (D)
		/// by assuming the matrix is double and symmetric
		/// See section 11.1 of Numerical Recipes in C for more information.
		int sym_eigen (mat9<T> & V, vec3<T> & D) const
		{
			int rot = 0;
			//double tresh;
			vec3<T> B;
			vec3<T> Z;

			// Save the input matrix in orig, use new matrix inp
			mat9<T> A = *this;
			// Set vectors to the identity matrix
			V.identity();
			// Set B and values to the diagonal of the input matrix
			for (uint i = 0; i < 3; i++) B[i] = D[i] = A[i*3+i];

			// Rotate until off diagonal elements of input matrix are zero
			for (int sweep = 0; sweep++ < 50;) {
				double sum = fabs(A[0*3+1]) + fabs(A[0*3+2]) + fabs(A[1*3+2]);
				double thresh;

				if (fabs(sum) < 1.0e-15) return rot;

				thresh = (sweep < 4) ? sum * 0.2 / 9.0 : 0.0;  // First three sweeps?

				for (int p =  0; p < 2; p++)
					for (int q = p+1; q < 3; q++) {
						double g = 100.0 * fabs(A[p*3+q]);

						// After 4 sweeps, skip the rotation if the
						// off-diagonal element is small.
						if ((sweep > 4) && (g < 1.0e-15)) A[p*3+q] = 0.0;
						else if (fabs(A[p*3+q]) > thresh) {
							double h = D[q] - D[p];
							double c, s, t;  // cosine, sine, tangent of rotation angle
							double tau;

							if (g < 1.0e-20) t = A[p*3+q] / h;
							else {
								double theta = 0.5 * h / A[p*3+q];
								t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
								if (theta < 0.0) t = -t;
							}

							c = 1.0 / sqrt(1.0 + t*t);   // cosine of rotation angle
							s = t*c;                     // sine of rotation angle
							tau = s / (1.0 + c);

							h = t * A[p*3+q];
							Z[p] -= h;
							Z[q] += h;
							D[p] -= h;
							D[q] += h;
							A[p*3+q] = 0.0;

							// case of rotations 0 <= j < p-1
							for (int j = 0; j <= p-1; j++) {
								g = A[j*3+p] ; h = A[j*3+q];
								A[j*3+p] = g - s*(h + g*tau);
								A[j*3+q] = h + s*(g - h*tau);
							}

							// case of rotations p < j < q
							for (int j = p + 1;  j < q; j++) {
								g = A[p*3+j] ; h = A[j*3+q];
								A[p*3+j] = g - s*(h - g*tau);
								A[j*3+q] = h + s*(g - h*tau);
							}

							// case of rotations q < j < 3
							for (int j = q + 1; j < 3; j++) {
								g = A[p*3+j] ; h = A[q*3+j];
								A[p*3+j] = g - s*(h + g*tau);
								A[q*3+j] = h + s*(g - h*tau);
							}

							// Set the eigen vectors
							for (int j = 0; j < 3; j++) {
								g = V[j*3+p];
								h = V[j*3+q];
								V[j*3+p] = g - s*(h + g*tau);
								V[j*3+q] = h + s*(g - h*tau);
							}
							rot++;
						}
					}

				// Set the eigen values
				B += Z;
				D = B;
				Z.set(0.0, 0.0, 0.0);
			}
			return -1;  // Non-normal return - too many rotations
		}
	*/

	// input/output
	friend std::ostream& operator<< (std::ostream& pStr, const mat9sym& M) {
		return (pStr <<  M.xx << CommBox().sep << M.xy << CommBox().sep << M.xz << CommBox().sep
		                                       << M.yy << CommBox().sep << M.yz << CommBox().sep
		                                                                << M.zz );
	}

	friend std::istream& operator>> (std::istream& pStr, mat9sym& M) {
		return (pStr >> M.xx >> M.xy >> M.xz >> M.yy >> M.yz >> M.zz);
	}

};

typedef mat9sym<double>       mat9symr;
typedef mat9sym<float>        mat9symf;
typedef mat9sym<int>          mat9symi;
typedef mat9sym<unsigned int> mat9symui;
typedef mat9sym<bool>         mat9symb;

#endif /* end of include guard: MATSYM9_HPP_929FBB93 */
