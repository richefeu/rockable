#ifndef COLORTABLE_HPP_B65A3038
#define COLORTABLE_HPP_B65A3038

#define ALL_BLACK 1
#define VIS5D 2
#define MATLAB_JET 3
#define SAMCET 4

#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib> // needed for rand()

using namespace std;

#define MIN(a,b) (((a)<(b))?(a):(b))

struct colorRGBA {
	int r, g, b, a;
	void set(int R, int G, int B, int A) {
		r = R; g = G; b = B; a = A;
	}
};

class ColorTable
{
private:

	bool big_endian;
	int size;
	float curvature, bias;
	int rotation;
	bool swap, invert;
	std::vector<unsigned int> table;
	float min, max;
	int tableID;
	float alpha, beta, alphapow;

	unsigned int PACK_COLOR(int R, int G, int B, int A) {
		if (big_endian) return ( (unsigned int)((R) << 24 | (G) << 16 | (B) << 8 | (A)) );
		else           return ( (unsigned int)((A) << 24 | (B) << 16 | (G) << 8 | (R)) );
	}
	
	int UNPACK_RED(unsigned int X) {
		if (big_endian) return ( ( (X) >> 24 ) & 0xff );
		else           return ( (X) & 0xff );
	}
	
	int UNPACK_GREEN(unsigned int X) {
		if (big_endian) return ( ( (X) >> 16 ) & 0xff );
		else           return ( ( (X) >> 8 ) & 0xff );
	}
	
	int UNPACK_BLUE(unsigned int X) {
		if (big_endian) return ( ( (X) >> 8 ) & 0xff );
		else           return ( ( (X) >> 16 ) & 0xff );
	}
	
	int UNPACK_ALPHA(unsigned int X) {
		if (big_endian) return ( (X) & 0xff );
		else           return ( ( (X) >> 24 ) & 0xff );
	}

	double gray(double s) {
		return s < 0. ? 0. : (s < 1. ? s : 1.);
	}
	
	double hot_r(double s) {
		return s < 0. ? 0. : (s < 3. / 8. ? 8. / 3. * s : 1.);
	}
	
	double hot_g(double s) {
		return s < 3. / 8. ? 0. : (s < 6. / 8. ? 8. / 3. * (s - 3. / 8.) : 1.);
	}
	
	double hot_b(double s) {
		return s < 6. / 8. ? 0. : (s < 1. ? 8. / 2. * (s - 6. / 8.) : 1.);
	}
	
	double cubic(double a, double b, double c, double d, double x) {
		return a + b * x + c * x * x + d * x * x * x;
	}

	void HSV_to_RGB(double H, double S, double V,
	                double *R, double *G, double *B) {
		if (S < 5.0e-6) {
			*R = *G = *B = V;
		}
		else {
			int i = (int)H;
			double f = H - (float)i;
			double p1 = V * (1.0 - S);
			double p2 = V * (1.0 - S * f);
			double p3 = V * (1.0 - S * (1.0 - f));
			switch (i) {
			case 0:
				*R = V;
				*G = p3;
				*B = p1;
				break;
			case 1:
				*R = p2;
				*G = V;
				*B = p1;
				break;
			case 2:
				*R = p1;
				*G = V;
				*B = p3;
				break;
			case 3:
				*R = p1;
				*G = p2;
				*B = V;
				break;
			case 4:
				*R = p3;
				*G = p1;
				*B = V;
				break;
			case 5:
				*R = V;
				*G = p1;
				*B = p2;
				break;
			}
		}
	}

	void RGB_to_HSV(double R, double G, double B,
	                double *H, double *S, double *V) {
		double maxv = R > G ? R : G;
		if (B > maxv) maxv = B;
		*V = maxv;
		if (maxv > 0) {
			double minv = R < G ? R : G;
			if (B < minv) minv = B;
			*S = 1.0 - double(minv) / maxv;
			if (maxv > minv) {
				if (maxv == R) {
					*H = (G - B) / double(maxv - minv);
					if (*H < 0) *H += 6.0;
				}
				else if (maxv == G) *H = 2.0 + (B - R) / double(maxv - minv);
				else *H = 4.0 + (R - G) / double(maxv - minv);
			}
		}
	}

public:

	ColorTable() {
		short int word = 0x0001;
		char *byte = (char*)&word;
		big_endian = (byte[0] ? false : true);

		rotation = 0.0;
		bias = 0.0;
		curvature = 0.0;
		min = 0.0;
		max = 1.0;
		size = 256;
		swap = false;
		tableID = 2;
		invert = false;
		alphapow = 0.0;
		alpha = 0.0;
		beta = 0.0;

		Rebuild();
	}

	void Rebuild() {
		double s, t, gamma;
		int r, g, b, a;

		if (!table.empty()) table.clear();
		table.reserve(size);

		for (int i = 0; i < size; i++) {

			if (size > 1) {
				if (i + rotation < 0)
					s = (double)(i + rotation + size) / (double)(size - 1);
				else if (i + rotation > size - 1)
					s = (double)(i + rotation - size) / (double)(size - 1);
				else
					s = (double)(i + rotation) / (double)(size - 1);
			}
			else
				s = 0.;

			if (swap) s = 1.0 - s;

			switch (tableID) {
			case 0:  // all black
				r = g = b = 0;
				break;
			case 1:  // vis5d
				t = (curvature + 1.4) * (s - (1. + bias) / 2.);
				r = (int)(128.0 + 127.0 * atan(7.0 * t) / 1.57);
				g = (int)(128.0 + 127.0 * (2 * exp(-7 * t * t) - 1));
				b = (int)(128.0 + 127.0 * atan(-7.0 * t) / 1.57);
				break;
			case 2: { // matlab "jet"
				double ii = (double)(s - bias) * 128.;
				if (ii < 0) ii = 0;
				if (ii > 128) ii = 128;
				double rr =
				    ii <= 46 ? 0. :
				    ii >= 111 ? -0.03125 * (ii - 111) + 1. :
				    ii >= 78 ? 1. :
				    0.03125 * (ii - 46);
				double gg =
				    ii <= 14 || ii >= 111 ? 0. :
				    ii >= 79 ? -0.03125 * (ii - 111) :
				    ii <= 46 ? 0.03125 * (ii - 14) :
				    1.;
				double bb =
				    ii >= 79 ? 0. :
				    ii >= 47 ? -0.03125 * (ii - 79) :
				    ii <= 14 ? 0.03125 * (ii - 14) + 1. :
				    1.;
				r = (int)(rr * 255.);
				g = (int)(gg * 255.);
				b = (int)(bb * 255.);
			}
			break;
			case 3: // lucie, samcef (?)
				if (s - bias <= 0.) {
					r = 0;
					g = 0;
					b = 255;
				}
				else if (s - bias <= 0.40) {
					r = 0;
					g = (int)((s - bias) * 637.5);
					b = (int)(255. - (s - bias) * 637.5);
				}
				else if (s - bias <= 0.60) {
					r = (int)(1275. * (s - bias - 0.4));
					g = 255;
					b = 0;
				}
				else if (s - bias <= 1.) {
					r = 255;
					g = (int)(255. - 637.5 * (s - bias - 0.6));
					b = 0;
				}
				else {
					r = 255;
					g = 0;
					b = 0;
				}
				break;
			case 4: // rainbow
				if (s - bias <= 0.) {
					r = 0;
					g = 0;
					b = 255;
				}
				else if (s - bias <= 0.25 + curvature) {
					curvature = (curvature == -0.25) ? -0.26 : curvature;
					r = 0;
					g = (int)((s - bias) * (255. / (0.25 + curvature)));
					b = 255;
				}
				else if (s - bias <= 0.50) {
					curvature = (curvature == 0.25) ? 0.26 : curvature;
					r = 0;
					g = 255;
					b = (int)(255. - (255. / (0.25 - curvature)) * (s - bias - 0.25 - curvature));
				}
				else if (s - bias <= 0.75 - curvature) {
					curvature = (curvature == 0.25) ? 0.26 : curvature;
					r = (int)((s - bias - 0.5) * (255. / (0.25 - curvature)));
					g = 255;
					b = 0;
				}
				else if (s - bias <= 1.) {
					curvature = (curvature == -0.25) ? -0.26 : curvature;
					r = 255;
					g = (int)(255. - (255. / (0.25 + curvature)) * (s - bias - 0.75 + curvature));
					b = 0;
				}
				else {
					r = 255;
					g = 0;
					b = 0;
				}
				break;
			case 5: // emc2000 (rainbow with black and white)
				if (s - bias <= 0.) {
					r = 0;
					g = 0;
					b = 0;
				}
				else if (s - bias <= 0.2) {
					r = (int)(57 * (1 - 100 * ((s - bias) - 0.1) * ((s - bias) - 0.1)));
					g = 0;
					b = (int)((s - bias) * (255. / 0.2));
				}
				else if (s - bias <= 0.3624) {
					r = 0;
					g = (int)((s - bias - 0.2) * (255. / 0.1624));
					b = 255;
				}
				else if (s - bias <= 0.50) {
					r = 0;
					g = 255;
					b = (int)(255. - (255. / 0.1376) * (s - bias - 0.3624));
				}
				else if (s - bias <= 0.6376) {
					r = (int)((s - bias - 0.5) * (255. / 0.1376));
					g = 255;
					b = 0;
				}
				else if (s - bias <= 0.8) {
					r = 255;
					g = (int)(255. - (255. / 0.1624) * (s - bias - 0.6376));
					b = 0;
				}
				else if (s - bias <= 1.0) {
					r = 255;
					g = (int)((255. / 0.2) * (s - bias - 0.8));
					b = (int)(-3187.66 * (s - bias) * (s - bias) + 7012.76 * (s - bias) - 3570.61);
				}
				else {
					r = 255;
					g = 255;
					b = 255;
				}
				break;
			case 6: // darkblue->red->yellow->white
				r = (int)(255. * cubic(-0.0506169, 2.81633, -1.87033, 0.0524573, s - bias));
				g = (int)(255. * cubic(0.0485868, -1.26109, 6.3074, -4.12498, s - bias));
				b = (int)(255. * cubic(0.364662, 1.50814, -7.36756, 6.51847, s - bias));
				break;
			case 7: // matlab "hot"
				r = (int)(255. * hot_r(s - bias));
				g = (int)(255. * hot_g(s - bias));
				b = (int)(255. * hot_b(s - bias));
				break;
			case 8: // matlab "pink"
				r = (int)(255. * sqrt((2.*gray(s - bias) + hot_r(s - bias)) / 3.));
				g = (int)(255. * sqrt((2.*gray(s - bias) + hot_g(s - bias)) / 3.));
				b = (int)(255. * sqrt((2.*gray(s - bias) + hot_b(s - bias)) / 3.));
				break;
			case 9: // grayscale
				if (s - bias <= 0.) {
					r = g = b = 0;
				}
				else if (s - bias <= 1.) {
					r = g = b = (int)(255 * (1. - curvature) * (s - bias));
				}
				else {
					r = g = b = (int)(255 * (1. - curvature));
				}
				break;
			case 10: // all white
				r = g = b = 255;
				break;
			case 11: { // matlab "hsv"
				double H = 6. * s + 1.e-10, R = 0.0, G = 0.0, B = 0.0;
				HSV_to_RGB(H, 1., 1., &R, &G, &B);
				r = (int)(255 * R);
				g = (int)(255 * G);
				b = (int)(255 * B);
			}
			break;
			case 12: { // spectrum (truncated hsv)
				double H = 5. * s + 1.e-10, R = 0.0, G = 0.0, B = 0.0;
				HSV_to_RGB(H, 1., 1., &R, &G, &B);
				r = (int)(255 * R);
				g = (int)(255 * G);
				b = (int)(255 * B);
			}
			break;
			case 13: // matlab "bone"
				r = (int)(255. * (7.0*gray(s - bias) + hot_b(s - bias)) / 8.0);
				g = (int)(255. * (7.0*gray(s - bias) + hot_g(s - bias)) / 8.0);
				b = (int)(255. * (7.0*gray(s - bias) + hot_r(s - bias)) / 8.0);
				break;
			case 14: // matlab "spring"
				r = (int)(255. * 1.0);
				g = (int)(255. * gray(s - bias));
				b = (int)(255. * (1.0 - gray(s - bias)));
				break;
			case 15: // matlab "summer"
				r = (int)(255. * gray(s - bias));
				g = (int)(255. * (0.5 + gray(s - bias) / 2.0));
				b = (int)(255. * 0.4);
				break;
			case 16: // matlab "autumn"
				r = (int)(255. * 1.0);
				g = (int)(255. * gray(s - bias));
				b = (int)(255. * 0.0);
				break;
			case 17: // matlab "winter"
				r = (int)(255. * 0.0);
				g = (int)(255. * gray(s - bias));
				b = (int)(255. * (0.5 + (1.0 - gray(s - bias)) / 2.0));
				break;
			case 18: // matlab "cool"
				r = (int)(255. * gray(s - bias));
				g = (int)(255. * (1.0 - gray(s - bias)));
				b = (int)(255. * 1.0);
				break;
			case 19: // matlab "copper"
				r = (int)(255. * MIN(1., gray(s - bias) * 1.25));
				g = (int)(255. * MIN(1., gray(s - bias) * 0.7812));
				b = (int)(255. * MIN(1., gray(s - bias) * 0.4975));
				break;
			case 20: // Random colors
				r = rand() / (double)RAND_MAX * 255.0;
				g = rand() / (double)RAND_MAX * 255.0;
				b = rand() / (double)RAND_MAX * 255.0;
				break;
			default:
				r = g = b = 0;
				break;
			}

			float aa = 1.0;
			if (alphapow) aa = pow(float(s ? s : 1.0e-10), float(alphapow));
			a = (int)(255.0 * aa * alpha);

			if (beta) {
				if (beta > 0.0) gamma = 1.0 - beta;
				else gamma = 1.0 / (1.001 + beta); // beta is thresholded to [-1,1]
				r = (int)(255.0 * pow((float)r / 255.0, gamma));
				g = (int)(255.0 * pow((float)g / 255.0, gamma));
				b = (int)(255.0 * pow((float)b / 255.0, gamma));
			}

			if (invert) {
				r = 255 - r;
				g = 255 - g;
				b = 255 - b;
			}

			// clamp to [0,255]
			r = r < 0 ? 0 : (r > 255 ? 255 : r);
			g = g < 0 ? 0 : (g > 255 ? 255 : g);
			b = b < 0 ? 0 : (b > 255 ? 255 : b);
			a = a < 0 ? 0 : (a > 255 ? 255 : a);

			table[i] = PACK_COLOR(r, g, b, a);
		}

	}
	
	void setTableID(int id) { tableID = id; }
	void setSwap(bool s) { swap = s; }
	void setInvert(bool i) { invert = i; }
	void setMinMax(float Min, float Max) {
		min = Min;
		max = Max;
	}
	void setSize(float Size) { size = Size; }
	void setBias(float Bias) { bias = Bias; }
	void setCurvature(float Curv) { curvature = Curv; }
	void setRotation(int Rot) { rotation = Rot; }

	int getSize() { return size; }
	float getMin() { return min; }
	float getMax() { return max; }
	void getRGB(float value, colorRGBA* col) {
		unsigned int i;

		float pos = (value - min) / (max - min); // in ]0.0 1.0[
		i = (unsigned int)(floor(pos * (float)(size - 2))) + 1; // in [1 size-2]
		//cout << "value = " << value << endl;
		//cout << "i = " << i << endl;

		if      (value <= min) i = 0;
		else if (value >= max) i = size - 1;

		col->r = UNPACK_RED   (table[i]);
		col->g = UNPACK_GREEN (table[i]);
		col->b = UNPACK_BLUE  (table[i]);
		col->a = UNPACK_ALPHA (table[i]);
	}

	void Print() {
		int i, r, g, b, a;

		for (i = 0; i < size; i++) {
			r = UNPACK_RED(table[i]);
			g = UNPACK_GREEN(table[i]);
			b = UNPACK_BLUE(table[i]);
			a = UNPACK_ALPHA(table[i]);
			cout << i << " [" << r << ", " << g << ", " << b << ", " << a << "]" << endl;
		}
	}
};


#endif /* end of include guard: COLORTABLE_HPP_B65A3038 */
