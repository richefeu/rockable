#ifndef FASTSORT3_HPP_29480B9B
#define FASTSORT3_HPP_29480B9B

/**
 @file
 Usage example:
 @code{.cpp}
 int a = 3, b = 1, c = 2;
 fastSort3<int>(a, b, c);
 @endcode
 The result is: a = 1, b = 2, c = 3
*/

template < typename T >
void fastSort3(T & a0, T & a1, T & a2)
{
	if (a0 < a1) {
		if (a1 < a2) { return; }
		else {
			if (a0 < a2) {
				std::swap(a1, a2);
			}
			else {
				std::swap(a0, a2);
				std::swap(a1, a2);
			}
		}
	}
	else {
		if (a0 < a2) {
			std::swap(a0, a1);
		}
		else {
			if (a1 < a2) {
				T temp = a0;
				a0 = a1;
				a1 = a2;
				a2 = temp;
			}
			else {
				std::swap(a0, a2);
			}
		}
	}
}

#endif /* end of include guard: FASTSORT3_HPP_29480B9B */
