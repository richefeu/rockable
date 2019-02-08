// Copyright (C) vec2 <vincent.richefeu@3sr-grenoble.fr>
// 
// This file is part of mbox.
// 
// vec2 can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note 
// Without a license, the code is copyrighted by default. 
// People can read the code, but they have no legal right to use it. 
// To use the code, you must contact the author directly and ask permission.

#ifndef VEC2_HPP_7A4024FD
#define VEC2_HPP_7A4024FD

/// @file   vec2.hpp
/// @brief  Template class for vectors with 2 components
/// @author Vincent Richefeu <Vincent.Richefeu@ujf-grenoble.fr>,
///         Lab 3SR, Grenoble University

#include <iostream>
#include <cmath>
#include "common.hpp"

/// @brief Vector with 2 components
template <typename T>
class vec2
{
public:
	T x, y;

	vec2(): x(0), y(0) { }
	vec2(T X, T Y) : x(X), y(Y) { }
	vec2(const vec2 &v) : x(v.x), y(v.y) { }
  
	vec2 & operator = (const vec2 &V) {
		x = V.x;
		y = V.y;
		return (*this);
	}

	static vec2 unit_x() { return vec2(1, 0); }
	static vec2 unit_y() { return vec2(0, 1); }
	static vec2 one() { return vec2(1, 1); }
	static vec2 zero() { return vec2(0, 0); }

	void reset() {
		x = y = 0;
	}

	void set(T X, T Y) {
		x = X;
		y = Y;
	}
	void set(T val) {
		x = y = val;
	}

	bool isnull(const T tol = 1e-20) const {
		return (fabs(x) < tol && fabs(y) < tol);
	}

	T* c_vec() { return &x; }

	T &operator [](int i) { return (&x)[i]; }

	const T &operator [](int i) const { return (&x)[i]; }

	// For local frames, the notation n,t and s is more appropriate than x and y
	const T n() const { return x; }
	const T t() const { return y; }

	// Arithmetic operations
	vec2& operator += (const vec2 &a) {
		x += a.x;
		y += a.y;
		return *this;
	}

	vec2& operator -= (const vec2 &a) {
		x -= a.x;
		y -= a.y;
		return *this;
	}

	vec2& operator *= (T k) {
		x *= k;
		y *= k;
		return *this;
	}

	vec2& operator /= (T k) {
		x /= k;
		y /= k;
		return *this;
	}

	friend vec2 operator + (const vec2 &a, const vec2 &b) {
		return vec2 (a.x + b.x, a.y + b.y);
	}

	friend vec2 operator - (const vec2 &a, const vec2 &b) {
		return vec2 (a.x - b.x, a.y - b.y);
	}

	friend vec2 operator - (const vec2 &a) {
		return vec2 (-a.x, -a.y);
	}

	friend vec2 operator * (const vec2 &a, T k) {
		return vec2 (a.x * k, a.y * k);
	}

	friend vec2 operator * (T k, const vec2 &a) {
		return vec2 (a.x * k, a.y * k);
	}

	friend vec2 operator / (const vec2 &a, T k) {
		return vec2 (a.x / k, a.y / k);
	}

	// --- Specific external operations ---

	/// Dot product
	friend T operator * (const vec2 &a, const vec2 &b) {
		return (a.x * b.x + a.y * b.y);
	}

	/// Multiply each component one another
	friend vec2<T> component_product(const vec2<T> &a, const vec2<T> &b) {
		return vec2<T>(a.x * b.x, a.y * b.y);
	}

	/// Find the smallest components
	friend vec2<T> component_min(const vec2<T> &a, const vec2<T> &b) {
		return vec2<T>(
			(a.x < b.x) ? a.x : b.x,
			(a.y < b.y) ? a.y : b.y
		);
	}

	/// Find the biggest components
	friend vec2<T> component_max(const vec2<T> &a, const vec2<T> &b) {
		return vec2<T>(
			(a.x > b.x) ? a.x : b.x,
			(a.y > b.y) ? a.y : b.y
		);
	}

	/// Absolut value of the components
	friend vec2<T> component_abs(const vec2<T> &a) {
		return vec2<T>(fabs(a.x), fabs(a.y));
	}

	/// Cross product
	friend T cross(const vec2<T> &a, const vec2<T> &b) {
		return (a.x * b.y - a.y * b.x);
	}

	/// Squared length of the vector
	friend T norm2(const vec2 &a) {
		return a * a;
	}

	/// Length of the vector
	friend T norm(const vec2 &a) { return sqrt(a * a); }

	T length() const { return norm(*this); }

	/// Normalize and return length (before being normalized)
	T normalize() {
		T n = norm(*this);
		if (n > 0.0) *this *= (1.0f / n);
		return n;
	}

	/// Normalize and return the normalized vector
	vec2 normalized() const {
		vec2 V = *this;
		V.normalize();
		return V;
	}

	/// Determinant
	friend T determinant(const vec2<T> a, const vec2<T> b) {
		return (a.x * b.y - b.x * a.y);
	}

	// Comparisons
    bool operator == (const vec2<T> & other) const {
		return (this->x == other.x && this->y == other.y);
	}

	bool operator != (const vec2<T> & other) const {
		return !(*this == other);
	}

	// input/output
	friend std::ostream & operator << (std::ostream& pStr, const vec2& pV) {
		return (pStr <<  pV.x << CommBox().sep << pV.y );
	}

	friend std::istream & operator >> (std::istream& pStr, vec2& pV) {
		return (pStr >> pV.x >> pV.y);
	}
};

typedef vec2<double>       vec2r;
typedef vec2<int>          vec2i;
typedef vec2<unsigned int> vec2ui;
typedef vec2<bool>         vec2b;

namespace std
{
template <class T>
struct less<vec2<T> > {
	bool operator() (const vec2<T> & lhs, const vec2<T> & rhs) const
	{
		if (rhs.x < lhs.x) return true;
		else if (lhs.x == rhs.x &&  rhs.y < lhs.y) return true;
		return false;
	}
};
}

#endif /* end of include guard: VEC2_HPP_7A4024FD */
