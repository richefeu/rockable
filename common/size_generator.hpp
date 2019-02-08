#ifndef SIZE_GENERATOR_HPP_F8C5824A
#define SIZE_GENERATOR_HPP_F8C5824A

/// @file
/// @date 2012 May
/// @brief classes to generate size distributions
/// @author Vincent Richefeu <Vincent.Richefeu@3sr-grenoble.fr>,
/// Lab 3SR, Grenoble University
/// @deprecated use generator from the std instead


class size_generator
{
public:
	virtual double get_size() = 0;
};
size_generator * SizeGenerator;


class constant_size : public size_generator
{
public:
	double size;
	static constant_size Instance;
	double get_size()
	{
		return size;
	}
private:
	constant_size(): size(1.0) { }
};
constant_size constant_size::Instance;


class uniform_size : public size_generator
{
public:
	double size_min, size_max;
	static uniform_size Instance;
	double get_size()
	{
		return (size_min + ran1(&seed1) * (size_max - size_min));
	}
private:
	uniform_size(): size_min(1.0), size_max(1.0) { }
};
uniform_size uniform_size::Instance;

#endif /* end of include guard: SIZE_GENERATOR_HPP_F8C5824A */
