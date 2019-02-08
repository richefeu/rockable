// Credits
// This class is issued from a tutorial course of Krishna Kumar 
// We acknowledge and are grateful to this developer for his contributions.

#ifndef FACTORY_HPP_8E5A4315
#define FACTORY_HPP_8E5A4315

#include <memory>
#include <string>
#include <map>
#include <functional>

/// The factory - implements singleton pattern!
template <
    class    BaseClass,
    typename Key = std::string
    >
class Factory
{
public:
	/// Get the single instance of the factory
	static Factory* Instance()
	{
		static Factory factory;
		return &factory;
	}

	/// register a factory function to create an instance of className
	void RegisterFactoryFunction(Key key, std::function<BaseClass*(void)> classFactoryFunction)
	{
		// register the class factory function
		factoryFunctionRegistry[key] = classFactoryFunction;
	}

	/// create an instance of a registered class
	BaseClass* Create(Key key)
	{
		// find name in the registry and call factory method.
		auto it = factoryFunctionRegistry.find(key);
		if (it != factoryFunctionRegistry.end()) return it->second();
		return nullptr;
	}

private:
	/// a private ctor
	Factory() { }

	/// the registry of factory functions
	std::map<Key, std::function<BaseClass*(void)> > factoryFunctionRegistry;
};

/// A helper class to register a factory function
template <
    class BaseClass,
    class DerivedClass,
    typename Key  = std::string
    >
class Registrar
{
public:
	explicit Registrar(Key key)
	{
		// register the class factory function
		Factory<BaseClass, Key>::Instance()->RegisterFactoryFunction(
		    key,
		    [](void) -> BaseClass* { return new DerivedClass(); }
		);
	}
};

#endif /* end of include guard: FACTORY_HPP_8E5A4315 */
