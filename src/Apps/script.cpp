//  Copyright or © or Copr. Rockable
//
//  vincent.richefeu@3sr-grenoble.fr
//
//  This software is a computer program whose purpose is
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces.
//  It is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include "Core/Rockable.hpp"
#include <chaiscript/chaiscript.hpp>

int main(int argc, char const* argv[]) {
  INIT_TIMERS();

  // clang-format off
  chaiscript::ModulePtr const vec3r_module = chaiscript::ModulePtr(new chaiscript::Module);
  chaiscript::utility::add_class<vec3r>(
      *vec3r_module, "vec3r", {chaiscript::constructor<vec3r()>()},
      {
				
				{chaiscript::fun(&vec3r::x), "x"}, 
				{chaiscript::fun(&vec3r::y), "y"}, 
				{chaiscript::fun(&vec3r::z), "z"}

      }
	);

  chaiscript::ModulePtr const Rockable_module = chaiscript::ModulePtr(new chaiscript::Module);
  chaiscript::utility::add_class<Rockable>(
      *Rockable_module, "Rockable", { chaiscript::constructor<Rockable()>() },
      {

          {chaiscript::fun(&Rockable::showBanner), "showBanner"},
          //{chaiscript::fun(&Rockable::initParser), "initParser"},
          {chaiscript::fun(&Rockable::setOpenMPThreads), "setOpenMPThreads"},
          {chaiscript::fun(static_cast<void (Rockable::*)(int)>(&Rockable::setVerboseLevel)), "setVerboseLevel"},
          {chaiscript::fun(static_cast<void (Rockable::*)(const std::string&)>(&Rockable::setVerboseLevel)),
           "setVerboseLevel"},
          {chaiscript::fun(&Rockable::initialChecks), "initialChecks"},
          {chaiscript::fun(&Rockable::integrate), "integrate"},
          {chaiscript::fun(static_cast<void (Rockable::*)(int)>(&Rockable::saveConf)), "saveConf"},
          {chaiscript::fun(static_cast<void (Rockable::*)(const char*)>(&Rockable::saveConf)), "saveConf"},
          {chaiscript::fun(static_cast<void (Rockable::*)(int)>(&Rockable::loadConf)), "loadConf"},
          {chaiscript::fun(static_cast<void (Rockable::*)(const char*)>(&Rockable::loadConf)), "loadConf"},
          {chaiscript::fun(&Rockable::console_run), "console_run"},
          {chaiscript::fun(&Rockable::dt), "dt"},
          {chaiscript::fun(&Rockable::t), "t"},
          {chaiscript::fun(&Rockable::tmax), "tmax"},

          {chaiscript::fun(&Rockable::interVerlet), "interVerlet"},
          {chaiscript::fun(&Rockable::interConf), "interConf"},
          {chaiscript::fun(&Rockable::DVerlet), "DVerlet"},
          {chaiscript::fun(&Rockable::dVerlet), "dVerlet"},
          {chaiscript::fun(&Rockable::numericalDampingCoeff), "numericalDampingCoeff"},
          {chaiscript::fun(&Rockable::velocityBarrier), "velocityBarrier"},
          {chaiscript::fun(&Rockable::angularVelocityBarrier), "angularVelocityBarrier"},
          {chaiscript::fun(&Rockable::velocityBarrierExponent), "velocityBarrierExponent"},
          {chaiscript::fun(&Rockable::angularVelocityBarrierExponent), "angularVelocityBarrierExponent"},

          {chaiscript::fun(&Rockable::gravity), "gravity"}

      }
	);

  // clang-format on

  chaiscript::ChaiScript chai;
  chai.add(vec3r_module);
  chai.add(Rockable_module);

#ifdef _OPENMP
  omp_set_num_threads(1);
  fmt::print("OpenMP number of threads enforced to 1 (if not chaiscript will be slow)\n");
#endif

  if (argc == 2) {
    chai.eval_file(argv[1]);
  } else {
    std::cerr << "usage = " << argv[0] << " <script_file>\n";
  }

  PRINT_TIMERS("rockable_chai");
  return 0;
}
