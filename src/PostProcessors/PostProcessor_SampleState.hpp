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

#ifndef POSTPROCESSOR_SAMPLESTATE_HPP
#define POSTPROCESSOR_SAMPLESTATE_HPP

#include <fstream>
#include <vector>

#include "PostProcessor.hpp"

// The state of a sample in each conf-file, on one line of SampleState.txt:
// time, sizes and volume of the sample, solid fraction, coordination number,
// and the stress tensor computed from the contact forces.
class SampleState : public PostProcessor {
 public:
  SampleState();
  void init() override;
  void end() override;
  void read(std::istream& is) override;

  void exec() override;

 private:
  std::vector<size_t> wallBox;  // xmin, xmax, ymin, ymax, zmin, zmax (body numbers), if the sample is in a box
  double Volume{0.0};           // imposed volume of the sample (0 = from the box of walls or the periodic cell)
  std::ofstream file;
};

#endif /* end of include guard: POSTPROCESSOR_SAMPLESTATE_HPP */
