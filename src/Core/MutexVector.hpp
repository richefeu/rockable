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

#pragma once

#include <mutex>
#include <vector>
#include <sstream>
#include <stdexcept>

/// Vector of Mutex
/// This compensate std::vector<std::mutex> is not resizable
class MutexVector {
 public:
  using Mutex = std::mutex;

  MutexVector() : data_{nullptr}, size_{0} {}
  MutexVector(std::size_t size) : size_{size} { data_ = new Mutex[size_]; }
  ~MutexVector() noexcept {
    if (data_) delete[] data_;
  }
  const Mutex& operator[](std::size_t index) const {
#ifndef NDEBUG
    checkBounds(index);
#endif
    return data_[index];
  }
  Mutex& operator[](std::size_t index) {
#ifndef NDEBUG
    checkBounds(index);
#endif
    return data_[index];
  }
  std::size_t size() const { return size_; }
  /// \warning resizing re-initialize the mutex
  void resize(unsigned long newSize) {
    if(data_) delete[] data_;
    data_ = new Mutex[newSize];
    size_ = newSize;
  }

 private:
  Mutex* data_;
  std::size_t size_;

  void checkBounds(std::size_t index) const {
    std::stringstream message;
    message << "MutexVector::checkBounds: index (which is  " << index << " ) >= this->size() (which is " << size_ << ")";
    if (index >= size_) throw std::out_of_range{message.str()};
  }
};
