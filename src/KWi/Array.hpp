//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License
//  version 3 or later (See COPYING file for the text of the licence)
//
#ifndef Array_H
#define Array_H

#include <cstdlib>
#include <iostream>
#include <stdint.h>

template <typename T> class Array {
private:
  // Stored data
  T *data;
  // Number of dimensions
  int ndim;
  // Lengths of the corresponding array dimensions
  int *shape;
  // Total number of elements
  int64_t size;

  int64_t getIndex(const int *indices);
  int64_t getPreviousIndex(const int *indices);

public:
  Array(int ndim, int *shape);
  int getDim();
  int *getShape();
  int getShape(int i);
  T *getData();
  int64_t getSize();
  void initialize(const T value);
  void print() const;
  void accumulate(const int *indices, const T value);
  void increment(const int *indices);
  T get(const int *indices);
  void set(const int *indices, const T value);
  ~Array();
};

#endif
