//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License
//  version 3 or later (See COPYING file for the text of the licence)
//

#include "Array.hpp"
#include <cstdlib>
#include <iostream>
#include <stdint.h>

template <typename T> int64_t Array<T>::getIndex(const int *indices) {
  int64_t index = 0;
  for (int i = 0; i < this->ndim; i++) {
    // std::cout << "indices[ " << i << "] = " << indices[i]
    //  	      << " shape[ " << i << "] = " << this->shape[i] <<
    //  std::endl;
    index = index * this->shape[i] + (indices[i] - 1);
  }
  // std::cout << "index : " << index << std::endl;
  return index;
}

template <typename T> int64_t Array<T>::getPreviousIndex(const int *indices) {
  return this->getIndex(indices) - 1;
}

template <typename T>
Array<T>::Array(int ndim, int *shape) : ndim(ndim), size(1) {
  this->shape = (int *)(calloc(ndim, sizeof(int)));
  for (int i = 0; i < ndim; i++) {
    this->shape[i] = shape[i];
    size *= this->shape[i];
  }
  data = static_cast<T *>(calloc(size, sizeof(T)));
}

template <typename T> int Array<T>::getDim() { return this->ndim; }

template <typename T> int *Array<T>::getShape() { return this->shape; }

template <typename T> int Array<T>::getShape(int i) { return this->shape[i]; }

template <typename T> T *Array<T>::getData() { return this->data; }

template <typename T> int64_t Array<T>::getSize() { return this->size; }

template <typename T> void Array<T>::initialize(const T value) {
  for (int64_t i = 0; i < this->size; i++) {
    this->data[i] = value;
  }
}

template <typename T> void Array<T>::print() const {
  for (int64_t i = 0; i < this->size; i++) {
    std::cout << "data[" << i << "] = " << this->data[i] << '\t';
    std::cout << "*(data+" << i << ") = " << *(this->data + i) << '\n';
  }
  std::cout << std::endl;
}

template <typename T>
void Array<T>::accumulate(const int *indices, const T value) {
  this->data[this->getIndex(indices)] =
      this->data[this->getPreviousIndex(indices)] + value;
}

template <typename T> void Array<T>::increment(const int *indices) {
  this->data[this->getIndex(indices)]++;
}

template <typename T> T Array<T>::get(const int *indices) {
  return (this->data)[this->getIndex(indices)];
}

template <typename T> void Array<T>::set(const int *indices, const T value) {
  this->data[this->getIndex(indices)] = value;
}

template <typename T> Array<T>::~Array() {
  free(this->shape);
  free(this->data);
}

template class Array<int>;
template class Array<float>;
template class Array<double>;
template class Array<char>;

// Note no need for trailing underscore.
extern "C" {
Array<int> *IntArray__new(int dim, int *size) {
  return new Array<int>(dim, size);
}

Array<float> *FloatArray__new(int dim, int *size) {
  return new Array<float>(dim, size);
}

Array<double> *DoubleArray__new(int dim, int *size) {
  return new Array<double>(dim, size);
}

Array<char> *CharArray__new(int dim, int *size) {
  return new Array<char>(dim, size);
}

void IntArray__delete(void *array) {
  Array<int> *intArray = static_cast<Array<int> *>(array);
  delete intArray;
}

void FloatArray__delete(void *array) {
  Array<float> *floatArray = static_cast<Array<float> *>(array);
  delete floatArray;
}

void DoubleArray__delete(void *array) {
  Array<double> *doubleArray = static_cast<Array<double> *>(array);
  delete doubleArray;
}

void CharArray__delete(void *array) {
  Array<char> *charArray = static_cast<Array<char> *>(array);
  delete charArray;
}

void IntArray__set(void *array, const int dim, const int *indices,
                   const int value) {
  Array<int> *a = static_cast<Array<int> *>(array);
  a->set(indices, value);
}

void FloatArray__set(void *array, const int dim, const int *indices,
                     const float value) {
  Array<float> *a = static_cast<Array<float> *>(array);
  a->set(indices, value);
}

void DoubleArray__set(void *array, const int dim, const int *indices,
                      const double value) {
  Array<double> *a = static_cast<Array<double> *>(array);
  a->set(indices, value);
}

void CharArray__set(void *array, const int dim, const int *indices,
                    const char value) {
  Array<char> *a = static_cast<Array<char> *>(array);
  a->set(indices, value);
}

int IntArray__get(void *array, const int dim, const int *indices) {
  Array<int> *a = static_cast<Array<int> *>(array);
  return static_cast<int>(a->get(indices));
}

float FloatArray__get(void *array, const int dim, const int *indices) {
  Array<float> *a = static_cast<Array<float> *>(array);
  return static_cast<float>(a->get(indices));
}

double DoubleArray__get(void *array, const int dim, const int *indices) {
  Array<double> *a = static_cast<Array<double> *>(array);
  return static_cast<double>(a->get(indices));
}

char CharArray__get(void *array, const int dim, const int *indices) {
  Array<char> *a = static_cast<Array<char> *>(array);
  return static_cast<char>(a->get(indices));
}

void IntArray__print(void *array) {
  Array<int> *a = static_cast<Array<int> *>(array);
  a->print();
}

void FloatArray__print(void *array) {
  Array<float> *a = static_cast<Array<float> *>(array);
  a->print();
}

void DoubleArray__print(void *array) {
  Array<double> *a = static_cast<Array<double> *>(array);
  a->print();
}

void CharArray__print(void *array) {
  Array<char> *a = static_cast<Array<char> *>(array);
  a->print();
}

void IntArray__accumulate(void *array, const int dim, int *indices,
                          const int value) {
  Array<int> *a = static_cast<Array<int> *>(array);
  a->accumulate(indices, value);
}

void FloatArray__accumulate(void *array, const int dim, int *indices,
                            const float value) {
  Array<float> *a = static_cast<Array<float> *>(array);
  a->accumulate(indices, value);
}

void DoubleArray__accumulate(void *array, const int dim, int *indices,
                             const double value) {
  Array<double> *a = static_cast<Array<double> *>(array);
  a->accumulate(indices, value);
}

void IntArray__increment(void *array, const int dim, const int *indices) {
  Array<int> *a = static_cast<Array<int> *>(array);
  a->increment(indices);
}

void FloatArray__increment(void *array, const int dim, const int *indices) {
  Array<float> *a = static_cast<Array<float> *>(array);
  a->increment(indices);
}

void DoubleArray__increment(void *array, const int dim, const int *indices) {
  Array<double> *a = static_cast<Array<double> *>(array);
  a->increment(indices);
}
}
