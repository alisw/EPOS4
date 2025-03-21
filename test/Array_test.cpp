#include "Array.hpp"

int main() {
  short int status = 1;
  
  int dim = 2;
  int shape[] = {4, 5};
  Array<int>* myIntArray = new Array<int>(dim, shape);

  int* array_shape = myIntArray->getShape();

  int number_of_elements = 1;  
  for(unsigned int i=0; i<dim; i++){
    number_of_elements *= array_shape[i];
  }
  
  if (number_of_elements == 20){
    status = 0;
  } else {
    status = 1;
  }
  
  delete myIntArray;
  return status;
}
