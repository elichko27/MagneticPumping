//Ais.cpp
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <set>
#include "Ais.h"

double Ais(int x, int y) {

  double n = double(x); 
  double l = double(y); 

  if (y == 0){
    return -3./2.*n*(pow(n,2.) - 3.*n + 2.)/(4.*pow(n,2.) - 8.*n + 3.);
  } else if (y == 1) {
    return -3./2.*n*(12.*pow(n,3.) + 32.*pow(n,2.) + 12.*n + 3.)/(16.*pow(n,4.) - 40.*pow(n,2.) + 9.);
  } else if (y == 2) {
    return  3./2.*(8.*pow(n,6.) + 36.*pow(n,5.) - 50.*pow(n,4.) - 357.*pow(n,3.) - 336.*pow(n,2.) + 39.*n + 63.)/ 
      (32.*pow(n,5.) + 80.*pow(n,4.) - 80.*pow(n,3.) - 200.*pow(n,2.) + 18.*n + 45.);
  } else if (y == 3) {
    return 3./2.*(-64.*pow(n,4.) - 264.*pow(n,3.) -260.*pow(n,2.) + 18.*n)/
      (32.*pow(n,5.) + 80.*pow(n,4.) - 80.*pow(n,3.) -200.*pow(n,2.) + 18.*n + 45.);
  } else {
    std::cout << "Invalid number (l) in Ais - Ais only integers, 0-3 \n";
  }

}
