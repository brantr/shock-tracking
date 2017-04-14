#ifndef COMPUTE_SHOCK_PROPERTIES_H
#define COMPUTE_SHOCK_PROPERTIES_H
#include "shock_data_types.hpp"
struct shock_properties
{
  int   isnap;
  long  l;
  float d;
  long  id;
  long  ishock;
  float t;
  float x[3];
  float v[3];

  float dcut;
  float rcut; 
  long  ncut;
};
bool shock_properties_sort(shock_properties sa, shock_properties sb);
shock_properties compute_shock_properties(shock s, vector<tracer> t, float time, long ishock, int isnap);
void print_shock_properties(shock_properties sp);
void save_shock_properties(char fname[200], vector<shock_properties> sp);
#endif // COMPUTE_SHOCK_PROPERTIES_H