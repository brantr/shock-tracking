#include "compute_shock_properties.hpp"

bool shock_properties_sort(shock_properties sa, shock_properties sb)
{
  //sort by increasing time
  if(sa.t==sb.t)
    return sa.isnap < sb.isnap;
  return sa.t < sb.t;
}
shock_properties compute_shock_properties(shock s, vector<tracer> t, float time, long ishock, int isnap)
{
  shock_properties sp;

  double dtot = 0;

  vector<tracer> tcut;

  //simple properties
  sp.l     = s.l;
  sp.d     = 0;
  sp.id    = s.id;
  sp.t     = time;
  sp.isnap = isnap;
  sp.rcut  = 0;


  //reset position and velocity
  for(int k=0;k<3;k++)
  {
    sp.x[k] = 0;
    sp.v[k] = 0;
  }


  //all particles
  //x       = 0.341561  0.736938  0.275418
  //vx      = -0.981087 -0.709419 -0.780719

  //dcut = 0.9 dmax


  //compute the density-weighted position and velocity

  sp.dcut = 0.9*t[0].d;
  sp.ncut = 0;
  double r;
  for(long tt=0;tt<t.size();tt++)
  {
    if(t[tt].d>sp.dcut)
    {
      if(t[tt].d>sp.d)
        sp.d = t[tt].d;
      tcut.push_back(t[tt]);
      sp.ncut++;
      dtot += t[tt].d;
      for(int k=0;k<3;k++)
      {
        sp.x[k] += t[tt].d*t[tt].x[k];
        sp.v[k] += t[tt].d*t[tt].v[k];
      }
    }
  }
  for(int k=0;k<3;k++)
  {
    sp.x[k] /= dtot;
    sp.v[k] /= dtot;
  }

  //find the extent of particles
  //used in computing the density
  //weighted properties
  for(long tt=0;tt<tcut.size();tt++)
  {
    r = 0;
    for(int k=0;k<3;k++)
      r += pow(tcut[tt].x[k] - sp.x[k], 2);
    r = sqrt(r);

    if(r>sp.rcut)
      sp.rcut = r;
  }

  //return the shock properties
  return sp;
}

void print_shock_properties(shock_properties sp)
{
  printf("Shock properties:\n");
  printf("Snapshot  = %04d\n",sp.isnap);
  printf("Time      = %f\n",sp.t);
  printf("Shock IDX = %ld\n",sp.ishock);
  printf("Length    = %ld\n",sp.l);
  printf("ID        = %ld\n",sp.id);
  printf("Density   = %f\n",sp.d);
  printf("dcut      = %f\n",sp.dcut);
  printf("rcut      = %f\n",sp.rcut);
  printf("ncut      = %ld\n",sp.ncut);
  printf("x         = % f\t% f\t% f\n",sp.x[0],sp.x[1],sp.x[2]);
  printf("v         = % f\t% f\t% f\n",sp.v[0],sp.v[1],sp.v[2]);

}

void save_shock_properties(char fname[200], vector<shock_properties> sp)
{

  std::sort(sp.begin(),sp.end(),shock_properties_sort);
  FILE *fp = fopen(fname,"w");
  for(long ss=0;ss<sp.size();ss++)
    fprintf(fp,"%04d\t%f\t%ld\t%ld\t%ld\t%f\t%f\t%f\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\n",sp[ss].isnap,sp[ss].t,sp[ss].ishock,sp[ss].l,sp[ss].id,sp[ss].d,sp[ss].dcut,sp[ss].rcut,sp[ss].ncut,sp[ss].x[0],sp[ss].x[1],sp[ss].x[2],sp[ss].v[0],sp[ss].v[1],sp[ss].v[2]);
  fclose(fp);
  for(long ss=0;ss<sp.size();ss++)
    printf("%04d\t%f\t%f\t%f\t%f\t%f\t% f\t% f\t% f\n",sp[ss].isnap,sp[ss].t,sp[ss].d,sp[ss].x[0],sp[ss].x[1],sp[ss].x[2],sp[ss].v[0],sp[ss].v[1],sp[ss].v[2]);
}