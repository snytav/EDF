#include "maxwell.h"
#include "load_data.h"
#include "particle.h"

#include <string>
#include <vector>

#define N 200

double f[N];

double get_maxv(std::vector<Particle> & p)
{
    double maxv = 0.0;
    
    for(int j = 0;j < p.size();j++)
    {
        double pu,pv,pw,ps,v;//,hv = 0.02/num;
	int n;
	
	pu = p[j].pu;
	pv = p[j].pv;
	pw = p[j].pw;
	
	ps = pu * pu + pv * pv + pw * pw;
	ps = pow(((pu * pu + pv * pv + pw * pw) + 1.0),-0.5);
	
	
	
	v = pu*ps;
	if(fabs(v) > maxv) maxv = fabs(v);
    }
    
    return maxv;
}

void getEDF(std::vector<Particle> & p,double *f,int num,double q,int nt)
{
    static int first  = 1;
    double maxv,hv;
//     e hv = 2*maxv/num;
    
    if(first == 1)
    {
       for(int i = 0;i < N;i++)
       {
	   f[i] = 0.0; 
	   
       }
      
       
       first  = 0;
    }
    maxv = get_maxv(p);
    hv = 2*maxv/num;
    
    FILE *f1;
    char fname[200];
    
    sprintf(fname,"electrons%05d_%10.3e.txt",nt,q);
    
    if((f1 = fopen(fname,"wt")) == NULL) return;
    
    
    
    for(int j = 0;j < p.size();j++)
    {
        double pu,pv,pw,ps,v;
	int n;
	
	pu = p[j].pu;
	pv = p[j].pv;
	pw = p[j].pw;
	
	ps = pu * pu + pv * pv + pw * pw;
	ps = pow(((pu * pu + pv * pv + pw * pw) + 1.0),-0.5);
	
	
	
	v = pu*ps;

	n = v/hv+num/2;
	f[n] += q;
	
	fprintf(f1,"%10d x %15.5e px %15.5e vx %15.5e n %5d \n",j,p[j].x,p[j].pu,v,n);
//     
    }
    fclose(f1);
    for(int i = 0;i < num;i++)
       {
	   f[i] /= p.size(); 
	   
       }
    sprintf(fname,"edf%05d.txt",nt,q);
    
    if((f1 = fopen(fname,"wt")) == NULL) return;  
    for(int i = 0;i < num;i++)
       {
	   fprintf(f1,"%5d %15.5e %15.5e\n",i,(i-num/2)*hv,f[i]); 
	   
       }   
       fclose(f1);
}


int main(int argc,char *argv[])
{
    double q = 1.0,maxv;  
    int    nt  = 128;//atoi(argv[1]);
    
  
    std::vector<Particle> ion_vp,el_vp,beam_vp;
    
    LoadParticleData(nt,ion_vp,el_vp,beam_vp,100,4,4);
    
    std::vector<Particle> all_el;
    all_el.reserve( el_vp.size() + beam_vp.size() );                // preallocate memory
    all_el.insert( all_el.end(), el_vp.begin(), el_vp.end() );        // add A;
    all_el.insert( all_el.end(), beam_vp.begin(), beam_vp.end() );        // add B;
    
    maxv = get_maxv(all_el);
    
    getEDF(all_el,f,N,q,nt);
//     q = 1.0/2000.0;
//     getEDF(beam_vp,f,N,q,nt);
//     FILE *f;
//     char fname[200];
//     
//     sprintf(fname,"electrons%05d.txt",1);
//     
//     if((f = fopen(fname,"wt")))
//     {
//        for(int j = 0;j < p.size();j++)
//        {
// 	   fprintf(f,"%10d x %15.5e px %15.5e vx %15.5e \n",j,el_vp[j].pu,);
//        }

    return 0;
}
