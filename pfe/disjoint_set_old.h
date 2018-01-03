/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef __DISJOINT_SET_H__
#define __DISJOINT_SET_H__
#include <vector>
#include <cstring>
#include <cstdint>
using namespace std;
// disjoint-set forests using union-by-rank and path compression (sort of).
typedef struct {
  float w;
  int a, b;
} edge;

typedef struct {
  int rank;
  int p;
  int size;
} uni_elt;

class universe {
public:
  universe(int elements);
  ~universe();
  int find(int x);  
  void join(int x, int y);
  int size(int x) const { return elts[x].size; }
  int set_size(int x, int num){ elts[x].size = num; }
  int num_sets() const { return num; }

private:
  uni_elt *elts;
  int num;
};

universe::universe(int elements) {
  elts = new uni_elt[elements];
  num = elements;
  for (int i = 0; i < elements; i++) {
    elts[i].rank = 0;
    elts[i].size = 1;
    elts[i].p = i;
  }
}
  
universe::~universe() {
  delete [] elts;
}

int universe::find(int x) {
  int y = x;
  while (y != elts[y].p)
    y = elts[y].p;
  elts[x].p = y;
  return y;
}

void universe::join(int x, int y) {
  if (elts[x].rank > elts[y].rank) {
    elts[y].p = x;
    elts[x].size += elts[y].size;
  } else {
    elts[x].p = y;
    elts[y].size += elts[x].size;
    if (elts[x].rank == elts[y].rank)
      elts[y].rank++;
  }
  num--;
}
struct Point
{
  int x;
  int y;
};
class region
{
private:
  vector<Point> pts;
  vector<double> values;
  size_t label;
public:
  void addpoint(int x, int y);
  void addpixel(int x, int y,double* image, int width, int height,int chs);
  bool finish();
  bool export2image(size_t* image,int width,int height);
  void set_label(int n);
  size_t get_label();
  void get_values(vector<double>& v);
  void get_pts(vector<Point>& ptstemp);
  size_t get_size() {return pts.size();}
  double dist(region& r);
  bool merge(region& r);
  region();
};
 region::region(){
  pts.clear();
  values.clear();
  label = 0;
 };
 void region::addpoint(int x, int y)
 {
    Point pt;
    pt.x = x;
    pt.y = y;
    pts.push_back(pt);
 }

  void region::addpixel(int x, int y,double* image, int width, int height,int chs)
 {
    int npixels = width*height;
    int idx = x*height+y;
    if(pts.size()==0)
    {
      for(int i=0;i<chs;i++,idx+=npixels)
      {
        values.push_back(image[idx]);
      }
    }
    else
    {
      for(int i=0;i<chs;i++,idx+=npixels)
      {
        values[i] += image[idx];
      }
    }
    Point pt;
    pt.x = x;
    pt.y = y;
    pts.push_back(pt);
 }

bool region::finish()
{
  size_t n = pts.size();
  if( n <= 0)
    return false;
  for(size_t i=0; i<values.size(); i++)
    values[i] /= n;
  return true;
}

bool region::export2image(size_t* image,int width,int height)
{
  for(size_t i=0;i<pts.size();i++)
  {
    if(pts[i].x>=width || pts[i].y>=height) return false;
    image[pts[i].x*height+pts[i].y] = label;
  }
  return true;
}

void region::set_label(int n)
{
  label = n;
}

size_t region::get_label()
{
  return label;
}

void region::get_values(vector<double>& v)
{
  v=values;
}
void region::get_pts(vector<Point>& ptstemp)
{
  ptstemp=pts;
}

double region::dist(region& r)
{
  vector<double> v;
  r.get_values(v);
  if (v.size() != values.size())
  {
    return -1;
  }
  else
  {
    double dist = 0;
    for(size_t i=0; i<v.size();i++)
    {
      dist += square(v[i]-values[i]);
    }
    return dist;
  }
}

bool region::merge(region& r)
{
  vector<double> values1;
  r.get_values(values1);
  if (values1.size() != values.size())
  {
    return false;
  }
  else
  {
    vector<Point> pts1;
    r.get_pts(pts1);
    size_t num1 = pts1.size();
    size_t num2 = pts.size();
    pts.insert(pts.end(), pts1.begin(),pts1.end());
    if(pts.size()>0){
      for(size_t i=0; i<values.size(); i++)
      {
        values[i] = (values[i]*num2 + values1[i]*num1)/pts.size();
      }
    }
    return true;
  }
}

class cluster 
{
  private:
  vector<int> pts;
  double* values;
  int dims;
  size_t label;
  public:
  void addpoint(int id){pts.push_back(id);}
  bool addpoint(int id,const double* data, int ndata,int nchs)
  {
    if(nchs!=dims || id>= ndata)
    {  
      printf("wrong point added to cluster!\n");
      return false;
    }
    int idxtp = id*nchs;
    for(size_t i=0; i<nchs; i++)
      values[i] += data[i+idxtp];
    pts.push_back(id);
    return true;
  }

  bool finish()
  {
    size_t n = pts.size();
    if( n <= 0)
      return false;
    for(size_t i=0; i<dims; i++)
      values[i] /= n;
  return true;
  }
  //bool export2image(size_t* image,int width,int height);
  
  void set_dims(int val)
  {
    if(values)
      delete[] values;
    dims=val;
    values = new double[dims]; 
    memset(values,0,dims*sizeof(double));
  }
  int get_dims(){return dims;}

  void set_label(int n){label=n;}
  size_t get_label(){return label;}

  void get_values(double** v) { *v = const_cast<double*> (values); }
  void get_pts(vector<int>* ptstemp) { *ptstemp = pts; }
  size_t get_size() { return pts.size();}
  double dist(cluster* r) 
  {
    if(r->get_dims() != dims )
      return -1;
    double temp=0;
    double* vtemp=NULL;
    r->get_values(&vtemp);

    for(int i=0;i<dims; i++)
    {
		//printf("%f %f ",values[i],vtemp[i]);
      temp+= square(values[i]-vtemp[i]);
    }
    return temp;
  }

  bool merge(cluster* r)
  {
    if(r->get_dims() != dims )
      return false;
    double* valuestp=NULL;
    r->get_values(&valuestp);

    vector<int> ptstemp;
    r->get_pts(&ptstemp);

    size_t numtemp = ptstemp.size();
    size_t num = pts.size();
    pts.insert(pts.end(), ptstemp.begin(),ptstemp.end());
    if(pts.size()>0){
      for(size_t i=0; i<dims; i++)
        values[i] = (values[i]*numtemp + valuestp[i]*num)/pts.size();
    }
    return true;
  }

  bool export2label(uint32_t* output,int ndata)
  {
    for(size_t i=0;i<pts.size();i++)
    {
      if(pts[i] > ndata) return false;
      output[pts[i]] = label;
    }
    return true;
  }
  cluster(){values=NULL;pts.clear();dims=0;label=0;}
  cluster(int val)
  {
    values = NULL;
    dims = 0;
    label=0;
    pts.clear();
    set_dims(val);
  }
  void clear()
  {
	  delete[] values;
	  pts.clear();
  }
  /*~cluster()
  {
	  delete[] values;
	  pts.clear();
  }*/
};
#endif
