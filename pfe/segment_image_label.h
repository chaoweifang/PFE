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

#ifndef SEGMENT_IMAGE
#define SEGMENT_IMAGE

#include <cstdlib>
#include <time.h>
#include <cstdint>
//#include "image.h"
//#include "misc.h"
//#include "filter.h"
#include "disjoint_set.h"
#include "segment_graph.h"
#include <vector>
#include <cstring>
#include <list>
using namespace std;
// random color
rgb random_rgb(){ 
  rgb c;
  //double r;
 
  c.r = (uchar)(rand()%256);
  c.g = (uchar)(rand()%256);
  c.b = (uchar)(rand()%256);

  return c;
}

// dissimilarity measure between pixels
static inline double diff(const double *im, const int nps,const int nchs,
			const int px1,const int px2) {
	double dist = 0;
	int npx=0;
	for (int i = 0; i < nchs; i++,npx+=nps)
		dist += square(im[npx+px1]-im[npx+px2]);
	return sqrt(dist);
}

bool compareEdge(edge e1, edge e2)
{
  return e1.w<e2.w;
}
/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 * num_ccs: number of connected components in the segmentation.
 */
void segment_image_merge(double *im, uint32_t* label, int height,int width,int chs,int radius,double c,int min_size, int fixsize) {
  //printf("entering!\n");
  int npixels = width*height;
  // build graph
  int neighborsize = ((2 * radius + 1)*(2 * radius + 1) - 1) / 2;
  edge *edges = new edge[width*height*neighborsize];
  int num = 0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
    for (int j = -radius; j <= radius; j++)
    {
      for (int i = 0; i <= radius; i++)
      {
        if (i > 0 || (i==0&&j>0))
        {
          int nx = x + i;
          int ny = y + j;
          if (nx < 0 || nx >= width || ny < 0 || ny >= height)
            continue;
          edges[num].a = x * height + y;
          edges[num].b = nx * height + ny;
          edges[num].w = (float) diff(im, npixels, chs, edges[num].a, edges[num].b);
          num++;
        }
      }
    }
  }}

  //printf("graph built: %d\n", num);
  // segment
  universe *u = segment_graph_limit(npixels, num, edges, c,fixsize);
  //printf("number of primary segmentation: %d\n", u->num_sets());
  // post process small components
  for (int i = 0; i < num && u->num_sets()>fixsize; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b)  && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  
  //printf("number of erase small segmentation: %d\n", u->num_sets());
  vector<region> vec_region;

  int count=1;
  for (int y = 0; y < height; y++) {
    int idx=y;
    for (int x = 0; x < width; x++,idx+=height) {
      int comp = u->find(idx);
      if(comp >= npixels) printf("out of bound comp!\n");
      if(label[comp] == 0)
      {
        label[comp] = count;
        region r;
        r.addpixel(x,y,im,width,height,chs);
        vec_region.push_back(r);
        count++;
      }
      else
      {
        vec_region[label[comp]-1].addpixel(x,y,im,width,height,chs);
      }

      label[idx] = label[comp]; 
    }
  } 
  for(size_t i=0; i<vec_region.size(); i++)
    vec_region[i].finish();
  //printf("number of regions: %d\n", vec_region.size());
  delete u;
  list<region> li_region(vec_region.begin(),vec_region.end());
  int n_region = (int)vec_region.size();
  vec_region.clear();
  while(n_region>fixsize)
  {
    //get the smallest region
    list<region>::iterator iter_li=li_region.begin(), smallest_it=iter_li;
    int smallest=(int) smallest_it->get_size();
    for(; iter_li!=li_region.end(); iter_li++)
    {
      if(smallest>=(int)iter_li->get_size())
      {
        smallest =(int) iter_li->get_size();
        smallest_it = iter_li;
      }
    }

    //merge the smallest region to the most similar
    list<region>::iterator similar_it;
    iter_li=li_region.begin();
    double similar_dst =100000000.0;
    for(;iter_li!=li_region.end(); iter_li++)
    {
      if(iter_li==smallest_it) continue;
      double dsttmp = iter_li->dist(*smallest_it);
      if(similar_dst > dsttmp)
      {
        similar_dst = dsttmp;
        similar_it = iter_li;
      }
    }

    similar_it->merge(*smallest_it);
    li_region.erase(smallest_it);
    n_region--;
  }

   //int nregs = vec_region.size();
  memset(label,0,npixels*sizeof(uint32_t));
  list<region>::iterator iter_li=li_region.begin();
  count=1;
  for (; iter_li!=li_region.end(); ++iter_li)
  {
    iter_li->set_label(count++);
    iter_li->export2image(label,width,height);
  }
}

void segmentintr_image_merge(double *im, double * intr_info, uint32_t* label, int height,int width,int chs,int radius,double c,int min_size, int fixsize) {
  //printf("entering!\n");
  int npixels = width*height;
  // build graph
  int neighborsize = ((2 * radius + 1)*(2 * radius + 1) - 1) / 2;
  edge *edges = new edge[width*height*neighborsize];
  int num = 0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
    for (int j = -radius; j <= radius; j++)
    {
      for (int i = 0; i <= radius; i++)
      {
        if (i > 0 || (i==0&&j>0))
        {
          int nx = x + i;
          int ny = y + j;
          if (nx < 0 || nx >= width || ny < 0 || ny >= height)
            continue;
          edges[num].a = x * height + y;
          edges[num].b = nx * height + ny;
          edges[num].w = (float) diff(im, npixels, chs, edges[num].a, edges[num].b);
          num++;
        }
      }
    }
  }}

  //printf("graph built: %d\n", num);
  // segment
  universe *u = segment_graph_limit(npixels, num, edges, c,fixsize);
  //printf("number of primary segmentation: %d\n", u->num_sets());
  // post process small components
  for (int i = 0; i < num && u->num_sets()>fixsize; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b)  && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  
  //printf("number of erase small segmentation: %d\n", u->num_sets());
  vector<region> vec_region;

  int count=1;
  for (int y = 0; y < height; y++) {
    int idx=y;
    for (int x = 0; x < width; x++,idx+=height) {
      int comp = u->find(idx);
      if(comp >= npixels) printf("out of bound comp!\n");
      if(label[comp] == 0)
      {
        label[comp] = count;
        region r;
        r.addpixel(x,y,im,width,height,chs);
        vec_region.push_back(r);
        count++;
      }
      else
      {
        vec_region[label[comp]-1].addpixel(x,y,im,width,height,chs);
      }

      label[idx] = label[comp]; 
    }
  } 
  for(size_t i=0; i<vec_region.size(); i++)
    vec_region[i].finish();
  //printf("number of regions: %d\n", vec_region.size());
  delete u;
  list<region> li_region(vec_region.begin(),vec_region.end());
  int n_region = (int)vec_region.size();
  vec_region.clear();
  while(n_region>fixsize)
  {
    //get the smallest region
    list<region>::iterator iter_li=li_region.begin(), smallest_it=iter_li;
    int smallest=(int) smallest_it->get_size();
    for(; iter_li!=li_region.end(); iter_li++)
    {
      if(smallest>=(int)iter_li->get_size())
      {
        smallest =(int) iter_li->get_size();
        smallest_it = iter_li;
      }
    }

    //merge the smallest region to the most similar
    list<region>::iterator similar_it;
    iter_li=li_region.begin();
    double similar_dst =100000000.0;
    for(;iter_li!=li_region.end(); iter_li++)
    {
      if(iter_li==smallest_it) continue;
      double dsttmp = iter_li->dist(*smallest_it);
      if(similar_dst > dsttmp)
      {
        similar_dst = dsttmp;
        similar_it = iter_li;
      }
    }

    similar_it->merge(*smallest_it);
    li_region.erase(smallest_it);
    n_region--;
  }

   //int nregs = vec_region.size();
  memset(label,0,npixels*sizeof(uint32_t));
  list<region>::iterator iter_li=li_region.begin();
  count=1;
  for (; iter_li!=li_region.end(); ++iter_li)
  {
    iter_li->set_label(count++);
    iter_li->export2image(label,width,height);
  }
}

#endif
