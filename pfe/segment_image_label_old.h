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
//#include "image.h"
//#include "misc.h"
//#include "filter.h"
#include "disjoint_set.h"
#include "segment_graph.h"
#include <vector>
#include <cstring>
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
static inline double diff(double *im, int nps, int nchs,
			int px1, int px2) {
	double dist = 0;
	int npx=0;
	for (int i = 0; i < nchs; i++,npx+=nps)
		dist += square(im[npx+px1]-im[npx+px2]);
	return sqrt(dist);
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
void segment_image(double *im, size_t* label, int height, int width,  int chs, int radius, double c, int min_size) {
 
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
    }
  }

  // segment
  universe *u = segment_graph(width*height, num, edges, c);
  
  // post process small components
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  delete [] edges;
  printf("%d\n", u->num_sets());
  int count=1;
  for (int y = 0; y < height; y++) {
    int idx=y;
    for (int x = 0; x < width; x++,idx+=height) {
      int comp = u->find(idx);
      if(label[comp] == 0)
      {
        label[comp] = count;
        count++;
      }
      label[idx] = label[comp]; 
    }
  } 
  printf("count=%d\n",count);
  //delete [] colors;  
  delete u;

  //return output;
}

void segment_image_fix1(double *im, size_t* label, int height, int width,  int chs, int radius, double c, int min_size,
        int fixsize) {
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
  //printf("number of regions: %d\n", vec_region.size());
  delete u;
  delete [] edges;
  int nregs =(int) vec_region.size();
  edges = new edge[nregs*(nregs-1)/2];
  num = 0;
  for(size_t i=0; i<nregs; i++)
  {
    for(size_t j=i+1; j<nregs; j++)
    {
      edges[num].a =(int) i;
      edges[num].b =(int) j;
      edges[num].w =(float) vec_region[i].dist(vec_region[j]);
      num++;
    }
  }
  //printf("number of edges: %d\n", num);
  //u = segment_graph(nregs, num, edges, c);
  u=new universe(nregs);
  //for(int i=0; i<count; i++)
  //  u.set_size(i,vec_region[i].get_size());
  sort(edges, edges + num);
  //int extra_regions = nregs - fixsize;
  //int num_erase = 0;
  for(int i=0; i<num && u->num_sets()>fixsize; i++ )
  {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if(a!=b)
    {
      u->join(a, b);
      //num_erase ++; 
    }
  }
  //printf("number of new segs: %d\n", u->num_sets());
  int count1=1;
  for (int i=0; i<nregs; i++){
      int comp = u->find(i);
      if(vec_region[comp].get_label() == 0)
      {
       vec_region[comp].set_label(count1);
        count1++;
      }
      vec_region[i].set_label((int) vec_region[comp].get_label());
    }
   //printf("label regrouped:%d!\n",count1);

  delete u;
  delete []edges;
  memset(label,0,npixels*sizeof(size_t));
  for (int i = 0; i < nregs; ++i)
  {
    vec_region[i].export2image(label,width,height);
  }

}

void segment_image_fix2(double *im, size_t* label, int height, int width,  int chs, int radius, double c, int min_size,
        int fixsize) {
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
  universe *u = segment_graph(npixels, num, edges, c);
  
  // post process small components
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  
  //printf("number of primary segmentation: %d\n", u->num_sets());
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
  //printf("number of regions: %d\n", vec_region.size());
  delete u;
  delete [] edges;
  int nregs =(int) vec_region.size();
  //edges = new edge[nregs*(nregs-1)/2];
  while(nregs > fixsize){
    double mindist=-1;
    size_t minidx_i=-1;
    size_t minidx_j=-1;
    for(size_t i=0; i<nregs; i++)
    {
      for(size_t j=i+1; j<nregs; j++)
      {
        //edges[num].a =(int) i;
        //edges[num].b =(int) j;
        //edges[num].w =(float) vec_region[i].dist(vec_region[j]);
        //num++;
        double temp = vec_region[i].dist(vec_region[j]);
        if(mindist < temp)
        {
          mindist =temp;
          minidx_i = i;
          minidx_j = j;
        }

      }
    }
    vec_region[minidx_i].merge(vec_region[minidx_j]);
    vec_region.erase(vec_region.begin()+minidx_j);
    nregs = (int)vec_region.size();
  }
  for(int i=0; i<nregs; i++)
  {
    vec_region[i].set_label(i+1);
  }
  //printf("number of edges: %d\n", num);
  /*u = segment_graph(nregs, num, edges, c);
  //for(int i=0; i<count; i++)
  //  u.set_size(i,vec_region[i].get_size());
  sort(edges, edges + num);
  int extra_regions = nregs - fixsize;
  int num_erase = 0;
  for(int i=0; i<num && num_erase<extra_regions;i++)
  {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if(a!=b)
    {
      u->join(a, b);
      num_erase ++; 
    }
  }
  //printf("number of new segs: %d\n", u->num_sets());
  int count1=1;
  for (int i=0; i<nregs; i++){
      int comp = u->find(i);
      if(vec_region[comp].get_label() == 0)
      {
       vec_region[comp].set_label(count1);
        count1++;
      }
      vec_region[i].set_label((int) vec_region[comp].get_label());
    }
   //printf("label regrouped:%d!\n",count1);

  delete u;
  delete []edges;*/
  memset(label,0,npixels*sizeof(size_t));
  for (int i = 0; i < nregs; ++i)
  {
    vec_region[i].export2image(label,width,height);
  }

}

void segment_image_fix(double *im, size_t* label, int height, int width,  int chs, int radius, double c, int min_size,
        int fixsize)
{
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
  }
  }

  // segment
  universe *u = segment_graph(width*height, num, edges, c);
  
  // post process small components
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if (a != b){
      if((u->size(a) < min_size) || (u->size(b) < min_size))
        u->join(a, b);
  }
  }

  int extra_region = u->num_sets()-fixsize;
  int count=0;
  for (int i = 0; i < num && count<=extra_region; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if (a != b){
        u->join(a, b);
        count++;
  }
  }
  delete [] edges;

  //hierarchical segmenation
  /*for(int i=0; i<depths; i++)
  {
      for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if (a != b){
      if((u->size(a) < min_size) || (u->size(b) < min_size))
        u->join(a, b);
  }
  }
}*/

  /**num_ccs = u->num_sets();*/

  count=1;
  for (int y = 0; y < height; y++) {
    int idx=y;
    for (int x = 0; x < width; x++,idx+=height) {
      int comp = u->find(idx);
      if(label[comp] == 0)
      {
        label[comp] = count;
        count++;
      }
      label[idx] = label[comp]; 
      //imRef(output, x, y) = colors[comp];
    }
  } 



  //printf("count=%d\n",count);
  //delete [] colors;  
  delete u;

}

#endif
