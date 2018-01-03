#include "mex.h"
#include "math.h"
#include "memory.h"
#include <vector>
using namespace std;

struct PointI
{
	int x;
	int y;
	bool flag;
	PointI();
	PointI(int vx,int vy, bool vf=true);
};
PointI::PointI()
{
	 x=0,y=0,flag = true;
}

PointI::PointI(int vx,int vy, bool vf)
{
	x = vx;
	y=vy;
	flag = vf;
}

class ConnectRegion
{
public:
	double numlabel;
	vector<vector<PointI> > regions;
	double* image;
	ConnectRegion();
	ConnectRegion(double* temp,int rs, int cs);
	~ConnectRegion();
	void Compute();
private:
	void labeling(int x, int y,double value);
	int rows;
	int cols;
};
ConnectRegion::ConnectRegion()
{
	image = NULL;
	rows = 0;
	cols = 0;
	numlabel = 0;
	regions.clear();
}
ConnectRegion::ConnectRegion(double* temp, int rs, int cs)
{
	double idxmax;
	int i;
	rows = rs;
	cols = cs;
	image = (double*) mxMalloc(rs*cs*sizeof(double));
	memcpy(image,temp,rs*cs*sizeof(double));
	numlabel = 0;
	regions.clear();
        idxmax = image[0];
	for(i=0; i<rows*cols; i++)
		 if(idxmax < image[i]) idxmax = image[i];
        //printf("idxmax=%f\n",idxmax);
	for(i=0; i<rows*cols; i++)
		 image[i] -= idxmax;
	Compute();
}
ConnectRegion::~ConnectRegion()
{
	if (image)
	{
		mxFree(image);
	}
}

//label connected regions with the same number from 1 to numlabel
void ConnectRegion::labeling(int x,int y, double value)
{
	int nx,ny,x1,y1;
	double pvalue;
	int label;
	//bool inner=true;
	image[y*rows+x] = numlabel;
	label = (int)(numlabel-1);
	PointI ptemp(x,y);
	vector<PointI> vec_pt_open;
	vec_pt_open.clear();
	vec_pt_open.push_back(ptemp);
	//ptemp.flag = true;
	//ptemp.x = x;
	//ptemp.y = y;
	while(vec_pt_open.size()>0)
	{
		ptemp = vec_pt_open.back();
                ptemp.flag = true;
		vec_pt_open.pop_back();
		for(nx=-1;nx<2;nx++)
		{
			for(ny=-1;ny<2;ny++)
			{
				if(nx==0 && ny==0) continue;
				x1 = ptemp.x+nx;
				y1 = ptemp.y+ny;
				if(x1 >= 0 && x1<rows && y1>=0 && y1<cols)
				{
					if(image[y1*rows+x1] == value) 
					{
						PointI ptnext(x1,y1);
						image[y1*rows+x1] = numlabel;
						vec_pt_open.push_back(ptnext);
					}
					else if(image[y1*rows+x1] != numlabel)
						ptemp.flag = false;
				}
				else
					ptemp.flag = false;
			}
		}
		regions[label].push_back(ptemp);
	}
	
	
	/*if( numlabel == 0)
	{
		vector<PointI> vecptemp;
		PointI ptemp;
		numlabel++;
		ptemp.flag = true;
		ptemp.x = x;
		ptemp.y = y;
		vecptemp.push_back(ptemp);
		regions.push_back(vecptemp);
		image[y*rows+x] = numlabel;
		for(nx=-1;nx<2;nx++)
		{
			for(ny=-1;ny<2;ny++)
			{
				if(nx==0 && ny==0) continue;
				x1 = x+nx;
				y1 = y+ny;
				if(x1 >= 0 && x1<rows && y1>=0 && y1<cols)
				  labeling(x1,y1,pvalue,numlabel-1);
			}
		}
	}
	else
	{
		if(pvalue == value)
		{
			PointI ptemp;
			ptemp.flag = true;
			ptemp.x = x;
			ptemp.y = y;
			image[y*rows+x] = idxlabel;
			regions[idxlabel].push_back(ptemp);
			for(nx=-1;nx<2;nx++)
			{
				for(ny=-1;ny<2;ny++)
				{
					if(nx==0 && ny==0) continue;
					x1 = x+nx;
					y1 = y+ny;
					if(x1 >= 0 && x1<rows && y1>=0 && y1<cols)
					{
						if(image[y1*rows+x1] < 1) labeling(x1,y1,pvalue,idxlabel);
					}
				}
			}
		}
		else
		{
			vector<PointI> vecptemp;
			PointI ptemp;
			numlabel++;
			ptemp.flag = true;
			ptemp.x = x;
			ptemp.y = y;
			vecptemp.push_back(ptemp);
			regions.push_back(vecptemp);
			image[y*rows+x] = numlabel;
			for(nx=-1;nx<2;nx++)
			{
				for(ny=-1;ny<2;ny++)
				{
					if(nx==0 && ny==0) continue;
					x1 = x+nx;
					y1 = y+ny;
					if(x1 >= 0 && x1<rows && y1>=0 && y1<cols)
					labeling(x1,y1,pvalue,numlabel-1);
				}
			}
		}
	}*/
}

void ConnectRegion::Compute()
{
	int i,j;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			if(image[j*rows+i] < 1)
			{
				vector<PointI> vecptemp;
				vecptemp.clear();
				regions.push_back(vecptemp);
				numlabel++;
				labeling(i,j,image[j*rows+i]);
			}
		}
	}
}

int calWidth(vector<PointI>& vecbound, vector<PointI> vecpt,double* image, int rows, int cols)
{
	double* erodeflag;
	size_t i;
	int ni,nj,x1,y1;
	vector<PointI> contour;
	double label;
	int step = 0;
	contour.clear();
	erodeflag = (double*) mxMalloc(rows*cols*sizeof(double));
	memcpy(erodeflag,image,rows*cols*sizeof(double));
	//printf("beginning cal width!\n");
	for(i=0; i<vecpt.size(); i++)
	{
		if(!vecpt[i].flag) 
		{
			contour.push_back(vecpt[i]);
			/*if(vecpt[i].y*rows + vecpt[i].x > rows*cols)
			{
				printf("x=%d, y=%d rows=%d cols=%d\n",vecpt[i].x,vecpt[i].y);
				mexErrMsgTxt("over boundary");
			}*/
			erodeflag[vecpt[i].y*rows + vecpt[i].x] = 0;
		}
	}
	//printf("finishing cal width!\n");
	vecbound = contour;
	if(vecpt.size()<=0)
	{
		printf("Empty region!\n");
		mxFree(erodeflag);
		return 0;
	}

	label = image[vecpt[0].y*rows + vecpt[0].x];
	//printf("the first point of contour:%d %d total=%d label=%f",contour[0].x,contour[0].y,contour.size(),label);
	//if(label==1.0)
		//printf("%f %f %f %f",erodeflag[(contour[0].y-1)*rows+contour[0].x], erodeflag[(contour[0].y)*rows+contour[0].x-1],erodeflag[(contour[0].y+1)*rows+contour[0].x],erodeflag[(contour[0].y)*rows+contour[0].x+1]);
	while(contour.size()>0)
	{
		step ++;
		vector<PointI> contourtemp;
		contourtemp.clear();
		for(i=0; i<contour.size(); i++)
		{
			for(ni=-1;ni<2;ni++)
			{
				for(nj=-1;nj<2;nj++)
				{
					if(ni==0 && nj==0) continue;
					x1 = contour[i].x + ni;
					y1 = contour[i].y + nj;
					if(x1>=0 && x1<rows && y1>=0 && y1<cols)
					{
						//if(erodeflag[y1*rows+ x1] == 1.0)
						//		printf("%f %f %d\n", erodeflag[y1*rows + x1],label, (int)(erodeflag[y1*rows + x1]==label));
						if(erodeflag[y1*rows+ x1] == label) 
						{
							PointI pttemp;
							pttemp.x = x1;
							pttemp.y = y1;
							erodeflag[y1*rows+x1] = 0;
							contourtemp.push_back(pttemp);
						}
					}
				}
			}
		}
		contour = contourtemp;
	}
	mxFree(erodeflag);
	return step;
}
void getNeigbourRegion(vector<double>& neibour_label, vector<double>& neibour_idx, const vector<PointI> boundary,  const double* image_label,const double* image_idx, const int rows, const int cols)
{
	neibour_label.clear();
	neibour_idx.clear();
	if(boundary.size()<=0) {
		printf("Empty boundary input!");
		return;
	}
	double label = image_label[boundary[0].y*rows + boundary[0].x];
	for(int i=0; i<boundary.size(); i++)
	{
		int ni,nj,x1,y1;
		double plabel = image_label[boundary[i].y*rows + boundary[i].x];
		for(ni=-1;ni<2;ni++)
			for(nj=-1;nj<2;nj++)
			{
				if(ni == 0 && nj==0) continue;
				x1 = boundary[i].x + ni;
				y1 = boundary[i].y + nj;
				if(x1>=0&&x1<rows && y1>=0 && y1<cols)
				{
					double tlabel = image_label[y1*rows+x1];
					if(tlabel != plabel)
					{
						bool existflag = false;
						for(int j=0;j<neibour_label.size();j++){
							if(tlabel == neibour_label[j]) {existflag = true;break;}
						}
						if(!existflag) 
						{
							neibour_label.push_back(tlabel);
							//if(image_idx[y1*rows+x1] < 1.0) mexErrMsgTxt("index changed!");
							neibour_idx.push_back(image_idx[y1*rows+x1]);
						}
					}
				}
			}
	}
}
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )     
{
	 int i,j,rows,cols;
	 double* idxcenter;
	 //double idxmax = 0;
	 int contourstep;
	 if(nrhs < 1)
		 mexErrMsgTxt("1 inputs are required!");
	 if(!mxIsDouble(prhs[0])) mexErrMsgTxt("center index matrix is required to be double");

	 if(nlhs != 1) mexErrMsgTxt("1 output is allowed!");
	 /*get input data*/
	 idxcenter = (double*)mxGetData(prhs[0]); //index matrix of centers
	 rows = mxGetM(prhs[0]);
	 cols = mxGetN(prhs[0]);
	 if(nlhs != 1)
		 mexErrMsgTxt("1 output is allowed!");
	 
	 
	 //printf("Data initialization finished!\n");

	 ConnectRegion cr(idxcenter,rows,cols);
	 plhs[0] = mxCreateNumericMatrix(rows,cols,mxDOUBLE_CLASS,mxREAL);
        
	 double* output = (double*)mxGetData(plhs[0]);
         memset(output,rows*cols*sizeof(double),0);
         if(cr.regions.size()<=0)
         {
            printf("Wrong input images or computation");
         }
	 for(i=0;i<cr.regions.size();i++)
	 {
		 if(cr.regions[i].size()>0)
		 {
			 for(j=0;j<cr.regions[i].size();j++)
			 {
                                 if(!cr.regions[i][j].flag)
				     output[cr.regions[i][j].y*rows + cr.regions[i][j].x] = 1;
			 }
		 }
	 }
}
