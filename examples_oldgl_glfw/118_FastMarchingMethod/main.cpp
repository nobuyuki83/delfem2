/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer3.h"
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <set>

namespace dfm2 = delfem2;

const unsigned int ROW_SIZE = 200;
const unsigned int COL_SIZE = 200;
const double GRID_LENGTH = 0.01;

#define MAX_VALUE 1.0e38

class CPOINT{
public:
	unsigned int x_co = 0;
	unsigned int y_co = 0;
	double data = 0.0;
public:
	CPOINT() = default;
	CPOINT(int i1, int i2, double d){
		x_co = i1; y_co = i2; data = d;
	}
	bool operator<(const CPOINT& p)const{
		if( x_co == p.x_co && y_co == p.y_co )	return false;
		else if( data != p.data )	return data < p.data;
		else if( x_co != p.x_co )		return x_co < p.x_co;
		else if( y_co != p.y_co )		return y_co < p.y_co;
		return	false;
	}
};

void add_point_to_band(
    double phi_fnc[COL_SIZE+1][ROW_SIZE+1],
    std::set<CPOINT>& narrow_band,
    unsigned int i,
    unsigned int j,
    int direction)
{
	double tantative_value;
	double less_value = MAX_VALUE;
	double from_value = MAX_VALUE;

	switch(direction){
	case 0:
		if( j == 0 ){	less_value = phi_fnc[i][1]; }
		else	if( j == ROW_SIZE ){	less_value = phi_fnc[i][ROW_SIZE-1]; }
		else{	less_value = (phi_fnc[i][j+1] < phi_fnc[i][j-1]) ? phi_fnc[i][j+1] : phi_fnc[i][j-1]; }
		from_value = phi_fnc[i-1][j];
		break;
	case 1:
		if( i == 0 ){	less_value = phi_fnc[1][j]; }
		else	if( i == COL_SIZE ){	less_value = phi_fnc[COL_SIZE-1][j]; }
		else{	less_value = (phi_fnc[i+1][j] < phi_fnc[i-1][j]) ? phi_fnc[i+1][j] : phi_fnc[i-1][j]; }
		from_value = phi_fnc[i][j-1];
		break;
	case 2:
		if( j == 0 ){	less_value = phi_fnc[i][1]; }
		else	if( j == ROW_SIZE ){	less_value = phi_fnc[i][ROW_SIZE-1]; }
		else{ less_value = (phi_fnc[i][j+1] < phi_fnc[i][j-1]) ? phi_fnc[i][j+1] : phi_fnc[i][j-1]; }
		from_value = phi_fnc[i+1][j];
		break;
	case 3:
		if( i == 0 ){ less_value = phi_fnc[1][j]; }
		else	if( i == COL_SIZE ){	less_value = phi_fnc[COL_SIZE-1][j]; }
		else{	less_value = (phi_fnc[i+1][j] < phi_fnc[i-1][j]) ? phi_fnc[i+1][j] : phi_fnc[i-1][j]; }
		from_value = phi_fnc[i][j+1];
		break;
	default:
		less_value = MAX_VALUE;
		from_value = MAX_VALUE;
	}
	if( from_value + GRID_LENGTH >= less_value){
		tantative_value = ( from_value + less_value + sqrt(2*GRID_LENGTH*GRID_LENGTH-(less_value - from_value)*(less_value - from_value) ) ) / 2;
	}
	else{
		tantative_value = from_value + GRID_LENGTH;
	}
	phi_fnc[i][j] = tantative_value;
	narrow_band.insert( CPOINT(i, j, tantative_value) );
}

void set_point(
    double phi_fnc[COL_SIZE+1][ROW_SIZE+1],
    std::set<CPOINT>& narrow_band,
    bool is_active[COL_SIZE+1][ROW_SIZE+1],
    unsigned int i,
    unsigned int j)
{
	phi_fnc[i][j] = 0.0;
	is_active[i][j] = true;
	if( !is_active[i+1][j] )	add_point_to_band(phi_fnc,narrow_band,i+1,     j, 0);
	if( !is_active[i][j+1] )	add_point_to_band(phi_fnc,narrow_band,    i, j+1, 1);
	if( !is_active[i-1][j] )	add_point_to_band(phi_fnc,narrow_band,i-1,     j, 2);
	if( !is_active[i][j-1] )	add_point_to_band(phi_fnc,narrow_band,     i,j-1, 3);
}

void display(
    double phi_fnc[COL_SIZE+1][ROW_SIZE+1],
    double thres)
{
   glBegin(GL_LINES);
   for(unsigned i=0;i<COL_SIZE;i++){
	   for(unsigned j=0;j<ROW_SIZE;j++){
		   unsigned int flag = 0;
		   const double v00 = phi_fnc[i+0][j+0] - thres;
       const double v01 = phi_fnc[i+0][j+1] - thres;
       const double v10 = phi_fnc[i+1][j+0] - thres;
       const double v11 = phi_fnc[i+1][j+1] - thres;
		   if( v00 > 0 ){ flag += 1; }
		   if( v10 > 0 ){ flag += 2; }
		   if( v11 > 0 ){ flag += 4; }
		   if( v01 > 0 ){ flag += 8; }

		   const double rx0 = v00 / (v00 - v10);
		   const double ry0 = v00 / (v00 - v01);
		   const double ry1 = v10 / (v10 - v11);
		   const double rx1 = v01 / (v01 - v11);

		   switch( flag ){
		   case 0:	//0000
			   break;
		   case 1:	//0001
         glVertex2d( GRID_LENGTH*(i + rx0), GRID_LENGTH*j ); // 01
				 glVertex2d( GRID_LENGTH*i, GRID_LENGTH*(j + ry0) );
			   break;
		   case 2:  //0010
				glVertex2d( GRID_LENGTH*(i + rx0), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*(i + 1), GRID_LENGTH*(j + ry1) );
			   break;
		   case 3:	//0011
				glVertex2d( GRID_LENGTH*(i + 1), GRID_LENGTH*(j + ry1) );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*(j + ry0) );
			   break;
		   case 4:	//0100
				glVertex2d( GRID_LENGTH*(i + 1), GRID_LENGTH*(j + ry1) );
				glVertex2d( GRID_LENGTH*(i + rx1), GRID_LENGTH*(j+1) );
				break;
		   case 5:	//0101
				glVertex2d( GRID_LENGTH*(i + 1), GRID_LENGTH*(j + ry1) );
				glVertex2d( GRID_LENGTH*(i + rx1), GRID_LENGTH*(j + 1) );

				glVertex2d( GRID_LENGTH*(i + rx0), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*(j + ry0) );
				break;
		   case 6:	//0110
				glVertex2d( GRID_LENGTH*(i + rx0), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*(i + rx1), GRID_LENGTH*(j + 1) );
				break;
		   case 7:	//0111
				glVertex2d( GRID_LENGTH*(i + rx1), GRID_LENGTH*(j + 1) );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*(j + ry0) );
				break;
		   case 8:	//0001
				glVertex2d( GRID_LENGTH*(i + rx1), GRID_LENGTH*(j + 1) );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*(j + ry0) );
			    break;
		   case 9:	//1001
				glVertex2d( GRID_LENGTH*(i + rx0), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*(i + rx1), GRID_LENGTH*(j + 1) );
				break;
		   case 10:	//1010
				glVertex2d( GRID_LENGTH*(i + rx1), GRID_LENGTH*(j + 1) );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*(j + ry0) );
				glVertex2d( GRID_LENGTH*(i + rx0), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*(i + 1), GRID_LENGTH*(j + ry1) );
				break;
		   case 11:	//1011
				glVertex2d( GRID_LENGTH*(i + 1), GRID_LENGTH*(j + ry1) );
				glVertex2d( GRID_LENGTH*(i + rx1), GRID_LENGTH*(j + 1) );
				break;
		   case 12:	//1100
				glVertex2d( GRID_LENGTH*(i + 1), GRID_LENGTH*(j + ry1) );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*(j + ry0) );
			   break;
		   case 13:	//1101
				glVertex2d( GRID_LENGTH*(i + rx0), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*(i + 1), GRID_LENGTH*(j + ry1) );
				break;
		   case 14:	//1110
				glVertex2d( GRID_LENGTH*(i + rx0), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*(j + ry0) );
				break;
		   case 15:	//1111
				break;
		   default:
				break;
		   }
	   }
   }
   glEnd();

   glFlush();
}

int main(int argc, char *argv[])
{
  double phi_fnc[COL_SIZE+1][ROW_SIZE+1];
  bool is_active[COL_SIZE+1][ROW_SIZE+1];
  for(unsigned i=0;i<COL_SIZE+1;i++){
    for(unsigned j=0;j<ROW_SIZE+1;j++){
      phi_fnc[i][j] = MAX_VALUE;
      is_active[i][j] = false;
    }
  }
  {
    std::set<CPOINT> narrow_band;
    set_point(phi_fnc,narrow_band,is_active,10, 30);
    set_point(phi_fnc,narrow_band,is_active,30, 30);
    set_point(phi_fnc,narrow_band,is_active,10, 10);
    set_point(phi_fnc,narrow_band,is_active,30, 10);
    set_point(phi_fnc,narrow_band,is_active,20, 20);
    std::set<CPOINT>::iterator itr;
    while( !narrow_band.empty() ){
      itr = narrow_band.begin();
      unsigned int i = (*itr).x_co;
      unsigned int j = (*itr).y_co;
      is_active[i][j] = true;
      narrow_band.erase(itr);
      if( i != COL_SIZE && !is_active[i + 1][j])	add_point_to_band(phi_fnc, narrow_band, i + 1, j, 0);
      if( j != ROW_SIZE && !is_active[i][j + 1])	add_point_to_band(phi_fnc, narrow_band, i, j + 1, 1);
      if( i != 0        && !is_active[i - 1][j])	add_point_to_band(phi_fnc, narrow_band, i - 1, j, 2);
      if( j != 0        && !is_active[i][j - 1]) add_point_to_band(phi_fnc, narrow_band, i, j - 1, 3);
    }
  }

  dfm2::opengl::CViewer3 viewer;
  viewer.Init_oldGL();

  while( !glfwWindowShouldClose(viewer.window) ){
    viewer.DrawBegin_oldGL();
    ::glLineWidth(5);
    ::glColor3d(0,0,0);
    display(phi_fnc,0.05);
    ::glColor3d(1,0,0);
    display(phi_fnc,0.10);
    ::glColor3d(0,1,0);
    display(phi_fnc,0.15);
    ::glColor3d(0,0,1);
    display(phi_fnc,0.20);
    ::glColor3d(0,1,1);
    display(phi_fnc,0.25);
    ::glColor3d(1,0,1);
    display(phi_fnc,0.30);
    ::glColor3d(1,1,0);
    display(phi_fnc,0.35);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  return 0;
}
