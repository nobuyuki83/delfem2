/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/points.h"
#include "delfem2/mshuni.h"
#include "delfem2/slice.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/color.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <cstdlib>
#include <set>
#include <iostream>

namespace dfm2 = delfem2;

#define ROW_SIZE  200
#define COL_SIZE  200
#define GRID_LENGTH 0.01

#define G_WIDTH GRID_LENGTH * COL_SIZE
#define G_HEIGHT GRID_LENGTH * ROW_SIZE

#define MAX_VALUE 1.0e38

double phi_fnc[COL_SIZE+1][ROW_SIZE+1];
bool is_active[COL_SIZE+1][ROW_SIZE+1];

class CPOINT{
public:
	int x_co;
	int y_co;
	double data;
public:
	CPOINT(){}
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

std::set<CPOINT> narrow_band;

void init_field(void){
	int i, j;
	for(i=0;i<COL_SIZE+1;i++){
		for(j=0;j<ROW_SIZE+1;j++){
			phi_fnc[i][j] = MAX_VALUE;
			is_active[i][j] = false;
		}
	}
}

void add_point_to_band(int i, int j, int direction){
	double tantative_value;
	double less_value;
	double from_value;

	switch(direction){
	case 0:
		if( j == 0 )	less_value = phi_fnc[i][1];
		else	if( j == ROW_SIZE )	less_value = phi_fnc[i][ROW_SIZE-1];
		else	less_value = (phi_fnc[i][j+1] < phi_fnc[i][j-1]) ? phi_fnc[i][j+1] : phi_fnc[i][j-1]; 
		from_value = phi_fnc[i-1][j];
		break;
	case 1:
		if( i == 0 )	less_value = phi_fnc[1][j];
		else	if( i == COL_SIZE )	less_value = phi_fnc[COL_SIZE-1][j];
		else	less_value = (phi_fnc[i+1][j] < phi_fnc[i-1][j]) ? phi_fnc[i+1][j] : phi_fnc[i-1][j];
		from_value = phi_fnc[i][j-1];
		break;
	case 2:
		if( j == 0 )	less_value = phi_fnc[i][1];
		else	if( j == ROW_SIZE )	less_value = phi_fnc[i][ROW_SIZE-1];
		else	less_value = (phi_fnc[i][j+1] < phi_fnc[i][j-1]) ? phi_fnc[i][j+1] : phi_fnc[i][j-1];
		from_value = phi_fnc[i+1][j];
		break;
	case 3:
		if( i == 0 )	less_value = phi_fnc[1][j];
		else	if( i == COL_SIZE )	less_value = phi_fnc[COL_SIZE-1][j];
		else	less_value = (phi_fnc[i+1][j] < phi_fnc[i-1][j]) ? phi_fnc[i+1][j] : phi_fnc[i-1][j];
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

void set_point(int i, int j){
	phi_fnc[i][j] = 0.0;
	is_active[i][j] = true;
	
	if( is_active[i+1][j] == false )	add_point_to_band(i+1, j, 0);
	if( is_active[i][j+1] == false )	add_point_to_band(i, j+1, 1);
	if( is_active[i-1][j] == false )	add_point_to_band(i-1, j, 2);
	if( is_active[i][j-1] == false )	add_point_to_band(i, j-1, 3);
}

void DoFMM(void){
	std::set<CPOINT>::iterator itr;
	int i, j;

	while( narrow_band.size() != 0 ){
		itr = narrow_band.begin();
		i = (*itr).x_co;
		j = (*itr).y_co;
		is_active[i][j] = true;
		narrow_band.erase(itr);

		if( is_active[i+1][j] == false && i != COL_SIZE )	add_point_to_band(i+1, j, 0);
		if( is_active[i][j+1] == false && j != ROW_SIZE )	add_point_to_band(i, j+1, 1);
		if( is_active[i-1][j] == false && i != 0 )	add_point_to_band(i-1, j, 2);
		if( is_active[i][j-1] == false && j != 0 )	add_point_to_band(i, j-1, 3);
	}
}

void point_init_phi(void){
	init_field();
	set_point(30, 30);
	set_point(10, 30);
	set_point(30, 10);
	set_point(10, 10);
	set_point(20, 20);
	DoFMM();
}

void cone_init_phi(void){
	int i, j;
	for(i=0;i<COL_SIZE+1;i++){
		for(j=0;j<ROW_SIZE+1;j++){
			phi_fnc[i][j] = sqrt( double( (ROW_SIZE/2-j)*(ROW_SIZE/2-j)+(COL_SIZE/2-i)*(COL_SIZE/2-i) ) ) - 10.0;
		}
	}
}

void init_phi(void){
//	cone_init_phi();
	point_init_phi();
}

void up_phi(void){
	int i, j;
	for(i=0;i<COL_SIZE+1;i++){
		for(j=0;j<ROW_SIZE+1;j++){
			phi_fnc[i][j] += 0.005;
		}
	}
}

void down_phi(void){
	int i, j;
	for(i=0;i<COL_SIZE+1;i++){
		for(j=0;j<ROW_SIZE+1;j++){
			phi_fnc[i][j] -= 0.005;
		}
	}
}

void display(void)
{
	int i,j;
	int flag;

   glClear(GL_COLOR_BUFFER_BIT);
/*
   glColor3f( 0.0, 0.0, 0.0);
   glLineWidth(1);
   glBegin(GL_LINES);
   for(i=0;i<=ROW_SIZE;i++){
	  glVertex2d( G_WIDTH,  GRID_LENGTH*i);
	  glVertex2d( 0, GRID_LENGTH*i);
	}

   for(i=0;i<=COL_SIZE;i++){
	  glVertex2d( GRID_LENGTH * i, G_HEIGHT);
	  glVertex2d( GRID_LENGTH * i, 0);
   }
   glEnd();
*/
   glColor3d(1.0, 0.0, 0.0);
   glBegin(GL_LINES);
   for(i=0;i<COL_SIZE;i++){
	   for(j=0;j<ROW_SIZE;j++){
		   flag = 0;
		   if( phi_fnc[i][j] > 0 )	flag += 1;
		   if( phi_fnc[i+1][j] > 0 ) flag += 2;
		   if( phi_fnc[i+1][j+1] > 0 ) flag += 4;
		   if( phi_fnc[i][j+1] > 0 ) flag += 8;
		   
		   switch( flag ){
		   case 0:	//0000
			   break;
		   case 1:	//0001	
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i+1][j]), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*j + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i][j+1]) );
			   break;
		   case 2:  //0010
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i+1][j]), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i + GRID_LENGTH, GRID_LENGTH*j + phi_fnc[i+1][j]*GRID_LENGTH/(phi_fnc[i+1][j] - phi_fnc[i+1][j+1]) );
			   break;
		   case 3:	//0011
				glVertex2d( GRID_LENGTH*i + GRID_LENGTH, GRID_LENGTH*j + phi_fnc[i+1][j]*GRID_LENGTH/(phi_fnc[i+1][j] - phi_fnc[i+1][j+1]) );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*j + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i][j+1]) );
			   break;
		   case 4:	//0100
				glVertex2d( GRID_LENGTH*i + GRID_LENGTH, GRID_LENGTH*j + phi_fnc[i+1][j]*GRID_LENGTH/(phi_fnc[i+1][j] - phi_fnc[i+1][j+1]) );
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j+1]*GRID_LENGTH/(phi_fnc[i][j+1] - phi_fnc[i+1][j+1]) , GRID_LENGTH*j + GRID_LENGTH);
				break;
		   case 5:	//0101
				glVertex2d( GRID_LENGTH*i + GRID_LENGTH, GRID_LENGTH*j + phi_fnc[i+1][j]*GRID_LENGTH/(phi_fnc[i+1][j] - phi_fnc[i+1][j+1]) );
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j+1]*GRID_LENGTH/(phi_fnc[i][j+1] - phi_fnc[i+1][j+1]) , GRID_LENGTH*j + GRID_LENGTH);

				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i+1][j]), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*j + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i][j+1]) );
				break;
		   case 6:	//0110
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i+1][j]), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j+1]*GRID_LENGTH/(phi_fnc[i][j+1] - phi_fnc[i+1][j+1]), GRID_LENGTH*j + GRID_LENGTH);
				break;
		   case 7:	//0111
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j+1]*GRID_LENGTH/(phi_fnc[i][j+1] - phi_fnc[i+1][j+1]) , GRID_LENGTH*j + GRID_LENGTH);
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*j + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i][j+1]));
				break;
		   case 8:	//1000
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j+1]*GRID_LENGTH/(phi_fnc[i][j+1] - phi_fnc[i+1][j+1]) , GRID_LENGTH*j + GRID_LENGTH);
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*j + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i][j+1]));
			    break;
		   case 9:	//1001
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i+1][j]), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j+1]*GRID_LENGTH/(phi_fnc[i][j+1] - phi_fnc[i+1][j+1]), GRID_LENGTH*j + GRID_LENGTH);
				break;
		   case 10:	//1010
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j+1]*GRID_LENGTH/(phi_fnc[i][j+1] - phi_fnc[i+1][j+1]) , GRID_LENGTH*j + GRID_LENGTH);
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*j + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i][j+1]));
		
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i+1][j]), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i + GRID_LENGTH, GRID_LENGTH*j + phi_fnc[i+1][j]*GRID_LENGTH/(phi_fnc[i+1][j] - phi_fnc[i+1][j+1]) );
				break;
		   case 11:	//1011
				glVertex2d( GRID_LENGTH*i + GRID_LENGTH, GRID_LENGTH*j + phi_fnc[i+1][j]*GRID_LENGTH/(phi_fnc[i+1][j] - phi_fnc[i+1][j+1]) );
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j+1]*GRID_LENGTH/(phi_fnc[i][j+1] - phi_fnc[i+1][j+1]) , GRID_LENGTH*j + GRID_LENGTH);
				break;
		   case 12:	//1100
				glVertex2d( GRID_LENGTH*i + GRID_LENGTH, GRID_LENGTH*j + phi_fnc[i+1][j]*GRID_LENGTH/(phi_fnc[i+1][j] - phi_fnc[i+1][j+1]) );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*j + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i][j+1]) );
			   break;
		   case 13:	//1101
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i+1][j]), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i + GRID_LENGTH, GRID_LENGTH*j + phi_fnc[i+1][j]*GRID_LENGTH/(phi_fnc[i+1][j] - phi_fnc[i+1][j+1]) );
				break;
		   case 14:	//1110
				glVertex2d( GRID_LENGTH*i + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i+1][j]), GRID_LENGTH*j );
				glVertex2d( GRID_LENGTH*i, GRID_LENGTH*j + phi_fnc[i][j]*GRID_LENGTH/(phi_fnc[i][j] - phi_fnc[i][j+1]) );
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

/*
void init(void)
{
  glClearColor(1.0, 1.0, 1.0, 0.0);
  init_phi();
}

void resize(int w, int h)
{
  glViewport(0, 0, w, h);
  glLoadIdentity();
  glOrtho( -0.1 , G_WIDTH + 0.1, -0.1, G_HEIGHT + 0.1, -1.0, 1.0);
}

void keyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033':
    exit(0);
  case 'u':
	  up_phi();
	  glutPostRedisplay();
	  break;
  case 'd':
	  down_phi();
	  glutPostRedisplay();
	  break;	  
  default:
    break;
  }
}
*/

int main(int argc, char *argv[])
{
  init_phi();
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();

  while( !glfwWindowShouldClose(viewer.window) ){
    for(unsigned int itr=0;itr<100;++itr) {
      down_phi();
      //
      viewer.DrawBegin_oldGL();
      display();
      viewer.SwapBuffers();
      glfwPollEvents();
    }
    for(unsigned int itr=0;itr<100;++itr) {
      up_phi();
      //
      viewer.DrawBegin_oldGL();
      display();
      viewer.SwapBuffers();
      glfwPollEvents();
    }
  }
  return 0;
}
