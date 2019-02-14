#ifdef GL_ES
precision highp float;
#endif
#define GX 6.0
#define GY 5.0

uniform float time;
uniform vec2 resolution;
vec2 position;
vec3 color;
uniform vec2 uv;
const int IMAX = int(GX * GY);
float a[IMAX];
vec4 rules;

// main
vec4 a1 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 a2 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 a3 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 a4 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 a5 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 a6 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 a7 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 a8 = vec4(0.0, 0.0, 0.0, 0.0);

// buffer
vec4 b1 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 b2 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 b3 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 b4 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 b5 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 b6 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 b7 = vec4(0.0, 0.0, 0.0, 0.0);
vec4 b8 = vec4(0.0, 0.0, 0.0, 0.0);

// initial
float fill() {
	a1 = vec4(1.0, 1.0, 1.0, 1.0);
        a2 = vec4(1.0, 1.0, 0.0, 1.0);
	a3 = vec4(1.0, 0.0, 1.0, 0.0);
        a4 = vec4(0.0, 1.0, 1.0, 1.0);
	a5 = vec4(1.0, 1.0, 1.0, 0.0);
        a6 = vec4(1.0, 1.0, 0.0, 1.0);
	a7 = vec4(1.0, 1.0, 1.0, 0.0);
        a8 = vec4(1.0, 1.0, 0.0, 1.0);
}

float moveab() {
	b1 = a1;
	b2 = a2;
	b3 = a3;
	b4 = a4;
	b5 = a5;
	b6 = a6;
	b7 = a7;
	b8 = a8;
}

float moveba() {
        a1 = b1;
        a2 = b2;
        a3 = b3;
        a4 = b4;
        a5 = b5;
        a6 = b6;
        a7 = b7;
        a8 = b8;
}


float mdl(in float a, in float b)
{
        return a - b * floor(a/b + 0.001);
}

float reada(in int ind) {

	float ret;

	if (ind == 0) {
	ret = a1.x;
}
if (ind == 1) {
	ret = a1.y;
}
if (ind == 2) {
	ret = a1.z;
}
if (ind == 3) {
	ret = a1.w;
}
if (ind == 4) {
	ret = a2.x;
}
if (ind == 5) {
	ret = a2.y;
}
if (ind == 6) {
	ret = a2.z;
}
if (ind == 7) {
	ret = a2.w;
}
if (ind == 8) {
	ret = a3.x;
}
if (ind == 9) {
	ret = a3.y;
}
if (ind == 10) {
	ret = a3.z;
}
if (ind == 11) {
	ret = a3.w;
}
if (ind == 12) {
	ret = a4.x;
}
if (ind ==13) {
	ret = a4.y;
}
if (ind == 14) {
	ret = a4.z;
}
if (ind == 15) {
	ret = a4.w;
}
if (ind == 16) {
	ret = a5.x;
}
if (ind == 17) {
	ret = a5.y;
}
if (ind == 18) {
	ret = a5.z;
}
if (ind == 19) {
	ret = a5.w;
}
if (ind == 20) {
	ret = a6.x;
}
if (ind == 21) {
	ret = a6.y;
}
if (ind == 22) {
	ret = a6.z;
}
if (ind == 23) {
	ret = a6.w;
}
if (ind == 24) {
	ret = a7.x;
}
if (ind == 25) {
	ret = a7.y;
}
if (ind == 26) {
	ret = a7.z;
}
if (ind == 27) {
	ret = a7.w;
}
if (ind == 28) {
	ret = a8.x;
}
if (ind == 29) {
	ret = a8.y;
}

	return ret;
}

float writea(in int ind,in float val) {
if (ind == 0) {
	a1.x = val;
}
if (ind == 1) {
	a1.y = val;
}
if (ind == 2) {
	a1.z = val;
}
if (ind == 3) {
	a1.w = val;
}
if (ind == 4) {
	a2.x = val;
}
if (ind == 5) {
	a2.y = val;
}
if (ind == 6) {
	a2.z = val;
}
if (ind == 7) {
	a2.w = val;
}
if (ind == 8) {
	a3.x = val;
}
if (ind == 9) {
	a3.y = val;
}
if (ind == 10) {
	a3.z = val;
}
if (ind == 11) {
	a3.w = val;
}
if (ind == 12) {
	a4.x = val;
}
if (ind ==13) {
	a4.y = val;
}
if (ind == 14) {
	a4.z = val;
}
if (ind == 15) {
	a4.w = val;
}
if (ind == 16) {
	a5.x = val;
}
if (ind == 17) {
	a5.y = val;
}
if (ind == 18) {
	a5.z = val;
}
if (ind == 19) {
	a5.w = val;
}
if (ind == 20) {
	a6.x = val;
}
if (ind == 21) {
	a6.y = val;
}
if (ind == 22) {
	a6.z = val;
}
if (ind == 23) {
	a6.w = val;
}
if (ind == 24) {
	a7.x = val;
}
if (ind == 25) {
	a7.y = val;
}
if (ind == 26) {
	a7.z = val;
}
if (ind == 27) {
	a7.w = val;
}
if (ind == 28) {
	a8.x = val;
}
if (ind == 29) {
	a8.y = val;
}

}
		
float writeb(in int ind, in float val) {

if (ind == 0) {
        b1.x = val;
}
if (ind == 1) {
        b1.y = val;
}
if (ind == 2) {
        b1.z = val;
}
if (ind == 3) {
        b1.w = val;
}
if (ind == 4) {
        b2.x = val;
}
if (ind == 5) {
        b2.y = val;
}
if (ind == 6) {
        b2.z = val;
}
if (ind == 7) {
        b2.w = val;
}
if (ind == 8) {
        b3.x = val;
}
if (ind == 9) {
        b3.y = val;
}
if (ind == 10) {
        b3.z = val;
}
if (ind == 11) {
        b3.w = val;
}
if (ind == 12) {
        b4.x = val;
}
if (ind ==13) {
        b4.y = val;
}
if (ind == 14) {
        b4.z = val;
}
if (ind == 15) {
        b4.w = val;
}
if (ind == 16) {
        b5.x = val;
}
if (ind == 17) {
        b5.y = val;
}
if (ind == 18) {
        b5.z = val;
}
if (ind == 19) {
        b5.w = val;
}
if (ind == 20) {
        b6.x = val;
}
if (ind == 21) {
        b6.y = val;
}
if (ind == 22) {
        b6.z = val;
}
if (ind == 23) {
        b6.w = val;
}
if (ind == 24) {
        b7.x = val;
}
if (ind == 25) {
        b7.y = val;
}
if (ind == 26) {
        b7.z = val;
}
if (ind == 27) {
        b7.w = val;
}
if (ind == 28) {
        b8.x = val;
}
if (ind == 29) {
        b8.y = val;
}

}

float counter(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

vec3 linex(in float xpos, in vec3 col)
{
	return int(mod(position.x*100.0,100.0)) == int(xpos*100.0) ? col : vec3(0.0);
}

vec3 liney(in float ypos, in vec3 col)
{
        return int(mod(position.y*100.0,100.0)) == int(ypos*100.0) ? col : vec3(0.0);
}

vec3 sqr(in float xpos, in float ypos, in float size, in vec3 col)
{
	float ysize = size * 1.0;
	vec3 cl = ((position.x > xpos) && (position.x < xpos + size) && (position.y > ypos) && (position.y < ypos + ysize)) ? col : vec3(0.0);
	return cl;
}


vec3 cell(in float ind)
{
	//parameters
	float XDST = 1.0/float(GX);
	float YDST = 1.0/float(GY);
	float SIZ = XDST - 0.01;
	
	//calculate position of cell in the grid
	float xp = mdl(ind,GX);
	float yp = floor((ind - xp) / GY);

	//this are the neighbors
        float aind = int(yp) == int(GY-1.0) ? (ind - 24.0) : (ind + 6.0);
        float bind = int(xp) == int(GX-1.0) ? (ind-5.0) : (ind+1.0);
        float cind = int(yp) == 0 ? (ind + 24.0) : (ind - 6.0);
        float dind = int(xp) == 0 ? (ind + 5.0) : (ind - 1.0);

	//determine color
	vec3 col = vec3(1.0,0.0,0.0);

	return sqr(xp*XDST,yp*YDST,SIZ,col);
	
}


vec3 rcell(in float ind)
{
        //parameters
        float XDST = 1.0/float(GX);
        float YDST = 1.0/float(GY);
        float SIZ = XDST - 0.01;

        //calculate position of cell in the grid
        float xp = mdl(ind,GX);
        float yp = floor((ind - xp) / GY);

        //determine color
        vec3 col = vec3(1.0);

        return sqr(xp*XDST,yp*YDST,SIZ,col);

}



vec3 dcell(in float ind)
{
        //parameters
        float XDST = 1.0/float(GX);
        float YDST = 1.0/float(GY);
        float SIZ = XDST - 0.01;

        //calculate position of cell in the grid
        float xp = mdl(ind,GX);
	float yp = floor((ind - xp) / GY);

	//this are the neighbors
        float aind = int(yp) == int(GY-1.0) ? (ind - 24.0) : (ind + 6.0);
	float bind = int(xp) == int(GX-1.0) ? (ind-5.0) : (ind+1.0);
	float cind = int(yp) == 0 ? (ind + 24.0) : (ind - 6.0); 
	float dind = int(xp) == 0 ? (ind + 5.0) : (ind - 1.0);

        //determine color
	vec3 col = vec3(1.0); //ignorant rule
	
        return sqr(xp*XDST,yp*YDST,SIZ,col);

}
	


void main( void ) {
vec3 col;
rules = vec4(1.0, 1.0, 1.0, 0.0);
position = ( gl_FragCoord.xy / resolution.xy );
int run = 1;

// fill up the array
fill();


//debug
col = vec3(1.0);
int kak = int(mod(30.0 * counter(10000),30.0));
color+=linex(float(kak)/100.0,col);
if (int(reada(kak)) == 1) {
	color += dcell(float(kak));
}


// MAIN LOGIC
if (run == 0) {
//	int i = int(mod(30.0 * counter(1000),30.0));
	int i = 0;
	float ind = float(i);

       	//calculate position of cell in the grid
       	float xp = mdl(ind,GX);
       	float yp = floor((ind - xp) / GY);

       	//this are the neighbors
	float aind = int(yp) == int(GY-1.0) ? (ind - 24.0) : (ind + 6.0);
       	float bind = int(xp) == int(GX-1.0) ? (ind-5.0) : (ind+1.0);
       	float cind = int(yp) == 0 ? (ind + 24.0) : (ind - 6.0);
       	float dind = int(xp) == 0 ? (ind + 5.0) : (ind - 1.0);
		
	//apply rule
	if ( (int(reada(int(aind))) == int(rules.x)) && (int(reada(int(bind))) == int(rules.y)) && (int(reada(int(cind))) == int(rules.z)) && (int(reada(int(dind))) == int(rules.w)) )
	{
		writeb(i,1.0);
	} else {
		writeb(i,0.0);
	}	

	//move to a
	moveba();

	// DRAWING
	for(int i = 0; i < IMAX; i++) {
		float ind = float(i);

		if (int(reada(i)) == 1) {

      	 	 	//parameters
			col = vec3(1.0);
       		 	float XDST = 1.0/float(GX);
       		 	float YDST = 1.0/float(GY);
        		float SIZ = XDST - 0.01;

        		//calculate position of cell in the grid
        		float xp = mdl(ind,GX);
        		float yp = floor((ind - xp) / GY);

			color += sqr(xp*XDST,yp*YDST,SIZ,col); 
		}
	}
}

gl_FragColor = vec4(color, 1.0 ); 

}

