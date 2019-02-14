#ifdef GL_ES

// TODO
// - test rand
// - test PD connection

// notes
// - fragcoord.xy / resolution - normalizes all coordinates between 0 and 1
//   it is fine when switching various resolutions, but not good for precision positioning

precision mediump float;
#endif
#define GX 6.0
#define GY 5.0

//////////////////// GLOBALS /////////////////

uniform float time;
uniform vec2 resolution;
// this should be rand
uniform float rand;
uniform float sa;
uniform float sb;
uniform float sc;
uniform float sd;
uniform float se;
uniform float sf;
uniform float sg;
uniform float sh;
uniform float si;
uniform float sj;

vec2 p; //relative position
vec3 c; //color
float mw = 0.001; //minimal width
float asp; //aspect ratio


//////////////////// END OF GLOBALS ////////////////////////

// sinusoid between 0 - 1
float zsin(in float a) {
	return (sin(a) + 1.0) / 2.0;
}

// cosinusoid between 0 - 1
float zcos(in float a) {
        return (cos(a) + 1.0) / 2.0;
}

// modulo
float mdl(in float a, in float b)
{
        return a - b * floor(a/b + 0.001);
}

float rnd(vec2 co){
        float k = sin(dot(co.xy ,vec2(12.9898,78.233)));
        return fract(k + k);
}

float noise2f( in vec2 p )
{
        vec2 ip = vec2(floor(p));
        vec2 u = fract(p);
        //u = u*u;
        u = u*u*u;

        float res = mix(
                                        mix(rnd(ip), rnd(ip+vec2(1.0,0.0)), u.x),
                                        mix(rnd(ip+vec2(0.0,1.0)), rnd(ip+vec2(1.0,1.0)), u.x),
                                        u.y);

        return res - 0.25;
        //return 2.0* (res-0.5);
}


// counter 0 to 1 
float cnt(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

// vertical line
vec3 olx(in float xp, in float width, in vec3 col)
{
        return ((p.x <= xp+width/2.) && (p.x >= xp-width/2.)) ? col : vec3(0.0);
}

// vertical line
vec3 lx(in float xp, in float width, in vec3 c, in vec3 col)
{
        return ((p.x <= xp+width/2.) && (p.x >= xp-width/2.)) ? -c+col : vec3(0.0);
}

// old horizontal line
vec3 oly(in float yp, in float width, in vec3 col)
{
        return ((p.y <= yp+width/2.) && (p.y >= yp-width/2.)) ? col : vec3(0.0);
}

// horizontal line
vec3 ly(in float yp, in float width, in vec3 c, in vec3 col)
{
        return ((p.y <= yp+width/2.) && (p.y >= yp-width/2.)) ? -c+col : vec3(0.0);
}

// old square
vec3 osqr(in float xp, in float yp, in float s, in vec3 col)
{
	float ys = s * asp;
	vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-s/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
	return cl;
}

// square
vec3 sqr(in float xp, in float yp, in float s, in vec3 c, in vec3 col)
{
        float ys = s * asp;
        vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-s/2.) && (p.y < yp + ys/2.)) ? -c+col : vec3(0.0);
        return cl;
}

// old circle
vec3 ocrc(in float xp, in float yp, in float r, in vec3 col)
{
        float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
        vec3 cl = dist < r ? col : vec3(0.);
        return cl;
}

// circle
vec3 crc(in float xp, in float yp, in float r, in vec3 c, in vec3 col)
{
	float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
	vec3 cl = dist < r ? -c+col : vec3(0.);
	return cl;
} 

// circle-test
vec3 crcx(in float xp, in float yp, in float r, in vec3 col)
{
        float dist = sqrt( (p.x*asp*tan(time) - xp) * (p.x*asp - xp) + (p.y*tan(time) - yp) * (p.y - yp*sin(time)) );
        vec3 cl = dist < r ? col : vec3(0.);
        return cl;
}




void main( void ) {
////////////////////// INIT ////////////////////////
vec3 col;
// position
p = ( gl_FragCoord.xy / resolution.xy );
// aspect ratio
asp = resolution.x / resolution.y;

////////////////////// END OF INIT /////////////////

// zilina styllie
vec3 c = vec3(0.5);
vec2 pp = p;
pp = 1. - 2.*pp;
float tt = 10.+sin(t)*50.;
float tx = cos(t);
//c+=noise2f(pow(log(vec2(sin(pp.x),sin(pp.y))),10.*t*log(pp)));
pp+= noise2f(pp*tt);
pp+= noise2f(pp*tt);
pp+= noise2f(pp*t);
c=c*c*c;

gl_FragColor = vec4(c, 1.0 );

}

