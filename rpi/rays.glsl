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
        // implementation found at: lumina.sourceforge.net/Tutorials/Noise.html
        float k = sin(dot(co.xy ,vec2(12.9898,78.233)));
        return fract(k + k);
}

float noise2f( in vec2 p )
{
        vec2 ip = vec2(floor(p));
        vec2 u = fract(p);
        // http://www.iquilezles.org/www/articles/morenoise/morenoise.htm
        u = u*u;
        //u = u*u*u*((6.0*u-15.0)*u+10.0);

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
vec3 linex(in float xp, in float width, in vec3 col)
{
        return ((p.x <= xp+width/2.) && (p.x >= xp-width/2.)) ? col : vec3(0.0);
}

// horizontal line
vec3 liney(in float yp, in float width, in vec3 col)
{
        return ((p.y <= yp+width/2.) && (p.y >= yp-width/2.)) ? col : vec3(0.0);
}

// square
vec3 sqr(in float xp, in float yp, in float s, in vec3 col)
{
	float ys = s * asp;
	vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-s/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
	return cl;
}

// circle
vec3 crc(in float xp, in float yp, in float r, in vec3 col)
{
	float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
	vec3 cl = dist < r ? col : vec3(0.);
	return cl;
} 

// circle-test
vec3 crcx(in float xp, in float yp, in float r, in vec3 col)
{
        float dist = sqrt( (p.x*asp*tan(time) - xp) * (p.x*asp - xp) + (p.y*tan(time) - yp) * (p.y - yp*sin(time)) );
        vec3 cl = dist < r ? col : vec3(0.);
        return cl;
}


float sphere(vec3 ray, vec3 dir, vec3 center, float radius)
{
 vec3 rc = ray-center;
 float c = dot(rc, rc) - (radius*radius);
 float b = dot(dir, rc);
 float d = b*b - c;
 float t = -b - sqrt(abs(d));
 float st = step(0.0, min(t,d));
 return mix(-1.0, t, st);
}

vec3 background(float t, vec3 rd)
{
 vec3 light = normalize(vec3(sin(t), zsin(t*10.), cos(t)));
 float sun = max(0.0, dot(rd, light));
 float sky = max(0.0, dot(rd, vec3(0.0, 1.0, 0.0)));
 float ground = max(0.0, -dot(rd, vec3(0.0, 1.0, 0.0)));
 return 
  (pow(sun, 256.0)+0.2*pow(sun, 2.0))*vec3(2.0, 1.6, 1.0) +
  pow(ground, 0.5)*vec3(0.4, 0.3, 0.2) +
  pow(sky, 1.0)*vec3(0.5, 0.6, 0.7);
}




void main( void ) {
////////////////////// INIT ////////////////////////
vec3 col;
// position
p = ( gl_FragCoord.xy / resolution.xy );
// aspect ratio
asp = resolution.x / resolution.y;

////////////////////// END OF INIT /////////////////

c+=liney(0.5,cnt(1000), vec3(1.0));
//c+=linex(0.5,0.5,vec3(1.));

///rays
vec2 uv = (-1.0 + 2.0*p) *  vec2(asp, 1.0);
vec3 ro = vec3(0.0, 0.0, -3.);
vec3 rd = normalize(vec3(uv, 1.0));
vec3 p = vec3(0.0, 0.0, 0.0);
float t = sphere(ro, rd, p, 1.0);
vec3 nml = normalize(p - (ro+rd*t));
vec3 bgCol = background(time, rd);
rd = reflect(rd, nml);
vec3 rcol = background(time, rd) * vec3(0.9, 0.8, 1.0);

gl_FragColor = vec4(c, 1.0 );
gl_FragColor = vec4( mix(bgCol, rcol, step(0.0, t)), 1.0 );

}

