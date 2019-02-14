#ifdef GL_ES

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

float rand(vec2 co){
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
                                        mix(rand(ip), rand(ip+vec2(1.0,0.0)), u.x),
                                        mix(rand(ip+vec2(0.0,1.0)), rand(ip+vec2(1.0,1.0)), u.x),
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
	vec3 cl = ((p.x > xp) && (p.x < xp + s) && (p.y > yp) && (p.y < yp + ys)) ? col : vec3(0.0);
	return cl;
}

// circle
vec3 crc(in float xp, in float yp, in float r, in vec3 col)
{
	float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
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

//c+=liney(0.5, cnt(500)*0.1, vec3(1.0));
//c+=linex(0.5, cnt(1000), vec3(1.0,0.0,0.0));
//c+=liney(0.5, cnt(1000), vec3(0.0,1.0,0.0));
//c+=sqr(0.0,0.0, cnt(1000), vec3(0.0,0.0,1.0));
//c+=crc(asp/2.,0.5,0.45, vec3(1.,0.,0.));
//c+=sqr(0.1,0.5, 0.3, vec3(1.,0.,0.));

vec2 ps = vec2(p.x+time, p.y);
//c-=vec3(rand(ps));
c+=vec3(noise2f(ps));

gl_FragColor = vec4(c, 1.0 );

}

