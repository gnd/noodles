// notes
// - fragcoord.xy / resolution - normalizes all coordinates between 0 and 1
//   it is fine when switching various resolutions, but not good for precision positioning

#ifdef GL_ES
precision mediump float;
#endif

//////////////////// GLOBALS /////////////////

uniform float time;
uniform vec2 resolution;
uniform vec2 rand;
uniform vec2 mouse;
vec2 p; //relative position
vec3 c; //color
float mw = 0.001; //minimal width
float asp; //aspect ratio

///////// PARAMS
/// m0 - m7 reserved for sound (8 bands)
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;
uniform float m5;
uniform float m6;
uniform float m7;


//////////////////// END OF GLOBALS ////////////////////////

// sin between 0 - 1
float zsin(in float a) {
	return (sin(a) + 1.0) / 2.0;
}

float a_sin(in float a) {
	return abs(sin(a));
}

// cos between 0 - 1
float zcos(in float a) {
        return (cos(a) + 1.0) / 2.0;
}

float a_cos(in float a) {
        return abs(cos(a));
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
        return ((p.y <= yp+width/2.0) && (p.y >= yp-width/2.0)) ? col : vec3(0.0);
	
}

// square
vec3 sqr(in float xp, in float yp, in float s, in vec3 col)
{
	float ys = s * asp;
	vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-s/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
	return cl;
}

// rect
vec3 rct(in float xp, in float yp, in float sx, in float sy, in vec3 col)
{
	float syy = sy*asp;
	vec3 cl = ((p.x > xp - sx/2.) && (p.x < xp + sx/2.) && (p.y > yp-syy/2.) && (p.y < yp + syy/2.)) ? col : vec3(0.0); 
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

//tst
float plot(vec2 p, float px){
  //return  smoothstep( pct-0.01, pct, p.y) - smoothstep( pct, pct+0.01, p.y);
	return smoothstep(px, px+1., p.y) + smoothstep(px, px-1., p.y);
}




void main( void ) {
////////////////////// INIT ////////////////////////
vec3 col;
// position upper left corner = 0,0
p = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
// aspect ratio
asp = 1.0;
vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cb = vec3(0.);

////////////////////// END OF INIT /////////////////
c=vec3(0.);

// SANDBOX
c+=smoothstep(0.,1.,p.x);
c+=crc(.5,.5,.1,cr);

// RED SHITE
int red = 1.;
if (red == 1) {
vec3 c_red;
	vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-0.5;
	vec4 d=p_ray; vec4 t; p_ray.z += time;
	for(float i=1.8; i>0.; i-=.1) {
		t = abs(mod(cos(p_ray),.3)-m1);
		float x = max(t.y*2., length(t.xz)*2.);
		c_red = vec3(i*.8,cos(i),i);
		if(x<.5) break;
		p_ray -= d*x;
	}
	for(float i=1.8; i>0.; i-=.1) {
		t = abs(mod(p_ray,.3)-0.34);
		float x = max(t.x*3., length(t.yz)*3.);
		c_red -= vec3(i*0.1,i,i);
		if(x<.5) break;
			p_ray -= d*x*time;
	}
c+=c_red;
}


// CHECK AUTONOM
int auto = 1;
if (auto == 1) {
	c+=plot(p, fract(p.x*sin(time/10.)));
	vec3 cz = vec3(zsin(time));
	vec3 ca = vec3(a_sin(time));
	c+=rct(0.,0.,cnt(5000),1.,ca);
	c+=rct(0.,1.,cnt(5000),1.,cz);
	c+=rct(1.,0.,cnt(5000),1.,cz);
	c+=rct(1.,1.,cnt(5000),1.,ca);
}

// CHECK RANDOMX/Y
int rnd = 0;
if (rnd == 1) {
	float once = cnt(1000);
	if (once > 0.9) {
		c+=rct(0.5,rand.x,1.-cnt(5000),.2,cw);
		c+=rct(0.5,rand.y,1.-cnt(5000),.2,cw);
	}
}

// CHECK MOUSE
//c+=sqr(0.5,0.5,mx.x,cw);

// CHECK NET
int cnet = 1;
if (cnet == 1) {
c+=sqr(0.075,.5,m0, cw);
c-=sqr(0.075,.5,m0-0.05, cw);
c+=sqr(0.2,.5,m1, cr);
c-=sqr(0.2,.5,m1-0.05, cw);
c+=sqr(0.325,.5,m2, cw);
c-=sqr(0.325,.5,m2-0.05, cw);
c+=sqr(0.45,.5,m3, cr);
c-=sqr(0.45,.5,m3-0.05, cw);
c+=sqr(0.575,.5,m4, cw);
c-=sqr(0.575,.5,m4-0.05, cw);
c+=sqr(0.7,.5,m5, cr);
c-=sqr(0.7,.5,m5-0.05, cw);
c+=sqr(0.825,.5,m6, cw);
c-=sqr(0.825,.5,m6-0.05, cw);
c+=sqr(0.95,.5,m7, cr);
c-=sqr(0.95,.5,m7-0.05, cw);
}

//c=vec3(cnt(50));
// DRAW
gl_FragColor = vec4(c, 1.0);

}

