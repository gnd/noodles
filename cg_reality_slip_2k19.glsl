// TODO
// - test rand
// - test PD connection

// notes
// - fragcoord.xy / resolution - normalizes all coordinates between 0 and 1
//   it is fine when switching various resolutions, but not good for precision positioning

#ifdef GL_ES
precision mediump float;
#endif

//////////////////// GLOBALS /////////////////

uniform float time;
uniform vec2 resolution;
uniform sampler2D backbuffer;
uniform vec2 rand;
uniform vec2 mouse;
vec2 p,ppp; //relative position
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

/////////////////// MIDIMIX //////////////////////
uniform float m8,m9,m10;
uniform float m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30;
uniform float m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50;
uniform float m51,m52,m53,m54,m55,m56,m57,m58,m59;
float mx11 = m8/127.;
float mx12 = m9/127.;
float mx13 = m10/127.;
float mx14 = m11/127.;
float mx15 = m12/127.;
float mx16 = m13/127.;
float mx21 = m14/127.;
float mx22 = m15/127.;
float mx23 = m16/127.;
float mx24 = m17/127.;
float mx25 = m18/127.;
float mx26 = m19/127.;
float mx31 = m20/127.;
float mx32 = m21/127.;
float mx33 = m22/127.;
float mx34 = m23/127.;
float mx35 = m24/127.;
float mx36 = m25/127.;
float mx41 = m26/127.;
float mx42 = m27/127.;
float mx43 = m28/127.;
float mx44 = m29/127.;
float mx45 = m30/127.;
float mx46 = m31/127.;
float mx51 = m32/127.;
float mx52 = m33/127.;
float mx53 = m34/127.;
float mx54 = m35/127.;
float mx55 = m36/127.;
float mx56 = m37/127.;
float mx61 = m38/127.;
float mx62 = m39/127.;
float mx63 = m40/127.;
float mx64 = m41/127.;
float mx65 = m42/127.;
float mx66 = m43/127.;
float mx71 = m44/127.;
float mx72 = m45/127.;
float mx73 = m46/127.;
float mx74 = m47/127.;
float mx75 = m48/127.;
float mx76 = m49/126.;
float mx81 = m50/127.;
float mx82 = m51/127.;
float mx83 = m52/127.;
float mx84 = m53/127.;
float mx85 = m54/127.;
float mx86 = m55/127.;
float mx91 = m56/127.;
float mx92 = m57/127.;
float mx93 = m58/127.;
float mx94 = m59/127.;
///////////////// END OF MIDIMIX ////////////////////



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
        //u = u*u;
        u = u*u*u*((6.0*u-15.0)*u+10.0);

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
	// correct aspect ratio transfer
	xp = xp * asp;
	float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
	vec3 cl = dist < r ? col : vec3(0.);
	return cl;
}

// circle ring
vec3 crcr(in float xp, in float yp, in float r, in float w, in vec3 col)
{
        // correct aspect ratio transfer
        xp = xp * asp;
        float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
        vec3 cl = (((dist < r) && (dist > r - w/2.)) || ((dist > r) && (dist < r + w/2.)))  ? col : vec3(0.);
        return cl;
}


// circle-test
vec3 crcx(in float xp, in float yp, in float r, in vec3 col)
{
        float dist = sqrt( (p.x*asp*tan(time) - xp) * (p.x*asp - xp) + (p.y*tan(time) - yp) * (p.y - yp*sin(time)) );
        vec3 cl = dist < r ? col : vec3(0.);
        return cl;
}

vec3 oscx(in float xp, in float w, in float sinm, in float sinm2, in vec3 c)
{
	float cl = 0.5*w / abs(4. * (p.x - xp) + sin(p.y * sinm * 100.)*sinm2);
	return cl * c;
}

vec3 oscxt(in float xp, in float w, in float sinm, in float sinm2, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.x - xp) + sin(p.y * 10. * sinm) + tan(p.y * sinm * 150.)*sinm2);
        return cl * c;
}

vec3 oscy(in float yp, in float w, in float sinm, in float sinm2, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.y - yp) + cos(p.x * sinm * 100.)*sinm2);
        return cl * c;
}

//tst
float plot(vec2 p, float px){
  //return  smoothstep( pct-0.01, pct, p.y) - smoothstep( pct, pct+0.01, p.y);
	return smoothstep(px, px+1., p.y) + smoothstep(px, px-1., p.y);
}




void main( void ) {
////////////////////// INIT ////////////////////////
vec3 col;
// position
//p = ( gl_FragCoord.xy / resolution.xy );
p = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
ppp = p;
p = vec2(1.-p.x, p.y);
//p.x = p.x / 1.8;
//p.y = (p.y + 0.5) / 1.5;
vec2 mx = ( gl_FragCoord.xy / mouse.xy );
// aspect ratio
asp = resolution.x / resolution.y;
//asp = 1.0;
float px, py, spd;
vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cg = vec3(0.,1.,0.);
vec3 cb = vec3(0.,0.,1.);
vec3 ck = vec3(0.);

vec3 ca,cc,cd,ce,cf;

////////////////////// END OF INIT /////////////////
c = vec3(0.);


float RAYY = 1.0;
if (RAYY == 1.0) {
 vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-m5;
 vec4 d=p_ray;
 vec4 t;
    p_ray.z += time;
    for(float i=4.8; i>0.; i-=.1)
    {
        t = sin(abs(mod(p_ray, 0.5)-.25)*1.8);
        float x = max(t.y*2., length(t.xz)*2.);
//        c = vec3(i*.3,0.,0.)*0.1;
	c = vec3(mod(sin(i)*20.,p.x*1.),rand.x,rand.y)*4.*m1;
	ca = vec3(rand.x, mod(sin(i)*20.,p.x*1.),m4); // texture use
        if(x<.2) break;
        p_ray -= d*x;
     }

}

// BW gradient
float BW = 0.0;
if (BW == 1.0) {
	//vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-vec4(m10);
	// adjusting for 1440x900
	vec4 p_ray = (vec4(p,p)-m10)*10.; // set size here
	vec4 d=p_ray;
	vec4 t;
	p_ray.z += time;
	for(float i=4.8; i>0.; i-=.1)
	{
		  // this loops between values
		  //t = abs(mod(p_ray*cnt(1000)*5., 0.5)-.25);
		  // this runs
		  //t = abs(mod(p_ray*time, 0.5)-.25);
		  // this bounces
		  t = abs(mod(p_ray*sin(cnt(500)*3.14)*.1, 0.5)-.25);
		  //t = abs(mod(p_ray*4., 0.5) -.25);
		  // forms working
		  float x = max(t.x*2., length(t.yy)*2.);
		  // new forms - dots
		  //float x = max(t.y*2., length(t.xy)*2.);
		  c += vec3(0.,noise2f( vec2( tan(i),i ) )*m0*100.,0.);
		  cf = vec3(noise2f( vec2( tan(i),i ) )*m0*100.);
		  if(x<.15) break;
		  // normal
		  //p_ray -= d*x;
		  // new forms - dost
		  p_ray -= d*x*x*x*x;
	}
}

// RAY2
float RAY = .0;
if (RAY == 1.0)
{
	vec4 p_ray = (vec4(p,p)-m10)*1.;
	vec4 d=p_ray;
	vec4 t;
	p_ray.z += time;
	for(float i=4.8; i>0.; i-=.1)
	{
		  // this loops between values
		  //t = abs(mod(p_ray*cnt(1000)*5., 0.5)-.25);
		  // this runs
		  //t = abs(mod(p_ray*time, 0.5)-.25);
		  // this bounces
		  //t = abs(mod(p_ray*a_sin(time/2.)*5., 0.5)-.25);
		  t = abs(mod(p_ray*0.5+time/m1*100., 0.5) -.25);
		  // forms working
		  float x = max(t.y*2., length(t.xx)*2.);
		  // new forms - dots
		  //float x = max(t.y*2., length(t.xy)*2.);
		  c += vec3(noise2f(vec2(time/100.,i))*.6, 0.,0.)*0.2;
		  if(x<.15) break;
		  // normal
		  //p_ray -= d*x;
		  // new forms - dost
		  p_ray -= d*x*x;
	}
}

// tapeta1
float TAP = 0.0;
if (TAP == 1.0) {
	px = p.x + 0.5;
	py = p.y + 0.3;
	spd = 0.1;
	spd = m1/100.;
	vec2 va = vec2( px*time*spd, tan(py*2.) );
	ca = vec3( noise2f(va*va*va), noise2f(va*va), noise2f(va*va) );
	c-=liney(0.5,m1*1.3, ca*3.)*m2;
}

// tapeta2
float TAPP = 0.0;
float iter = 20.;
if (TAPP == 1.0)
{
	for(float i=0.; i<iter; i+=1.)
	{
		px = p.x + 0.5;
		py = p.y + 0.3;
		spd = 0.008*i;
		//spd = m1;
		vec2 va = vec2( sin(px)*time*spd*m1, floor(tan(py*time*spd)) );
		ca = vec3( noise2f(va*va), noise2f(va*va), noise2f(va*va) );
		c+=liney(i/iter,1./iter, ca*3.)*m0;
	}
}
cc = c;
//c += crc(0.5,0.5,0.3+cnt(80)*m3*0.3,cc)*0.8;

// OSCES I @ mx5*
//
float OSC = mx52;
if (OSC == 1.)
{
	float oscmult = mx56*20.;
	float step = 1./oscmult;
	float alpha = mx51;
	for (float i=step; i <= oscmult; i+=step)
	{
		c+=oscxt(i*step,step,m1,m2,cr)*alpha;
	}
}

// OSCES II @ mx6*
// vertical sinusoids traveling at *period
float OSCC = mx62;
if (OSCC == 1.)
{
	int period = (int)(mx66*10000.);
	float weight = mx65*10.;
	float alpha = mx61*100.;
	c+=oscx(cnt(period),weight,m1,m2,ca)*m1*alpha;
}

// BACKBUFFER1 @ mx7*
// not really a feedback, this is a treshold filter
// it needs something happening on c (like mx6* for example)
float BCKK = mx72;
if (BCKK == 1.)
{
	ppp = ppp.xy;
	float xs = mx76;
	float ys = mx75;
	float size = mx74;
	float treshold = mx71;
	vec3 ccc = texture2D(backbuffer,((ppp-xs)*(1.-ys)+0.5)*size).rgb;
	c*=floor(mix(ccc,c,0.1)*treshold*2.);
}

// BACKBUFFER2 @ mx8*
// this is like the classical feedback
// also its somehow flipped like (1-x,1-y) (2 lazy 2 deal w rn)
// TODO fix:
//		-- the weird flip on ppp
//		-- when changing size dont change position
float BCK = mx82;
if (BCK == 1.)
{
	ppp = ppp.xy;
	float xs = mx86-.5;
	float ys = (mx85-.5)*asp;
	float size = (1.-mx84)+.5;
	float alpha = mx81;
	float invert = mx83;
	vec3 ccc = texture2D(backbuffer, (ppp-vec2(xs,ys))*size).rgb;
	ccc = vec3(abs(invert-ccc.r), abs(invert-ccc.g), abs(invert-ccc.b));
	ccc *= alpha;
	c+=ccc;
}

// this is like the final strobe
//c*=cnt(50);

// DRAW
gl_FragColor = vec4(c, 1.0);

}
