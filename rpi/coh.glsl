#ifdef GL_ES

// TODO
// - test rand
// - test PD connection

// notes
// - fragcoord.xy / resolution - normalizes all coordinates between 0 and 1
//   it is fine when switching various resolutions, but not good for precision positioning

precision mediump float;
#endif
#define PI 3.141592

//////////////////// GLOBALS /////////////////

uniform float time;
uniform vec2 resolution;
uniform sampler2D backbuffer;
uniform vec2 rand;
uniform float s1;
uniform float s2;
uniform float s3;
uniform float s4;
uniform float s5;
uniform float s6;
uniform float s7;
uniform float s8;
uniform vec4 bng;
uniform vec4 knt1;
uniform vec4 knt2;
uniform vec4 glb;

// MEMORY
uniform vec4 scn0_a;
uniform vec4 scn0_b;
uniform vec4 scn1_a;
uniform vec4 scn1_b;
uniform vec4 scn2_a;
uniform vec4 scn2_b;
uniform vec4 scn3_a;
uniform vec4 scn3_b;
uniform vec4 scn4_a;
uniform vec4 scn4_b;

vec2 p; //relative position
vec3 c; //color
float mw = 0.001; //minimal width
float asp; //aspect ratio
float t = time;

//////////////////// END OF GLOBALS ////////////////////////

// counter 0 to 1
float cnt(in int m)
{
        float divider = 1000.0 / float(m)*10.0;
        return mod(time*divider, 10.0) / 10.0;
}

// sin between 0 - 1
float zsin(in float a) {
	return (sin(a) + 1.0) / 2.0;
}

// cos between 0 - 1
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

// old vertical line
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
        vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? -c+col : vec3(0.0);
        return cl;
}

//stripe
vec3 ostr(in float xp, in float yp, in float x, in float y, in vec3 col)
{
        //float ys = y * asp;
        float ys = y;
        vec3 cl = ((p.x > xp- x/2.) && (p.x < xp + x/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
        return cl;
}

//stripe 
vec3 str(in float xp, in float yp, in float x, in float y, in vec3 c, in vec3 col)
{
        float ys = y * asp;
        vec3 cl = ((p.x > xp- x/2.) && (p.x < xp + x/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? -c+col : vec3(0.0);
        return cl;
}


// old circle
vec3 ocrc(in float xp, in float yp, in float r, in vec3 col)
{
        float dist = sqrt( (p.x - xp) * (p.x - xp) + (p.y - yp) * (p.y - yp) );
        vec3 cl = dist < r ? col : vec3(0.);
        return cl;
}

// circle
vec3 crc(in float xp, in float yp, in float r, in vec3 c, in vec3 col)
{
	float dist = sqrt( (p.x - xp) * (p.x - xp) + (p.y - yp)/asp * (p.y - yp)/asp );
	vec3 cl = dist < r ? -c+col : vec3(0.);
	return cl;
} 

// circle-test
vec3 crcx(in float xp, in float yp, in float r, in vec3 col)
{
        float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y- yp) * (p.y - yp) );
        vec3 cl = dist < r ? col : col * (max(0.,r-dist+0.1)*10.);
        return cl;
}

// fake time
float faket()
{
	return cnt(3600*1000)*3600.;
}

vec3 oscx(in float xp, in float w, in float sinm, in float sinm2, in vec3 c)
{
	float cl = 0.5*w / abs(4. * (p.x - xp) + sin(p.y * sinm * 100.)*sinm2);
	return cl * c;
} 

vec3 oscy(in float yp, in float w, in float sinm, in float sinm2, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.y - yp) + cos(p.x * sinm * 100.)*sinm2);
        return cl * c;
}

vec3 oscq(in float clstr, in float sinm1, in float sinm2)
{
	vec2 uv = p - 0.5;
	vec3 col = vec3( 0.0 );
	float tt;
	tt = .01 / abs( uv.x + .5 * sin(uv.y*100.*p.x*sinm1));
	col += vec3( tt );
	tt = .001 / abs( uv.x + 0.5 * cos(uv.y*100.*p.y*sinm2));
	col += vec3( tt );
	tt = .001 / abs( uv.x + .5 * sin(uv.y*50.*p.x*sinm1));
        col += vec3( tt );
        tt = .001 / abs( uv.x + 0.5 * cos(uv.y*10.*p.y*sinm2));
        col += vec3( tt );
        return col;
}
	
vec3 feedb(in float xpos, in float ypos, in float xsiz, in float ysiz, in float bsiz, in float blnd, in float brg, in vec3 c) 
{
        vec2 mm = (p.x > xpos-xsiz/2.) && (p.x < xpos+xsiz/2.) && (p.y > ypos-ysiz/2.) && (p.y < ypos+ysiz/2.) ? (vec2(p.x,p.y)-0.5)*bsiz+0.5 : vec2(0.);
        vec3 ccc = texture2D(backbuffer, mm).xyz;
	if (blnd == 1.) {
        	// when too much light use mix
        	c=mix(str(xpos,ypos,xsiz,ysiz,c,ccc),c,brg);
	} else {
        	// when not much light write directly
        	c+=str(xpos,ypos,xsiz,ysiz,c,ccc);
	}
	return c;
}

float tacol(in float mod1) {
	float ta = 1.-noise2f(tan(p)*mod1*100.);
        ta = ta*ta*ta;
	return ta;
}

float powcol(in float mod1) {
	vec2 pp = p;
	float tt = 30.+sin(t)*10.*mod1;
	float cx = cos(t)/100.;
	pp+= noise2f(pp*tt);
	cx+= noise2f(pp*tt);
	float cxx=cx*cx*10.4;
	return 1.-cxx; 
}

vec3 tata(in float mod1, in float mod2, in float w, in vec3 c) {
	float ta = 1.-noise2f(sin(p)*mod1*100.);
	ta = ta*ta*ta*ta;
	c+= str(0.5,mod2,1.2, w, c, vec3(1.,1.,1.)) * ta;
	return c;
}

vec3 tbtb(in float mod1, in float mod2, in float w, in vec3 c) {
        float ta = 1.-noise2f(cos(p)*mod1*100.);
        ta = ta*ta*ta*ta;
        c+= str(0.5,mod2,1.2, w, c, vec3(1.,1.,1.)) * ta;
        return c;
}




void main( void ) {
////////////////////// INIT ////////////////////////
// position
p = ( gl_FragCoord.xy / resolution.xy );
// aspect ratio
asp = resolution.x / resolution.y;
const vec3 cw = vec3(1.);
const vec3 cb = vec3(0.);
const vec3 cr = vec3(1.,0.,0.);

vec3 ccc;

////////////////////// END OF INIT /////////////////
vec3 c;
c = cb;

////////////////////// CONTROL MAPPER ///////////////
// GLOBAL
//g1 - unused
//g2 - reserved
// backbuffer
float BCKX = glb.z;
float BCKY = glb.w;

// SCENES 
// scene0 - globals
float glb_r = scn0_a.x;
float glb_g = scn0_a.y;
float glb_b = scn0_a.z;
float glb_bw = scn0_a.w;
float glb_mix = scn0_b.x*5.-.1;
// scene1 - scenes
float scene2 = scn1_a.x;
float scene3 = scn1_a.y;
// scene2 - osces
float scn2_w1 = scn2_a.x; 
float scn2_w2 = scn2_a.y;
float scn2_x1 = scn2_a.z;
float scn2_y1 = scn2_a.w;
float scn2_x2 = scn2_b.x;
float scn2_y2 = scn2_b.y;
// scene3 - feedback
float scn3_xp = scn3_a.x;
float scn3_yp = scn3_a.y;
float scn3_xs = scn3_a.z*1.5;
float scn3_ys = scn3_a.w*1.5;
float scn3_bs = scn3_b.x*10.-5.;
float scn3_bl = scn3_b.y;
float scn3_size_control = scn3_b.z;
float scn3_bright = scn3_b.w*5.-0.1;


////////////////////// SCENES /////////////////////////


// BANGS
if (bng.x == 1.) {
	c=vec3(1.);
}
//if (bng.y == 1.) {
	 c+=crc(0.5,0.5,cnt(1000)/5.,c,cw)*tacol(rand.x*10.);
	//c+=crc(0.5,0.5,cnt(1000)*5.,c,cw)*powcol(s2);
//}
if (bng.z == 1.) {
	c+=tata(rand.x*10.,rand.x,rand.y,c);
}


// SCENE 3 - FEEDBACK
if (scene3 == 1.) {
        //float xp = scn3_xp;
        //float yp = scn3_yp;
	float xp = .5;
	float yp = .5;
        float xs = scn3_xs;
        float ys = scn3_ys;
	float bs = 0.5;
	float brg = 1.;
	if (scn3_size_control == 1.) {
        	bs = scn3_bs;
	} else {
		bs = s3;
	}	
        float bl = scn3_bl;
	if (bl == 1.) {
		brg = scn3_bright;
	}
        c = feedb(xp, yp, xs, ys, bs, bl, brg, c);
}

if (OSC == 1.) {
        float cent = 0.5;
        float w1 = mx46;
	float w2 = mx45;
		
	if (mx62 > 1.) {
        	c+=oscx(cent, w1, m1, m2, vec3(2.,1.,.5));
	}
	if (mx62 > 3.) {
        	c+=oscy(cent, w1, m1, m2, vec3(2.,1.,.5));
	}
	if (mx62 > 3.)) {
        	c+=oscx(m3, w2, m3, m4, vec3(.5,1.,2.));
	}
	if (mx62 > 4.) {
        	c+=oscy(cent, w2, m3, m4, vec3(.5,1.,2.));
	}
}

// BACKBUFFER
// BS1 or BM1
if ((BCKX == 1.) || (BCKY == 1.)) {
	if (BCKY == 1.) {
		p = p.yx;
	}
	ccc = texture2D(backbuffer, (p-s2)*(1.-s1)+0.5).xyz;
	c=mix(ccc,c,0.1);
}

// COLOR STUFF
ccc = c;
if (glb_bw == 1.) {
	ccc = vec3(dot(c.rgb, vec3(0.299, 0.587, 0.114)));
	ccc = mix(ccc,vec3(glb_r,glb_g,glb_b),glb_mix);
} else {
	ccc = vec3(ccc.r*glb_r,ccc.g*glb_g,ccc.b*glb_b);
}

gl_FragColor = vec4( ccc, 1.0 );

}

