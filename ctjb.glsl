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

// RAYMARCH
#define PI 3.1415926535898
const float eps = 0.005;
const int maxIterations = 256;
const float stepScale = 0.5;
const float stopThreshold = 0.005;
const float clipNear = 0.0;
const float clipFar = 4.0;

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

// others
uniform float m8, m9, m10, m11;


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


////////////////////////////////////////////////////////////////////////// GEOM
mat4 rmat(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s, oc * axis.z * axis.x + axis.y * s, 0.0,
                oc * axis.x * axis.y + axis.z * s, oc * axis.y * axis.y + c, oc * axis.y * axis.z - axis.x * s, 0.0,
                oc * axis.z * axis.x - axis.y * s, oc * axis.y * axis.z + axis.x * s, oc * axis.z * axis.z + c, 0.0,
                0.0, 0.0, 0.0, 1.0);
}

//vec3 rot( vec3 p, mat4 m )
//{
    //vec3 q = invert(m)*p;
    //return primitive(q);
//}


///////////////////////////////////////////////////////////////////// raymarch
float sphere(in vec3 p, in vec3 centerPos, float radius) {
        return length(p-centerPos)*m1 - radius*m0*2.;
}

float box( vec3 p, vec3 b ) {
  return length(max(abs(p)-b,0.0));
}

float sinusoidBumps(in vec3 p){
        float frq = 10.;
        return 3.*sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31);
}

float scene(in vec3 p) {
	p.xy *= m1;
	p.x = mod(p.x,m0)-m0;
	p.y = mod(p.y,m0)-m0;
        return sphere(p, vec3(0., 0. , 2.), 1.) + 0.04*sinusoidBumps(p)*m2*10.;
}

vec3 getNormal(in vec3 p) {
        return normalize(vec3(
                scene(vec3(p.x+eps,p.y,p.z))-scene(vec3(p.x-eps,p.y,p.z)),
                scene(vec3(p.x,p.y+eps,p.z))-scene(vec3(p.x,p.y-eps,p.z)),
                scene(vec3(p.x,p.y,p.z+eps))-scene(vec3(p.x,p.y,p.z-eps))
        ));
}

float rayMarching( vec3 origin, vec3 dir, float start, float end ) {
        float sceneDist = 1e4;
        float rayDepth = start;
        for ( int i = 0; i < maxIterations; i++ ) {
                sceneDist = scene( origin + dir * rayDepth );
                if (( sceneDist < stopThreshold ) || (rayDepth >= end)) {
                        break;
                }
                rayDepth += sceneDist * stepScale;
        }
        if ( sceneDist >= stopThreshold ) rayDepth = end;
        else rayDepth += sceneDist;
        return rayDepth;
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
c = ck;

float BLOB = 1.0;
if (BLOB == 1.0) {
	vec2 aspect = vec2(resolution.x/resolution.y, 1.0);
	vec2 screenCoords = (2.0*gl_FragCoord.xy/resolution.xy - 1.0)*aspect*.5;
	vec3 lookAt = vec3(0.,0.,0.);
	vec3 camPos = vec3(0., 0., -1.); 
	vec3 forward = normalize(lookAt-camPos);
	vec3 right = normalize(vec3(forward.z, 0., -forward.x ));
	vec3 up = normalize(cross(forward,right));
	float FOV = .9;

	vec3 ro = camPos; 
	vec3 rd = normalize(forward + FOV*screenCoords.x*right + FOV*screenCoords.y*up);
	vec3 bgcolor = vec3(1.,0.97,0.92)*0.15;

	float dist = rayMarching(ro, rd, clipNear, clipFar );
	vec3 BACK = vec3(1.1); 
	if ( dist >= clipFar ) {
	    c += vec3(bgcolor);
	    //gl_FragColor = vec4(bgcolor, .3);
	    //BACK = vec3(bgcolor);
	    return;
	}

	vec3 sp = ro + rd*dist;
	vec3 surfNormal = getNormal(sp);
	vec3 lp = vec3(1.5*sin(time*0.5), 0.75+0.25*cos(time*0.5), -1.0);
	vec3 ld = lp-sp;
	vec3 lcolor = vec3(1.,1.,1.);
	float len = length( ld );
	ld /= len; 
	float lightAtten = min( 1.0 / ( 0.25*len*len ), 1.0 ); // Keeps things between 0 and 1.

	vec3 ref = reflect(-ld, surfNormal); 
	vec3 sceneColor = vec3(0.0);
	vec3 objColor = vec3(1., .0, 0.);
	float bumps =  sinusoidBumps(sp*1.);
    	objColor += clamp(objColor-vec3(0.3, 0.4, 0.)*bumps, 0.0, 1.0);

	float ambient = 1.3; 
	float specularPower = 100.0;
	float diffuse = max( 0.0, dot(surfNormal, ld) ); 
	diffuse = pow(diffuse, 1000.);
	float specular = max( 0.0, dot( ref, normalize(camPos-sp)) ); 
	specular = pow(specular, specularPower);
	sceneColor = vec3(0.);
	sceneColor += (objColor*(diffuse*0.8+ambient)+specular*0.5)*lcolor*lightAtten;
	c += sceneColor*m4;
}

float RAYY = 1.0;
if (RAYY == 1.0) {
	vec4 p_ray = mod(vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-0.5,1.)-m1/1.;
	vec4 d=p_ray;
	vec4 t;
	p_ray.z += time*m0;
	for(float i=4.8; i>0.; i-=.1)
	{
        	t = mod(p_ray*time/50., 0.5)-.25;
        	float x = max(t.y*2., length(t.xx)*2.);
		c = vec3(i/3.,0.,0.)*m3;
		//c = vec3(mod(sin(i)*20.,p.x*1.),i,i)*4.;
		//ca = vec3(rand.x, mod(i*20.,p.x*1.),m4); // texture use
	       	if(x<.01) break;
        	p_ray -= d*x;
     	}
}

// BW gradient
float BW = 1.0;
if (BW == 1.0) {
	vec4 p_ray = (vec4(p,p)); // set size here
	vec4 d=p_ray;
	vec4 t;
	p_ray.z += time/1000.;
	for(float i=4.8; i>0.; i-=.1)
	{
		// this loops between values
		t = abs(mod(p_ray*cnt(500)*10., 0.5)-.25*m2);
		// this runs
		//t = abs(mod(p_ray*time, 0.5)-.25);
		//t = abs(mod(p_ray., 0.5) -.25);
		// forms working 
		float x = max(t.x*2., length(t.yy)*2.);
		// new forms - dots
		//float x = max(t.y*2., length(t.xy)*2.);
		c /= vec3(noise2f( vec2( tan(i),i ) ))/m1;
		cf = vec3(noise2f( vec2( tan(i),i ) ));
		if(x<.15) break;
		// normal
		//p_ray -= d*x;
		// new forms - dost
		p_ray -= d*x*x*x*x;
	}
}

// RAY2
float RAY = 1.0;
if (RAY == 1.0) {
	vec4 p_ray = (vec4(p,p))*2.;
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
		t = abs(mod(p_ray*0.5+time/10., 0.5) -.25);
		// forms working
		float x = max(t.y*2., length(t.xx)*2.);
		// new forms - dots
		//float x = max(t.y*2., length(t.xy)*2.);
		c += vec3(noise2f(vec2(time/100.,i))*.6, 0.,0.)*m2;
		if(x<.15) break;
		// normal
		//p_ray -= d*x;
		// new forms - dost
		p_ray -= d*x*x;
	}
}

// tapeta1
float TAP = 1.0;
if (TAP == 1.0) {
	px = p.x-m4; 
	py = p.y+0.3; 
	spd = 0.7;
//	spd = 1.;
	vec2 va = vec2( px*cnt(500)*6.+4., tan(py*2.) );
	ca = vec3( noise2f(va*va*va), noise2f(va*va), noise2f(va*va) );
	c*=liney(.6,1.3, ca*3.)*0.1;
}

// tapeta2
float TAPP = 1.0;
float iter = 20.;
if (TAPP == 1.0) {
	for(float i=0.; i<iter; i+=1.) {
		px = p.x + m0;
		py = p.y + m2;
		spd = 0.008*i;
		//spd = m1;
		vec2 va = vec2( sin(px)*time*spd, floor(tan(py*time*spd)) );
		ca = vec3( noise2f(va*va), noise2f(va*va), noise2f(va*va) );
		c/=liney(i/iter,1./iter, ca*3.)*m4;
	}
}

cc = c;
// pulzar
c += crc(0.5,0.5,m2/10.+cnt(80)*0.3,ca)*m3/10.;

// osces
float OSC = 1.0;
if (OSC == 1.) {
	float oscmult = 10.;
	for (float i=0.1; i <= oscmult; i+=1.) {
		c=oscx(i*0.1,0.01,rand.x,0.1,cr);
	}
}
c+=vec3(1.0);

// blinkroll osc
c+=oscx(cnt(100),0.1,0.1,0.2,cw)*m3;

// crcs
// concentric circs
float CRCS = 1.0;
if (CRCS == 1.0) {
	for (float i = 0.; i < .9; i+=.1) {
		c=crcr(.5,.5,.1+i*m1,0.1,cw);
	}
}

// BACKBUFFER
float BCK = .0;
if (BCK == 1.) {
	ppp = ppp.yx;
	float xs = m8;
	float ys = m9;
	float s1 = m8;
	float s2 = m9;
	//ys = 0.5;
	//vec3 ccc = texture2D(backbuffer, vec2(ppp.x,1.0-ppp.y)*1.01).rgb;
	vec3 ccc = texture2D(backbuffer,((ppp-s2)*(1.-s1)+0.5)*m11).rgb;
	//vec3 ccc = vec3(texture2D(backbuffer,(ppp-s2)*(1.-s1)+0.5).r,0.,0.) * 10.;
	c+=mix(ccc,c,0.9);
}


// DRAW
gl_FragColor = vec4(c, 1.0);



}

