//precision mediump float;

//////////////////// GLOBALS /////////////////
uniform float time;
uniform vec2 resolution;
uniform sampler2D backbuffer;
uniform vec2 rand;
vec2 p; //relative position
vec3 c; //color
float mw = 0.001; //minimal width
float asp; //aspect ratio

/// RAYMARCH
#define PI 3.1415926535898
const float eps = 0.005;
const int maxIterations = 256;
const float stepScale = 0.5;
const float stopThreshold = 0.005;
const float clipNear = 0.0;
const float clipFar = 4.0;

/// PARAMS
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;
uniform float m5;
uniform float m6;
uniform float m7;

/// COLORS
vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cg = vec3(0.,1.,0.);
vec3 cb = vec3(0.,0.,1.);
vec3 ck = vec3(0.);
vec3 ca,cc,cd,ce,cf,cx;
//////////////////// END OF GLOBALS //////////////

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
float mx76 = m49/127.;
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

/////////////////////////////////////// SIMPLE FUNCTIONS

// Sin between 0 - 1
float zsin(in float a) {
	return (sin(a) + 1.0) / 2.0;
}

// Sin between 0 - 1 - bouncing
float a_sin(in float a) {
	return abs(sin(a));
}

// Cos between 0 - 1
float zcos(in float a) {
        return (cos(a) + 1.0) / 2.0;
}

// Cos between 0 -1 - bouncing
float a_cos(in float a) {
        return abs(cos(a));
}

// Modulo
float mdl(in float a, in float b)
{
        return a - b * floor(a/b + 0.001);
}

// Random number generator
float rnd(vec2 co){
        // implementation found at: lumina.sourceforge.net/Tutorials/Noise.html
        float k = sin(dot(co.xy ,vec2(12.9898,78.233)));
        return fract(k + k);
}

// Noise2f
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

// Counter 0 to 1 
float cnt(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

////////////////////////////////////// SIMPLE PRIMITIVES

// Vertical line
vec3 linex(in float xp, in float width, in vec3 col)
{
        return ((p.x <= xp+width/2.) && (p.x >= xp-width/2.)) ? col : vec3(0.0);
}

// Horizontal line
vec3 liney(in float yp, in float width, in vec3 col)
{
        return ((p.y <= yp+width/2.0) && (p.y >= yp-width/2.0)) ? col : vec3(0.0);
	
}

// Square
vec3 sqr(in float xp, in float yp, in float s, in vec3 col)
{
	float ys = s * asp;
	vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-s/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
	return cl;
}

// Rectangle
vec3 rect(in float xp, in float yp, in float sx, in float sy, in vec3 col)
{
	float syy = sy*asp;
	vec3 cl = ((p.x > xp - sx/2.) && (p.x < xp + sx/2.) && (p.y > yp-syy/2.) && (p.y < yp + syy/2.)) ? col : vec3(0.0); 
	return cl;
}

// Stripe 
vec3 str(in float xp, in float yp, in float x, in float y, in vec3 c, in vec3 col)
{
        float ys = y * asp;
        vec3 cl = ((p.x > xp- x/2.) && (p.x < xp + x/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? -c+col : vec3(0.0);
        return cl;
}

// Circle
vec3 circ(in float xp, in float yp, in float r, in vec3 col)
{
	// correct aspect ratio transfer
	xp = xp * asp;
	float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
	vec3 cl = dist < r ? col : vec3(0.);
	return cl;
} 

// Circle - remove background
vec3 circ_del(in float xp, in float yp, in float r, in vec3 c, in vec3 col)
{
        // correct aspect ratio transfer
        xp = xp * asp;
        float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
        vec3 cl = dist < r ? -c+col : vec3(0.);
        return cl;
}

// Circle ring
vec3 circr(in float xp, in float yp, in float r, in float w, in vec3 col)
{
        // correct aspect ratio transfer
        xp = xp * asp;
        float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
        vec3 cl = (((dist < r) && (dist > r - w/2.)) || ((dist > r) && (dist < r + w/2.)))  ? col : vec3(0.);
        return cl;
}

// Oscillator on the x-axis
vec3 oscx(in float xp, in float w, in float sinm, in float sinm2, in vec3 c)
{
	float cl = 0.5*w / abs(4. * (p.x - xp) + sin(p.y * sinm * 100.)*sinm2);
	return cl * c;
} 

// Oscillator ??
vec3 oscxt(in float xp, in float w, in float sinm, in float sinm2, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.x - xp) + sin(p.y * 10. * sinm) + tan(p.y * sinm * 150.)*sinm2);
        return cl * c;
}

// Oscillator on the y-axis
vec3 oscy(in float yp, in float w, in float sinm, in float sinm2, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.y - yp) + cos(p.x * sinm * 100.)*sinm2);
        return cl * c;
}


/////////////////////////////////////////// TEXTURES AND PAINTING
vec3 feedb(in float xpos, in float ypos, in float xsiz, in float ysiz, in float bsiz, in float blnd, in float brg, in vec3 c) 
{
	if (mx62 == 1.) {
		float dist = sqrt( (p.x - xpos) * (p.x - xpos) + (p.y - ypos) * (p.y - ypos) );
        	vec2 mm = dist < xpos*asp ? (vec2(p.x,p.y)-0.5)*bsiz+0.5 : vec2(0.);
        	vec3 ccc = texture2D(backbuffer, mm).xyz;
		c+=circ_del(xpos,ypos,xsiz,c,ccc);
	} else {
	        vec2 mm = (p.x > xpos-xsiz/2.) && (p.x < xpos+xsiz/2.) && (p.y > ypos-ysiz/2.) && (p.y < ypos+ysiz/2.) ? (vec2(p.x,p.y)-0.5)*bsiz+0.5 : vec2(0.);
        	vec3 ccc = texture2D(backbuffer, mm).xyz;
        	c+=str(xpos,ypos,xsiz,ysiz,c,ccc);
	}
	return c;
}


/////////////////////////////////////////// GEOMETRIC MODIFIERS
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
 //   vec3 q = inverse(m)*p;
  //  return q;
//}


///////////////////////////////////////////// RAYMARCHING STUFF

// Sphere field
float sphere(in vec3 p, in vec3 centerPos, float radius) {
        return length(p-centerPos) - radius*0.8;
}

// Box field
float box( vec3 p, vec3 b ) {
  return length(max(abs(p)-b*m1,0.0));
}

float hex( vec3 p, vec2 h )
{
    vec3 q = abs(p);
    return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
}

// Sinusoid bumps
float sinusoidBumps(in vec3 p){
        float frq = 10.;
        return 3.*sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31);
}

// The whole scene
float scene(in vec3 p) {
	if (mx12 == 0.) {
	return box(p, vec3(0.3,0.3,0.3)) + mx16*sinusoidBumps(p)*m0;	
	}
	if (mx12 == 1.) {
	return sphere(p, vec3(0., 0. , 2.), 1.) + mx16*sinusoidBumps(p);
        }
	if (mx12 == 2.) {
	return hex(p*mx12, vec2(0.1,0.2)) + mx16*sinusoidBumps(p);
        }
}

// Gets the surface normal for p
vec3 getNormal(in vec3 p) {
        return normalize(vec3(
                scene(vec3(p.x+eps,p.y,p.z))-scene(vec3(p.x-eps,p.y,p.z)),
                scene(vec3(p.x,p.y+eps,p.z))-scene(vec3(p.x,p.y-eps,p.z)),
                scene(vec3(p.x,p.y,p.z+eps))-scene(vec3(p.x,p.y,p.z-eps))
        ));
}

// Raymarches
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

// lighting
vec3 lights() {
        vec2 aspect = vec2(asp, mx14*70.);
        vec2 screenCoords = (2.0*gl_FragCoord.xy/resolution.xy - 1.0)*aspect*.5;
        vec3 lookAt = vec3(0.,0.,0.);
        //rotate by hand
        mx11*=10.;
        vec3 camPos = vec3(sin(mx11), sin(mx11), cos(mx11));
        //static
        //vec3 camPos = vec3(0., 0., -1.);
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
            //c += vec3(bgcolor);
            //gl_FragColor = vec4(bgcolor, .3);
            //BACK = vec3(bgcolor);
               // c=ck;
            return;
        }
        vec3 sp = ro + rd*dist*mx15*m2*10.;
        vec3 surfNormal = getNormal(sp*sp);
        //vec3 lp = vec3(1.5*sin(time*.5), 0.75+0.25*cos(time*0.5), -1.0);
        vec3 lp = vec3(noise2f(p),noise2f(p*p),noise2f(p*p*p)*m0);
        vec3 ld = lp-sp;
        vec3 lcolor = vec3(1.,1.,1.);
        float len = length( ld );
        ld /= len;
        float lightAtten = min( 1.0 / ( 0.25*len*len ), 1.0 ); // Keeps things between 0 and 1.
        vec3 ref = reflect(-ld, surfNormal);
        vec3 sceneColor = vec3(0.0);
        vec3 objColor = vec3(1., .0, 0.);
        float bumps =  sinusoidBumps(sp);
        objColor += clamp(objColor-vec3(0.3, 0.4, 0.)*bumps, 0.0, 1.0);
        float ambient = 1.3;
        float specularPower = 100.0;
        float diffuse = max( 0.0, dot(surfNormal, ld) );
        diffuse = pow(diffuse, 1000.);
        float specular = max( 0.0, dot( ref, normalize(camPos-sp)) );
        specular = pow(specular, specularPower);
        sceneColor = vec3(0.);
        sceneColor += (objColor*(diffuse*0.8+ambient)+specular*0.5)*lcolor*lightAtten;
	return sceneColor;
} 



void main(void) {
	/// INIT 
	p = vec2( gl_FragCoord.x / resolution.x, gl_FragCoord.y / resolution.y);
	asp = resolution.x / resolution.y;

	c= ck;

float BLOB = mx13;
if (BLOB == 1.0) {
	c += lights();
}

float RAY = 1.0;
if (RAY == 1.) {
 vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-0.4;
	p_ray = vec4(p, vec2(mx55,0.))-0.5;
 vec4 d=p_ray;
 vec4 t;
    p_ray.z += time;
    for(float i=mx56*5.; i>0.; i-=.1)
    {
        t = abs(mod(p_ray, 0.5)-.25*m0);
        float x = max(t.y*2., length(t.xz)*2.*m1);
	if (mx52 == 0.) {
		c = vec3(i)*mx51;
	}
	if (mx52 == 1.) {
        	c = vec3(sin(mod(i*p.y*i*10.,2.)), tan(i*i*p.y), tan(i/p.x))*mx51;
	}
	if (mx52 == 2.) {
		c = vec3(sin(100.*abs(noise2f(i*p.x*10./p.y))))*mx51;
        }
        if(x<.2) break;
        p_ray -= d*x;
     }
}

float OSC=1.;
if (OSC == 1.) {
        float cent = 0.5;
        float w1 = mx46;
	float w2 = mx45;
		
	if (mx42 >= 1.) {
        	c+=oscx(cent, w1, m1, m2, vec3(2.,1.,.5))*mx41;
	}
	if (mx42 >= 2.) {
        	c+=oscy(cent, w1, m1, m2, vec3(2.,1.,.5))*mx41;
	}
	if (mx42 >= 3.) {
        	c+=oscx(m3, w2, m3, m4, vec3(.5,1.,2.))*mx41;
	}
	if (mx42 >= 4.) {
        	c+=oscy(cent, w2, m3, m4, vec3(.5,1.,2.))*mx41;
	}
}

// FEEDBACK
if (mx63 == 1.) {
	float xp = .5;
	float yp = .5;
        float xs = mx66*m2;
        float ys = mx65*m1;
	float bs = mx64*10.*m0;
	for (int i=0; i<10; i++) {
       		cx= feedb(i/10.+0.05, yp, xs, ys, bs, 0., 1., c);
		 c = mix(c,cx,mx61);
	}
//	c = mix(c,cx,mx61);
}


// BACKBUFFER
mx75 = mx75*2.-1.; 
if (mx73 == 1.) {
	// normal / inverse y / inverse x
	if (mx72 == 1.) {
		p = vec2(p.x,1.-p.y);
	}
	if (mx72 == 2.) {
                p = vec2(1.-p.x,p.y);
        }
	if (mx72 == 3.) {
		p = vec2(p.y,p.x);
	}
	cx = texture2D(backbuffer, (p-mx76)*(1.-mx75)+0.5).xyz;
	//threshold
	c=mix(c,cx,mx71);
}



// COLOR STUFF
cx = c;
if (mx83 == 1.) { 
	cx = vec3(dot(c.rgb, vec3(0.299, 0.587, 0.114)));
	cx = mix(cx,vec3(mx86,mx85,mx84),mx81);
} else {
	cx = vec3(cx.r*mx86,cx.g*mx85,cx.b*mx84);
}

cx+=vec3(cnt(50))*mx32;


// DRAW
gl_FragColor = vec4(cx, 1.0);

}
