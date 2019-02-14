#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
uniform vec2 resolution;
uniform vec4 mouse;
uniform vec2 rand;
uniform sampler2D backbuffer;

///
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;
uniform float m5;
uniform float m6;
uniform float m7;
uniform float m8;
uniform float m9;
uniform float m10;
uniform float m11;


#define PI 3.1415926535898 // That was from memory, so if things start flying off the screen...
const float eps = 0.005; 
const int maxIterations = 256;
const float stepScale = 0.5;
const float stopThreshold = 0.005; // I'm not quite sure why, but thresholds in the order of a pixel seem to work better for me... most times.

float sphere(in vec3 p, in vec3 centerPos, float radius) {
	return length(p-centerPos) - radius*0.8;
}

float box( vec3 p, vec3 b ) {
  return length(max(abs(p)-b,0.0));
}

float sinusoidBumps(in vec3 p){
	float frq = 10.;
    	//return 15.*sin(p.x*frq+time*0.57)*cos(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31) + 0.5*sin(p.x*32.+time*0.07)*cos(p.y*32.+time*2.11)*sin(p.z*32.-time*1.23);
	return 3.*sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31);
}

float scene(in vec3 p) {
	//mult
	p.xy *= 0.5;
	//p.x = mod(3.*p.x,1.);
	//p.y = mod(3.*p.y,1.);
	//rot
	//vec3 q = invert(m)*p;
	//p *= 0.4;
	return sphere(p, vec3(0., 0. , 2.), 1.) + 0.04*sinusoidBumps(p);
	//return smoothstep(box(p-vec3(0.,0.,1.3+tan(time)-1.), vec3(1.5,0.2,0.1)),sphere(p, vec3(0., 0. , 2.), 1.) + 0.04*sinusoidBumps(p), 100.);
	p = p - vec3(0.,-1.6,1.6);
	float c = cos(20.*p.x);
    float s = sin(20.*p.x);
    mat2  m = mat2(c,-s,s,c);
    vec3  q = vec3(m*p.xz,p.y);
	p = q;
	//return min(box(p-vec3(0.,0.,1.3+sin(time/2.)), vec3(1.5,0.2,0.1)),sphere(p, vec3(0., 0. , 2.), 1.) + 0.04*sinusoidBumps(p));
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

void main(void) {
	vec3 c = vec3(0.);

	vec2 aspect = vec2(resolution.x/resolution.y, 1.0);
	vec2 screenCoords = (2.0*gl_FragCoord.xy/resolution.xy - 1.0)*aspect;
	vec2 p = screenCoords;
	vec3 lookAt = vec3(0.,0.,0.);
	vec3 camPos = vec3(0., 0., -1.); 
    vec3 forward = normalize(lookAt-camPos);
    vec3 right = normalize(vec3(forward.z, 0., -forward.x ));
    vec3 up = normalize(cross(forward,right));
    float FOV = .9;

    vec3 ro = camPos; 
    vec3 rd = normalize(forward + FOV*screenCoords.x*right + FOV*screenCoords.y*up);

    vec3 bgcolor = vec3(1.,0.97,0.92)*0.15;

	const float clipNear = 0.0;
	const float clipFar = 4.0;
	float dist = rayMarching(ro, rd, clipNear, clipFar );
	vec3 BACK = vec3(1.1); 
	if ( dist >= clipFar ) {
	    //c += vec3(bgcolor);
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
	c = sceneColor;

	// BACKBUFFER
float BCK = 0.;
if (BCK == 1.) {
	vec2 ppp = p.xy;
	float xs = m8;
	float ys = m9;
float s1 = m1;
float s2 = m2;
	//ys = 0.5;
	//vec3 ccc = texture2D(backbuffer, vec2(ppp.x,1.0-ppp.y)*1.01).rgb;
	vec3 ccc = texture2D(backbuffer,((ppp-s2)*(1.-s1)+0.5)*m11).rgb;
	//vec3 ccc = vec3(texture2D(backbuffer,(ppp-s2)*(1.-s1)+0.5).r,0.,0.) * 10.;
	c+=mix(ccc,c,0.1);
}

	        	
	gl_FragColor = vec4(clamp(c, 0.0, 1.0), 1.9);
}
