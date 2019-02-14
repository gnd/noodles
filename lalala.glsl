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
uniform float m8, m9, m10, m11;

/// COLORS
vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cg = vec3(0.,1.,0.);
vec3 cb = vec3(0.,0.,1.);
vec3 ck = vec3(0.);
vec3 ca, cc, cd, ce, cf;
//////////////////// END OF GLOBALS //////////////

// //#include include.glsl

void main(void) {
	/// INIT 
	p = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
	asp = resolution.x / resolution.y;
	
	c= cr;

	// DRAW
	gl_FragColor = vec4(c, 1.0);
}
