#version 330
uniform vec2 resolution;
uniform float time;
uniform vec2 rand;
uniform sampler2D backbuffer;
out vec4 PixelColor;

/// PARAMS
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;

/// COLORS
vec3 cx;
vec3 ck = vec3(0.);
vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cg = vec3(0.,1.,0.);
vec3 cb = vec3(0.,0.,1.);
//////////////////// END OF GLOBALS //////////////

/////////////////// MIDIMIX //////////////////////
uniform float m5,m6,m7,m8,m9,m10;
uniform float m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30;
uniform float m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50;
uniform float m51,m52,m53,m54,m55,m56,m57,m58,m59;
float mx11 = m5/127.;
float mx12 = m6/127.;
float mx13 = m7/127.;
float mx14 = m8/127.;
float mx15 = m9/127.;
float mx16 = m10/127.;
float mx21 = m11/12.7;
float mx22 = m12/127.;
float mx23 = m13/127.;
float mx24 = m14/127.;
float mx25 = m15/127.;
float mx26 = m16/127.;
float mx31 = m17/127.;
float mx32 = m18/127.;
float mx33 = m19/127.;
float mx34 = m20/127.;
float mx35 = m21/127.;
float mx36 = m22/127.;
float mx41 = m23/127.;
float mx42 = m24/127.;
float mx43 = m25/127.;
float mx44 = m26/127.;
float mx45 = m27/127.;
float mx46 = m28/127.;
float mx51 = m29/127.;
float mx52 = m30/127.;
float mx53 = m31/127.;
float mx54 = m32/127.;
float mx55 = m33/127.;
float mx56 = m34/127.;
float mx61 = m35/127.;
float mx62 = m36/127.;
float mx63 = m37/127.;
float mx64 = m38/127.;
float mx65 = m39/127.;
float mx66 = m40/127.;
float mx71 = m41/127.;
float mx72 = m42/127.;
float mx73 = m43/127.;
float mx74 = m44/127.;
float mx75 = m45/127.;
float mx76 = m46/127.;
float mx81 = m47/127.;
float mx82 = m48/127.;
float mx83 = m49/127.;
float mx84 = m50/127.;
float mx85 = m51/127.;
float mx86 = m52/127.;
float mx91 = m53/127.;
float mx92 = m54/127.;
float mx93 = m55/127.;
float mx94 = m56/127.;
///////////////// END OF MIDIMIX ////////////////////
#define FLY

struct march { float step; vec3 ro; vec3 rd; };
struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material mat1;
shading s1;

vec2 p2;
vec3 p3;
float asp, spec;
vec3 sky, color; // 'sky color'
bool refracted = false;

///////////////////////////// NUMERICAL FUNCTIONS ///////////////////////////////////////////////////
mat4 inverse(mat4 m) {
  float
      a00 = m[0][0], a01 = m[0][1], a02 = m[0][2], a03 = m[0][3],
      a10 = m[1][0], a11 = m[1][1], a12 = m[1][2], a13 = m[1][3],
      a20 = m[2][0], a21 = m[2][1], a22 = m[2][2], a23 = m[2][3],
      a30 = m[3][0], a31 = m[3][1], a32 = m[3][2], a33 = m[3][3],

      b00 = a00 * a11 - a01 * a10,
      b01 = a00 * a12 - a02 * a10,
      b02 = a00 * a13 - a03 * a10,
      b03 = a01 * a12 - a02 * a11,
      b04 = a01 * a13 - a03 * a11,
      b05 = a02 * a13 - a03 * a12,
      b06 = a20 * a31 - a21 * a30,
      b07 = a20 * a32 - a22 * a30,
      b08 = a20 * a33 - a23 * a30,
      b09 = a21 * a32 - a22 * a31,
      b10 = a21 * a33 - a23 * a31,
      b11 = a22 * a33 - a23 * a32,

      det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

  return mat4(
      a11 * b11 - a12 * b10 + a13 * b09,
      a02 * b10 - a01 * b11 - a03 * b09,
      a31 * b05 - a32 * b04 + a33 * b03,
      a22 * b04 - a21 * b05 - a23 * b03,
      a12 * b08 - a10 * b11 - a13 * b07,
      a00 * b11 - a02 * b08 + a03 * b07,
      a32 * b02 - a30 * b05 - a33 * b01,
      a20 * b05 - a22 * b02 + a23 * b01,
      a10 * b10 - a11 * b08 + a13 * b06,
      a01 * b08 - a00 * b10 - a03 * b06,
      a30 * b04 - a31 * b02 + a33 * b00,
      a21 * b02 - a20 * b04 - a23 * b00,
      a11 * b07 - a10 * b09 - a12 * b06,
      a00 * b09 - a01 * b07 + a02 * b06,
      a31 * b01 - a30 * b03 - a32 * b00,
      a20 * b03 - a21 * b01 + a22 * b00) / det;
}

mat4 rmat(vec3 axis, float angle) {
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat4(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s, oc * axis.z * axis.x + axis.y * s, 0.0,
                oc * axis.x * axis.y + axis.z * s, oc * axis.y * axis.y + c, oc * axis.y * axis.z - axis.x * s, 0.0,
                oc * axis.z * axis.x - axis.y * s, oc * axis.y * axis.z + axis.x * s, oc * axis.z * axis.z + c, 0.0,
                0.0, 0.0, 0.0, 1.0);
}

// Random number generator
// lumina.sourceforge.net/Tutorials/Noise.html
float rnd(vec2 co){
    float k = sin(dot(co.xy ,vec2(12.9898,78.233)));
    return fract(k + k);
}

// http://www.iquilezles.org/www/articles/morenoise/morenoise.htm
float noise2f( in vec2 p ) {
    vec2 ip = vec2(floor(p));
    vec2 u = fract(p);
    u = u*u*u*((6.0*u-15.0)*u+10.0);
    float res = mix( mix(rnd(ip), rnd(ip+vec2(1.0,0.0)), u.x),
                     mix(rnd(ip+vec2(0.0,1.0)), rnd(ip+vec2(1.0,1.0)), u.x),
                    u.y);
    return res - 0.25;
}

// Counter 0 to 1
float cnt(int m) {
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

float smin( float a, float b, float k ) {
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

vec3 blend(in vec3 d1, in vec3 d2, in float k) {
	return k * d1 + (1.-k) * d2;
}

float blend(in float d1, in float d2, in float k) {
	return k * d1 + (1.-k) * d2;
}

float vmax(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

vec3 rot(vec3 p, mat4 m) {
    vec4 q = inverse(m)*vec4(p,1.);
    return vec3(q.x,q.y,q.z);
}

// Sinusoid bumps
float sinbumps(in vec3 p){
    float frq = 1.7;
    return sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31)*sin(time)*1.7;
}

////////////////////////////////////////////// SCENE SETUP /////////////////////////////////////////////
// precompute stuff
float speed = pow(1.01,100.);
float speed_z = 0.1+(100.*mx15);
float speed_x = speed_z * .5;

    float camera_rot = 0.;
	vec2 path2d = vec2(cos(time/3.)*(speed_x+camera_rot), time*speed_z);
	vec2 path2d_next = vec2(cos(time+2./3.)*speed_x+camera_rot/2., (time+2.)*speed_z);

	// s1 / heightmap
	vec3 path = vec3( path2d.x, 2.-noise2f(path2d.xy/10.)*10., path2d.y );
	vec3 path_next = vec3(path2d_next.x, 2.-noise2f(path2d_next.xy/10.)*10., path2d_next.y );
	int MAXSTEPS = 32;
	int MAXDIST = 20;
    float eps = 0.0001;

	//vec3 ro = path;
	//vec3 lookat = path_next;
    vec3 ro = vec3(0.,m0*10., 5.);
    //float lookat_y = 2.-noise2f(path2d_next.xy/10.)*10.;
    vec3 lookat = vec3(0.);
	vec3 light_position = vec3(path.x,0.,path.z + sin(time)*4);

float plane(vec3 p) {
    // bend
    float k = .2; // or some other amount
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2 m = mat2(c,-s,s,c);
    p = vec3(m*p.xy,p.z);
	return p.y + 10.+ noise2f(p.xz/10.)*10. + noise2f(p.xz) + noise2f(p.xz*5.);
}

float sky_plane(vec3 p) {
    // bend
    float k = -.2; // or some other amount
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2 m = mat2(c,-s,s,c);
    p = vec3(m*p.xy,p.z);
	return -p.y + 10. - noise2f(p.xz/10.)*10. + noise2f(p.xz) + noise2f(p.xz*5.);
}

float scene(vec3 p) {
    p -= path;
   	return min(plane(p), sky_plane(p));
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p) {
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*scene( p + k.xyy*eps ) +
                      k.yyx*scene( p + k.yyx*eps ) +
                      k.yxy*scene( p + k.yxy*eps ) +
                      k.xxx*scene( p + k.xxx*eps ) );
}

void main() {
	// raymarching & camera setup
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);
    l1.color = vec3(1.,1.,1.);
	l1.position = light_position;
    color = vec3(.0);
    sky = color;

    // raymarch a scene
    float step = .0;
    vec3 p_refr;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p3 = ro + (rd * step);
        float d = scene(p3);
        if (d > MAXDIST) {
            break;
        }
        if (d < eps) {
            break;
        }
		color += vec3(float(i)*(1./float(MAXSTEPS*7)),0.,0.);
        step += d;
    }

    // send to screen
    PixelColor = vec4(color+spec, 1.);
}
