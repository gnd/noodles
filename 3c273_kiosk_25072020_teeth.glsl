#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define MAXSTEPS 256
#define MAXDIST 25.0
#define eps 0.001

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

vec3 rot(vec3 p, vec3 axis, float angle) {
    axis = normalize(axis);
    angle = radians(angle);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    mat4 m = inverse(mat4(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s, oc * axis.z * axis.x + axis.y * s, 0.0,
                oc * axis.x * axis.y + axis.z * s, oc * axis.y * axis.y + c, oc * axis.y * axis.z - axis.x * s, 0.0,
                oc * axis.z * axis.x - axis.y * s, oc * axis.y * axis.z + axis.x * s, oc * axis.z * axis.z + c, 0.0,
                0.0, 0.0, 0.0, 1.0));

    vec4 q = m*vec4(p,1.);
    return vec3(q.x,q.y,q.z);
}

///////////////////////// END OF NUMERICAL FUNCTIONS ////////////////////////////////////

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material m1;
shading s1;

vec3 sky, color; // 'sky color'

// Counter 0 to 1
float cnt(int m) {
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

// Sin between 0 - 1 - bouncing
float a_sin(in float a) {
	return abs(sin(a));
}

float sinbumps(in vec3 p, in float frq) {
    return sin(p.x*frq) + sin(p.y*frq) + sin(p.z*frq);
}

float bumps(in vec3 p, in float frq) {
    //return sinbumps(p, frq);
    return sinbumps(p, frq);
}

vec3 dom_distort_iq(in vec3 q) {
        q.xyz += 1.000*sin(  2.0*q.xzy )*.2;
        q.xyz += 0.500*sin(  4.0*q.yzx )*.2;
        q.xyz += 0.250*sin(  8.0*q.yzx )*.1;
        q.xyz += 0.050*sin( 26.0*q.zyx )*.3;
        return q;
}

float smin( float a, float b, float k ) {
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

float maxcomp(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

float box(vec3 p, vec3 b) { //cheap box
	return maxcomp(abs(p) - b);
}

float capsule(vec3 p, float r, float c) {
	return mix(length(p.xz) - r, length(vec3(p.x, abs(p.y) - c, p.z)) - r, step(c, abs(p.y)));
}

float plane(vec3 p) {
    // slow (but has 'substance')
    //return box(p-vec3(.0,-1.*.5,.0), vec3(15.,.01,15.), .0);
    return p.y+.5;
}

float sphere( vec3 p, float r ) {
    return length(p) - r;
}

float cylinder(vec3 p, float r, float height) {
	float d = length(p.zy) - r;
	d = max(d, abs(p.x) - height);
	return d;
}

float bone (in vec3 p) {
    p.y*=1.7;
    p = dom_distort_iq(p);
    // joints
    float s1 = sphere(p-vec3(-1.9, 1.5, .28), 0.19);
    float s2 = sphere(p-vec3(-1.9, 1.5, -.28), 0.09);
    float s3 = sphere(p-vec3(1.9, 1.5, .28), 0.28);
    float s4 = sphere(p-vec3(1.9, 1.5, -.28), 0.09);
    // bone trunk
    float t = cylinder(p-vec3(.0, 1.5, .0), .2, 1.3);

    return smin(t, smin(smin(s1, s2, 20.), smin(s3, s4, 20.), 1.), 2.5);
}


// TODO
// make both tooths react to p0 - p4
// add a slight sound jitter on final smin
// make two more teeth, one of them walking
float tooth1(in vec3 p) {
    float mov = sin(time*10.)+1.0;
    float x = .0;
    float z = 0.;

    vec3 pp = p;
    p.x += sin(pp.y*3.)*.07;
    p.z += sin(pp.y*4.)*.04;
    float b1 = capsule(p-vec3(x-.4,.7,z+.3), .1, .6);
    p.x += sin(pp.y*4.)*.05;
    p.z += sin(pp.y*2.)*.01;
    float b2 = capsule(p-vec3(x-.4,.7,z-.3), .1,.5);
    p.x += sin(pp.y*2.)*.03;
    p.z += sin(pp.y*8.)*.01;
    float b3 = capsule(p-vec3(x+.4,.7,z+.25), .1*mov, .8);
    p.x += sin(pp.y*7.)*.05;
    p.z += sin(pp.y*6.)*.03;
    float b4 = capsule(p-vec3(x+.4,.7,z-.4), .09, .8);
    p = pp;

    p.x -= sin(p.y*mov)*.4;
    p.z += cos(p.y*5.*p.x)*.1;
    float s = sphere(p-vec3(.0,2.4,.0), .8);
    p = pp;

    return smin(s, min(b1,min(b2,smin(b3,b4, 22.))), 3.);
}

float tooth2(in vec3 p) {
    float mov = sin(time*5.)+1.0;
    float x = .5;
    float z = 2.;

    vec3 pp = p;
    p.x += sin(pp.y*3.)*.07;
    p.z += sin(pp.y*4.)*.04;
    float b1 = capsule(p-vec3(x-.4,.7,z+.3), .1, .6);
    p.x += sin(pp.y*4.)*.05;
    p.z += sin(pp.y*2.)*.01;
    float b2 = capsule(p-vec3(x-.4,.7,z-.3), .1,.5);
    p.x += sin(pp.y*2.)*.03;
    p.z += sin(pp.y*8.)*.01*mov;
    float b3 = capsule(p-vec3(x+.4,2.-mov,z+.25), .1, .8*mov);
    p.x += sin(pp.y*7.)*.05;
    p.z += sin(pp.y*6.)*.03;
    float b4 = capsule(p-vec3(x+.4,1.7-mov/2.,z-.4), .09, .8*mov);
    p = pp;

    p.x -= sin(p.y*mov)*.1;
    p.z += cos(p.y*7.*p.x)*.1;
    float s = sphere(p-vec3(x,2.4,z), .8);
    p = pp;

    //return smin(s, smin(b, min(b1,min(b2,smin(b3,b4, 22.))), 30.), 5.);
    return smin(s, min(b1,min(b2,smin(b3,b4, 22.))), 3.);
}

float scene(vec3 p) {
    return min(plane(p), min(tooth1(p), tooth2(p)));
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p) {
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*scene( p + k.xyy*eps ) +
                      k.yyx*scene( p + k.yyx*eps ) +
                      k.yxy*scene( p + k.yxy*eps ) +
                      k.xxx*scene( p + k.xxx*eps ) );
}

// from https://alaingalvan.tumblr.com/post/79864187609/glsl-color-correction-shaders
vec3 brightcon(vec3 c, float brightness, float contrast) {
    return (c - 0.5) * contrast + 0.5 + brightness;
}

shading get_shading(material m, light l, vec3 p, vec3 n, vec3 ld, vec3 ed) {
    shading s;
    float step;
    s.diffuse = clamp(dot(n,ld), 0., 1.);
    s.specular = clamp(m.reflection_ratio * pow(dot(normalize(reflect(-ld, n)), ed), m.shininess), 0., 1.);

    s.shadow = 1.;
    step = 0.01;
    float ph = 1e10;
    s.shadow = clamp(s.shadow, 0., 1.);
    s.aoc = 1.;
    s.amb = clamp(0.4*n.y+0.5, 0., 1.);

    return s;
}

void main() {
    float orbit = 5.;
    vec3 ro = vec3(cos(time/4.)*orbit, 2., sin(time/4.)*orbit); // ray origin, here also known as 'eye'
    //vec3 ro = vec3(1.,1.8,1.);
    vec3 lookat = vec3(-0.9, 1.4,2.);
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    l1.position = vec3(4., 2.5, 3.);
    color = vec3(.0);
    sky = color;

    // raymarch a scene
    float step = 1.;
    float d = 1.;
    vec3 p;
    for (int i = 0; i < MAXSTEPS; i++) {
        if( (d<eps) || (step > MAXDIST) ) break;
        p = ro+rd*step;
        d = scene( p );
        step += d;
    }
    if( step>MAXDIST ) step=-1.;

    // apply materials
    if (step > 0.) {
        // precompute stuff
        vec3 n = normal(p);
        vec3 ld = normalize(l1.position-p);
        vec3 ed = normalize(ro-p);

        if (plane(p) < eps) {
            float checker = mod(floor(p.x)+floor(p.z)-.5, 2.0);
            m1.color = vec3(checker);
            m1.reflection_ratio = 0.01;
            m1.shininess = 3.;
            s1 = get_shading(m1, l1, p, n, ld, ed);
            color = m1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += m1.color * s1.amb *.2;
        }

        if (tooth1(p) < eps) {
            m1.color = vec3(-1.,0.0,.0);
            m1.color += vec3(tan((p+cnt(1000))*5.).x,.0,.0);
            m1.reflection_ratio = 7.5;
            m1.shininess = 20.1;
            s1 = get_shading(m1, l1, p, n, ld, ed);
            color = m1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color += m1.color * s1.amb *.2;
        }

        if (tooth2(p) < eps) {
            m1.color = vec3(.9);
            m1.reflection_ratio = 0.01;
            m1.shininess = 3.;
            s1 = get_shading(m1, l1, p, n, ld, ed);
            color = m1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += m1.color * s1.amb *.2;
        }
    }

    // Exponential distance fog
    color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * step * step));

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(color, 1.);
}
