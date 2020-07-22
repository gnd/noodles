#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
uniform vec2 rand;
uniform sampler2D backbuffer;
#define MAXSTEPS 256
#define MAXDIST 25.0
#define eps 0.001

/////////////////// MIDIMIX //////////////////////
uniform float m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10;
uniform float m11,m12,m13,m14,m15,m16,m17,m18,m19,m20;
uniform float m21,m22,m23,m24,m25,m26,m27,m28,m29,m30;
uniform float m31,m32,m33,m34,m35,m36,m37,m38,m39,m40;
uniform float m41,m42,m43,m44,m45,m46,m47,m48,m49,m50;
uniform float m51,m52,m53,m54,m55,m56;

float mx11 = m5/127.;
float mx12 = m6/127.;
float mx13 = m7/127.;
float mx14 = m8/127.;
float mx15 = m9/127.;
float mx16 = m10/127.;
float mx21 = m11/127.;
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
///////////////// NANOKONTROL //////////////////////
uniform float m57,m58,m59,m60,m61,m62,m63,m64,m65,m66,m67,m68;
uniform float m69,m70,m71,m72,m73,m74,m75,m76;
uniform float m77,m78,m79,m80,m81,m82,m83,m84;

float nkk1 = m57/127;
float nkk2 = m58/127;
float nkk3 = m59/127;
float nkk4 = m60/127;
float nkk5 = m61/127;
float nkk6 = m62/127;
float nkk7 = m63/127;
float nkk8 = m64/127;
float nkk9 = m65/127;
float nkk10 = m66/127;
float nkk11 = m67/127;
float nkk12 = m68/127;
float nk11 = m69/127;
float nk12 = m70/127;
float nk21 = m71/127;
float nk22 = m72/127;
float nk31 = m73/127;
float nk32 = m74/127;
float nk41 = m75/127;
float nk42 = m76/127;
float nk51 = m77/127;
float nk52 = m78/127;
float nk61 = m79/127;
float nk62 = m80/127;
float nk71 = m81/127;
float nk72 = m82/127;
float nk81 = m83/127;
float nk82 = m84/127;
///////////////// END OF NANOKONTROL ////////////////////

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

// Counter 0 to 1
float cnt(int m) {
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

// Sin between 0 - 1 - bouncing
float a_sin(in float a) {
	return abs(sin(a));
}

float smin( float a, float b, float k ) {
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

float maxcomp(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

///////////////////////// END OF NUMERICAL FUNCTIONS ////////////////////////////////////

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };
vec3 sky, color, glow;
vec2 p2;

light l1;
material mat,mt1,mt2,mt3,mt4,mt5;
shading s1,s2,s3;

// precompute various rhytms
float mov0 = sin(time)+1.0;
float mov1 = sin(time*5.)+1.0;
float mov2 = sin(time/2.)+1.0;
float mov3 = sin(time*20.)+1.0;
float mov4 = sin(time*50.)+1.0;

float sinbumps(in vec3 p, in float frq) {
    return sin(p.x*frq) + sin(p.y*frq) + sin(p.z*frq);
}

float cosbumps(in vec3 p, in float frq) {
    return cos(p.x*frq) + cos(p.y*frq) + cos(p.z*frq);
}

vec3 dom_distort_iq(in vec3 q) {
        q.xyz += 1.000*sin(  2.0*q.xzy )*.2;
        q.xyz += 0.500*sin(  4.0*q.yzx )*.2;
        q.xyz += 0.250*sin(  8.0*q.yzx )*.1;
        q.xyz += 0.050*sin( 26.0*q.zyx )*.3;
        return q;
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
float tooth1(in vec3 p) {
    float mov = sin(time*10.)+1.0;
    float x = 0;
    float z = 0.;
    float y = sinbumps(p-vec3(time*10.,.0,.0), 0.25)+.25;

    vec3 pp = p;
    pp.y += y;

    p.x += sin(pp.y*3.)*.07;
    p.z += sin(pp.y*4.)*.04;
    float b1 = capsule(p-vec3(x-.4,.7-y,z+.3), .1, .6);
    p.x += sin(pp.y*4.)*.05;
    p.z += sin(pp.y*2.)*.01;
    float b2 = capsule(p-vec3(x-.4,.7-y,z-.3), .1,.5);
    p.x += sin(pp.y*2.)*.03;
    p.z += sin(pp.y*8.)*.01;
    float b3 = capsule(p-vec3(x+.4,.7-y,z+.25), .1*mov, .8);
    p.x += sin(pp.y*7.)*.05;
    p.z += sin(pp.y*6.)*.03;
    float b4 = capsule(p-vec3(x+.4,.7-y,z-.4), .09, .8);
    p = pp;

    p.x -= sin(p.y*mov)*.4;
    p.z += cos(p.y*5.*p.x)*.1;
    float s = sphere(p-vec3(x,2.4,z), .8);
    p = pp;

    return smin(s, min(b1,min(b2,smin(b3,b4, 22.))), 3.);
}

float tooth2(in vec3 p) {
    float x = -2.;
    float y = -4.;
    float z = .0;

    float mov = mov1;
    p = rot(p, vec3(0.,-0.05,1.), 120);
    vec3 pp = p;

    p.x += sin(pp.y*3.)*.07;
    p.z += sin(pp.y*4.)*.04;
    float b1 = capsule(p-vec3(x-.4,y+.7,z+.3), .1, .6);
    p.x += sin(pp.y*4.)*.05;
    p.z += sin(pp.y*2.)*.01;
    float b2 = capsule(p-vec3(x-.4,y+.7,z-.3), .1,.5);
    p.x += sin(pp.y*2.)*.03;
    p.z += sin(pp.y*8.)*.01*mov;
    float b3 = capsule(p-vec3(x+.4,y+2.-mov,z+.25), .1, .8*mov);
    p.x += sin(pp.y*7.)*.05;
    p.z += sin(pp.y*6.)*.03;
    float b4 = capsule(p-vec3(x+.4,y+1.7-mov/2.,z-.4), .09, .8*mov);
    p = pp;

    p.x -= sin(p.y*mov)*.1;
    p.z += cos(p.y*7.*p.x)*.1;
    float s = sphere(p-vec3(x,y+2.4,z), .8);
    p = pp;

    //return smin(s, smin(b, min(b1,min(b2,smin(b3,b4, 22.))), 30.), 5.);
    return smin(s, min(b1,min(b2,smin(b3,b4, 22.))), 3.);
}

float tooth3(in vec3 p) {
    float x = 2.;
    float y = -4.;
    float z = 0.;
    float mov = mov1;
    p = rot(p, vec3(0.,0.1,1.), -120);
    vec3 pp = p;

    p.x += sin(pp.y*3.)*.07;
    p.z += sin(pp.y*4.)*.1*mov0;
    float b1 = capsule(p-vec3(x-.4,y+.7-5.*mov0,z+.3), .1, .6+5.*mov0);
    p.x += sin(pp.y*4.)*.05;
    p.z += sin(pp.y*2.)*.01;
    float b2 = capsule(p-vec3(x-.4,y+.7,z-.3), .1,.5);
    p.x += sin(pp.y*2.)*.03;
    p.z += sin(pp.y*8.)*.01*mov;
    float b3 = capsule(p-vec3(x+.4,y+.9,z+.25), .1, .8);
    p.x += sin(pp.y*7.)*.05;
    p.z += sin(pp.y*6.)*.03;
    float b4 = capsule(p-vec3(x+.4,y+.9,z-.4), .09, .8);
    p = pp;

    p.x -= sin(p.y)*mov/2.-.25;
    p.z += cos(p.y*7.*p.x*mov1/2.)*.1;
    float s = sphere(p-vec3(x,y+2.4-mov/4.,z), .8);
    p = pp;

    //return smin(s, smin(b, min(b1,min(b2,smin(b3,b4, 22.))), 30.), 5.);
    return smin(s, min(b1,min(b2,smin(b3,b4, 10.))), 3.);
}

float org1(in vec3 p) {
    p = p + vec3(0.,22.*mx15-5.,.0);
    return smin(tooth1(p), smin(tooth2(p), tooth3(p), 5*mx16), 5*mx16);
}

float org2(in vec3 p) {
    float x = 4.;
    float y = 6.;
    float z = 0.;
    vec3 pp = p;

    p = rot(pp, vec3(-.7,0.5,1.), 360. * m1);
    p.x += sin(pp.y*30.*m3)*.09;
    p.z += sin(pp.y*4.)*.1*mov0;
    float b1 = capsule(p-vec3(x-.4,y+.7-5.*mov0,z+.3), .1, (.6+5.*mov0)*m0);

    p = rot(pp, vec3(1.,0.5,0.), 360. * m1);
    p.x += sin(pp.y*30.*m3)*.09;
    p.z += sin(pp.y*4.)*.1*mov0;
    float b2 = capsule(p-vec3(x-.8,y-3-2.*mov1,z+.3), .1+.1*m2, (.5+2.*mov1)*m0);

    p = rot(pp, vec3(1.,0.5,0.), 180. * m2);
    p.x += sin(pp.y*7.)*.06;
    p.z += sin(pp.y*40.*m2)*.1*mov0;
    float b3 = capsule(p-vec3(x-2.,y-3-2.*mov0,z+.3), .1+.1*m2, (.6+2.*mov0)*m0);

    p = rot(pp, vec3(0.,1.,0.2), 360. * m3);
    p.x += sin(pp.y*20.)*.06;
    p.z += sin(pp.y*15.*m1)*.1*mov0;
    float s1 = capsule(p-vec3(x-2.,y-3-5.*mov2,z-3.), .1+.1*m2, (.6+5.*mov0)*m0);

    return smin(smin(min(b1, b2),b3, 1.*mx14), s1, 1.*mx14);
}

float scene(vec3 p) {
    //return min(tooth1(p), min (tooth2(p), tooth3(p)));
    return org2(p);
    //return min(org1(p), org2(p));

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

vec3 str(in float xp, in float yp, in float x, in float y, in vec3 c, in vec3 col) {
        vec3 cl = ((p2.x > xp-x*.51) && (p2.x < x*.51+xp) && (p2.y > yp-y*.51) && (p2.y < y*.51+yp)) ? -c+col : vec3(0.0);
        return cl;
}

// Feedback
vec3 feedb_sqr(in float xpos, in float ypos, in float xsiz, in float ysiz, in float bsiz, in vec3 c) {
    vec3 ccc = texture2D(backbuffer, (vec2(p2.x,p2.y)-0.51)*bsiz*1.5+0.51).xyz;
    c+=str(xpos,ypos,xsiz,ysiz,c,ccc);
	return c;
}

vec3 stripe(in vec3 p, in float xp, in float yp, in float x, in float y, in vec3 c, in vec3 col) {
    vec3 cl = ((p.x > xp-x*.51) && (p.x < x*.51+xp) && (p.y > yp-y*.51) && (p.y < y*.51+yp)) ? -c+col : vec3(0.0);
    return cl;
}

// new feedback
vec3 int_feedback(in vec3 p, in float xpos, in float ypos, in float xsiz, in float ysiz, in float bsiz, in vec3 c) {
    vec3 texture = texture2D(backbuffer, (vec2(p.x,p.y)-0.51)*bsiz*1.5+0.51).xyz;
    c += stripe(p, xpos, ypos, xsiz, ysiz, c, texture);
	return c;
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
    float orbit = 7.;
    vec3 ro = vec3(cos(time/4.)*orbit, 2., sin(time/4.)*orbit); // ray origin, here also known as 'eye'
    //ro = vec3(0.,1.8,-8.);
    vec3 lookat = vec3(-0.9, 0.,0.);
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    color = vec3(.0);
    sky = color;
    glow = color;

    // raymarch a scene
    float step = 1.;
    float d = 1.;
    vec3 p;
    for (int i = 0; i < MAXSTEPS; i++) {
        if( (d<eps) || (step > MAXDIST) ) break;
        p = ro+rd*step;
        d = scene( p );
        step += d;

        // red glow
        glow = vec3(float(i)/float(MAXSTEPS),0.,0.);
    }
    if( step>MAXDIST ) step=-1.;

    // apply materials
    if (step > 0.) {
        // precompute stuff
        vec3 n = normal(p);
        vec3 ld = normalize(l1.position-p);
        vec3 ed = normalize(ro-p);

        // setup materials
        // tooth1
        mt1.color = vec3(-1.,0.0,.0);
        mt1.color += vec3(tan((p+cnt(1000))*5.).x,.0,.0);
        mt1.reflection_ratio = 7.5;
        mt1.shininess = 20.1;
        // tooth2
        mt2.color = vec3(.9) + glow*20.;
        mt2.reflection_ratio = 0.01;
        mt2.shininess = 3.;
        // tooth3
        mt3.color = glow*50.;
        mt3.reflection_ratio = 1.01;
        mt3.shininess = 3.;
        // iridescent1
        vec3 pp = mod(p*100., 2.) - .5;
        vec3 i = normalize(ro - pp);
        vec3 ir = refract(i, n, 1./1.);
        vec3 ir2 = refract(i, -n, 1./2.);
        vec3 ir_color2 = vec3(cross(ir2, -n));
        mt4.color = ir_color2 * clamp(dot(n,ld), 0., 1.);
        mt4.color += ir_color2*1.3 + clamp(2.9 * pow(dot(normalize(reflect(-ld, n)), ed), 1.9), 0., 1.) * .1;
        mt4.color += clamp(.5 * pow(dot(normalize(reflect(-ld, n)), ed), .1), 0., 1.) * 0.1;
        mt4.color += vec3(cross(ir, n)) * clamp(.9 * pow(dot(normalize(reflect(-ld, n)), ed), .9), 0., 1.) * 2.5;
        mt4.color += ir_color2 * clamp(0.4*n.y+0.5, 0., 1.) * 1.1;
        mt4.reflection_ratio = .01;
        mt4.shininess = 0.;

        if (org2(p) < eps) {
            mat = mt3;
            s1 = get_shading(mat, l1, p, n, ld, ed);
            color = mat.color * s1.diffuse * s1.shadow/.5;
            color += s1.specular;
            color += mat.color * s1.amb *.2;
        }

        if (org1(p) < eps*20.) {
            if (mx12 == 0.) {
                mat = mt1;
            }
            if (mx12 == 1.) {
                mat = mt2;
            }
            if (mx12 == 2.) {
                mat = mt3;
            }
            if (mx12 == 3.) {
                mat = mt4;
            }
            s1 = get_shading(mat, l1, p, n, ld, ed);
            color = mat.color * s1.diffuse * s1.shadow/.5;
            color += s1.specular;
            color += mat.color * s1.amb *.2;
        }

        if (tooth1(p) < eps*20.) {
            if (mx13 < 1.) {
                mat = mt1;
            }
            s1 = get_shading(mat, l1, p, n, ld, ed);
            color = mat.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color += mat.color * s1.amb *.2;
        }

        if (tooth2(p) < eps*20.) {
            if (mx13 < 1.) {
                mat = mt2;
            }
            s1 = get_shading(mat, l1, p, n, ld, ed);
            color = mat.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color += mat.color * s1.amb *.2;
        }

        if (tooth3(p) < eps*20.) {
            if (mx13 < 1.) {
                mat = mt3;
            }
            s1 = get_shading(mat, l1, p, n, ld, ed);
            color = mat.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color += mat.color * s1.amb *.2;
        }
    }

    // Exponential distance fog
    //color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * step * step));

    // feedback normal
    float asp = resolution.x / resolution.y;
    p2 = gl_FragCoord.xy / resolution;
    if (nkk11 > 0.) {
        color += feedb_sqr(.5, .5, 1., 1., mx74, color)*mx71*1.1;
    }

    // feedback internal
    if (nkk12 > 0.) {
        color -= int_feedback(p, .0, .0, 20., 20., mx76*.05, color)*mx75*10.;
    }

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // channel mix
	color = vec3(color.r*mx86,color.g*mx85,color.b*mx84);

    // strobe
    if (nkk1 > 0.) {
        color += cnt(100);
    }

    // send to screen
    PixelColor = vec4(color, 1.);
}
