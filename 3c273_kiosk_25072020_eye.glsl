#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define MAXSTEPS 128
#define MAXDIST 30.0
#define eps 0.0001
#define PI 3.14159265

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material m1;
shading s1;

vec3 sky, color; // 'sky color'
//vec3 ro = vec3(3., 2., 3.);
vec3 ro = vec3(cos(time/4.)*3., 2., sin(time/4.) *3.);
vec3 lookat = vec3(0.,0.,0.);
vec3 fwd = normalize(lookat-ro);
vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
vec3 up = normalize(cross(fwd, right));
vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
float aspect = resolution.x/resolution.y;
vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

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

float sphere( vec3 p, float r ) {
    return length(p) - r;
}

float maxcomp(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

vec3 opRep( in vec3 p, in vec3 c )
{
    return mod(p+0.5*c,c)-0.5*c;
}

float box(vec3 p, vec3 b) { //cheap box
	return maxcomp(abs(p) - b);
}

float bbox(vec3 p, vec3 b) {
	vec3 d = abs(p) - b;
	return length(max(d, vec3(0))) + maxcomp(min(d, vec3(0)));
}

float cbox(in vec3 p) {
    p.x = p.x * sin(p.x)*.4;
    p = rot(p, vec3(0.,1.,0.), mod(time*50., 360));
    p = opRep(p, vec3(1., 0., .3));
    p = rot(p, vec3(1.,1.,.0), mod(time*50., 360));
    return bbox(p - vec3(0.,.9,.0), vec3(.1,.7,.2));
}

float plane(vec3 p) {
    return p.y;
}

float sclera(vec3 p, float r) {
    return sphere(p, r);
}

float iris(vec3 p, float r) {
    return sphere(p, r);
}

float pupil(vec3 p, float r) {
    return sphere(p, r);
}

float eye(vec3 p) {
    vec3 center = vec3(0.,.9,.0);
    float sclera_dst = sclera(p - center, 1.);
    vec3 iris_c = normalize(ro) + vec3(.0,.65,.0);
    float iris_dst = iris(p-iris_c, .4);
    float pupil_dst = pupil(p - (ro)/2. - center/2., .1);
    return min(pupil_dst, min(sclera_dst, iris_dst));
}

float scene(vec3 p) {
    return eye(p);
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
    for (int i = 0; i < MAXSTEPS/3; i++ ) {
        vec3 shadow_p = ld * step + p;
        float shadow_d = scene(shadow_p);
        if (shadow_d < eps) {
            s.shadow = 0.;
            break;
        }
        if (step > MAXDIST) {
            break;
        }
        float y = pow(shadow_d, 2)/(2.0*ph);
        float d = sqrt(shadow_d*shadow_d-y*y);
        s.shadow = min( s.shadow, 32.*d/max(.0,step-y));
        ph = shadow_d;
        step += shadow_d;
    }
    s.shadow = clamp(s.shadow, 0., 1.);
    s.aoc = 1.;
    s.amb = clamp(0.4*n.y+0.5, 0., 1.);

    return s;
}


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

float eyp(vec2 uv) {
	return noise2f( vec2(20.0*uv/pow(length(uv),.9) ) );
}

float sinbumps(in vec3 p){
    float frq = 1.7;
    return sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31)*sin(time)*1.7;
}

void main() {
    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    color = vec3(.2);
    sky = color;

    // raymarch a scene
    float ray_step = .1;
    float d = .1;
    vec3 p;
    for (int i = 0; i < MAXSTEPS; i++) {
        if ((d<eps) || (ray_step > MAXDIST)) break;
        p = ro+rd*ray_step;
        d = scene( p );
        ray_step += d;
    }
    if( ray_step>MAXDIST ) ray_step=-1.;

    // apply materials
    if (ray_step > .0) {
        // precompute stuff
        vec3 n = normal(p);
        vec3 ld = normalize(l1.position-p);
        vec3 ed = normalize(ro-p);

        if (sclera(p - vec3(0.,.9,.0), 1.) < eps) {
            m1.color = vec3(1.)*2.;
            //m1.color = vec3(abs(sin(p.x*10.+cos(p.z*20.)))-.01);
            //vec3 lol = ro - rd;
            //m1.color += eyp(p.xz*1000.*fwd.xx)*20.;

            //m1.color *= sinbumps(m1.color);
            //m1.color = vec2( length(m1.color.xz) - .01, m1.color.y).xyy;
            m1.reflection_ratio = 10.01;
            m1.shininess = 30.;
            s1 = get_shading(m1, l1, p, n, ld, ed);
            color = m1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += m1.color * s1.amb *.2;
        }

        if (iris(p- normalize(ro) - vec3(.0,.65,.0), .4) < eps) {
            m1.color = vec3(10.1, 0.2, 2.9);
            m1.color *= eyp((uv+vec2(0.,-.35))*100.)*5.;
            m1.color *= eyp((uv+vec2(0.,-.35))*10.)*15.;


            m1.reflection_ratio = 10.01;
            m1.shininess = 350.;
            s1 = get_shading(m1, l1, p, n, ld, ed);
            color = m1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += m1.color * s1.amb *.2;
        }

        if (pupil(p - (ro-lookat)/2. - vec3(0.,.45,.0), .1) < eps) {
            color = vec3(0.);
        }


    }

    // Exponential distance fog
    color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * ray_step * ray_step));

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(color, 1.);
}