#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define eps 0.00001
#define MAXSTEPS 128
#define MAXDIST 12

// http://glslsandbox.com/e#39898.3

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material m1,m2,m3;
shading s1,s2,s3;

vec3 sky, color; // 'sky color'

// Sinusoid bumps
float sinbumps(in vec3 p, in float frq) {
    return sin(p.x*frq) * sin(p.y*frq+time*4.) * sin(p.z*frq) * .17;
}

float smin( float a, float b, float k ) {
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

// Sphere field
float sphere(in vec3 p, in vec3 centerPos, float radius) {
        return length(p-centerPos) - radius;
}

float ball(vec3 p, vec3 pos, float bumps) {
   return sphere(p, pos, 0.1) + bumps;
}

float ball2(vec3 p, vec3 pos, float bumps) {
   return sphere(p, pos, 0.1) + bumps;
}

float structure(vec3 p) {
    float d = 0;
    vec3 a;
    float s = .2;
    for (int j = 0; j < 6; j++)
        a = max(a = abs(mod(p*s+.3,.97)-.9), a.zxy),
        d = max(d,(0.22-length(a*0.51)*.35)/s),
        s *= 3.2;

    return d;
}

float scene(vec3 p, vec3 pos) {
    float bumps = sinbumps(p, 10.);
    return smin (
            smin (
                structure(p),
                ball(p, pos, bumps),
                1.7),
            ball2(p, pos - vec3(3.,0.,3.), bumps),
            1.7
        );
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p, vec3 pos) {
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*scene( k.xyy*eps + p , pos) +
                      k.yyx*scene( k.yyx*eps + p , pos) +
                      k.yxy*scene( k.yxy*eps + p , pos) +
                      k.xxx*scene( k.xxx*eps + p , pos) );
}

// from https://www.shadertoy.com/view/llcXWM
vec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ) {
    return b*cos( 6.28318*(c*t+d) ) + a;
}

vec3 spectrum(float n) {
    return pal( n, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );
}

shading get_shading(material m, light l, vec3 p, vec3 n, vec3 ld, vec3 ed, vec3 pos) {
    shading s;
    float step;
    s.diffuse = clamp(dot(n,ld), 0., 1.);
    s.specular = clamp(m.reflection_ratio * pow(dot(normalize(reflect(-ld, n)), ed), m.shininess), 0., 1.);

    s.shadow = 1.;
    step = 0.01;
    float ph = 1e10;
    for (int i = 0; i < MAXSTEPS/3; i++ ) {
        vec3 shadow_p = ld * step + p;
        float shadow_d = scene(shadow_p, pos);
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

void main() {

    // classic scene instantiate
    vec3 lookat = vec3(time*2.,0.,0.);
    vec3 ro = lookat + vec3(cos(time/4.)*3., sin(time/5)*2.5, sin(time/4.)*3.);
    
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.1*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = lookat + vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    sky = vec3(.0);
    color = sky;

    vec3 p, pp;
    float step = .0;
    
    // RAYMARCH A SCENE
    pp = p;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = rd * step + ro;
    
        float d = scene(p, lookat);

        if (d > MAXDIST) {
            break;
        }

        step += d;
        
        if (d < eps) {
            // precompute stuff
            vec3 n = normal(p, lookat);
            vec3 ld = normalize(l1.position-p);
            vec3 ed = normalize(ro-p);
            float bumps = sinbumps(p, 10.);

            PixelColor += vec4(0.005 * (MAXSTEPS-i));
            break;
        }

    }
}
