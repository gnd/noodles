#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define MAXSTEPS 128
#define MAXDIST 20.0
#define eps 0.0001

/////////////////// MIDIMIX //////////////////////
/// PARAMS
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;
uniform float m5;
uniform float m6;
uniform float m7;

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

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material mat1;
shading s1;

vec3 sky, color; // 'sky color'

float sphere( vec3 p, float r ) {
    return length(p) - r;
}

float maxcomp(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

float box(vec3 p, vec3 b) { //cheap box
	return maxcomp(abs(p) - b);
}

float plane(vec3 p) {
    return p.y;
}

float box1(in vec3 p) {
    return box(p-vec3(-1.,1.,0.), vec3(.1,m0,.1));
}

float box2(in vec3 p) {
    return box(p-vec3(-.5,1.,0.), vec3(.1,m1,.1));
}

float box3(in vec3 p) {
    return box(p-vec3(0.,1.,0.), vec3(.1,m2,.1));
}

float box4(in vec3 p) {
    return box(p-vec3(.5,1.,0.), vec3(.1,m3,.1));
}

float box5(in vec3 p) {
    return box(p-vec3(1.,1.,0.), vec3(.1,m4,.1));
}

float scene(vec3 p) {
    float b1 = box1(p);
    float b2 = box2(p);
    float b3 = box3(p);
    float b4 = box4(p);
    float b5 = box5(p);

    return min(min(b1, min(b2, min(b3, min(b4,b5)))), plane(p));
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

void main() {
    vec3 ro = vec3(cos(time/4.)*3., 2., sin(time/4.)*3.); // ray origin, here also known as 'eye'
    vec3 lookat = vec3(0.,0.,0.);
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    //l1.position = vec3(1.2, 2.0, -1.9);
    color = vec3(.0);
    sky = color;

    // raymarch a scene
    float step = 1.;
    float d = 1.;
    vec3 p;
    for (int i = 0; i < MAXSTEPS; i++) {
        if ((d<eps) || (step > MAXDIST)) break;
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
            mat1.color = vec3(checker);
            mat1.reflection_ratio = 0.01;
            mat1.shininess = 3.;
            s1 = get_shading(mat1, l1, p, n, ld, ed);
            color = mat1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += mat1.color * s1.amb *.2;
        }

        if (box1(p) < eps) {
            mat1.color = vec3(mx16, mx15, mx14);
            mat1.reflection_ratio = 0.01;
            mat1.shininess = 3.;
            s1 = get_shading(mat1, l1, p, n, ld, ed);
            color = mat1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += mat1.color * s1.amb *.2;
        }

        if (box2(p) < eps) {
            mat1.color = vec3(mx26, mx25, mx24);
            mat1.reflection_ratio = 0.01;
            mat1.shininess = 3.;
            s1 = get_shading(mat1, l1, p, n, ld, ed);
            color = mat1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += mat1.color * s1.amb *.2;
        }

        if (box3(p) < eps) {
            mat1.color = vec3(mx36, mx35, mx34);
            mat1.reflection_ratio = 0.01;
            mat1.shininess = 3.;
            s1 = get_shading(mat1, l1, p, n, ld, ed);
            color = mat1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += mat1.color * s1.amb *.2;
        }

        if (box4(p) < eps) {
            mat1.color = vec3(mx46, mx45, mx44);
            mat1.reflection_ratio = 0.01;
            mat1.shininess = 3.;
            s1 = get_shading(mat1, l1, p, n, ld, ed);
            color = mat1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += mat1.color * s1.amb *.2;
        }

        if (box5(p) < eps) {
            mat1.color = vec3(mx56, mx55, mx54);
            mat1.reflection_ratio = 0.01;
            mat1.shininess = 3.;
            s1 = get_shading(mat1, l1, p, n, ld, ed);
            color = mat1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += mat1.color * s1.amb *.2;
        }
    }

    // Exponential distance fog
    color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * step * step));

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(color, 1.);
}
