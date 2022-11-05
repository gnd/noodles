#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
uniform sampler2D backbuffer;
#define MAXSTEPS 128
#define MAXDIST 20.0
#define eps 0.0001
//#define AA // comment out to turn of anti-aliasing

/* with AA one frame takes cca 190ms on geforce 940M
   without AA its cca 50ms */

struct march { float step; vec3 ro; vec3 rd; };
struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

// PARAMS
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;

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


bool reflected = false;
bool refracted = false;

light l1;
material mt1,mt2;
shading s1,s2;

vec3 color; // 'sky color'

float sphere( vec3 p, float r ) {
    return length(p) - r;
}

float maxcomp(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

float box(vec3 p, vec3 b, float r) {
    vec3 d = abs(p) - b;
    return min(maxcomp(d),0.0) - r + length(max(d,0.0));
}

float neon_sphere(vec3 p) {
    return sphere(p - vec3(0., 2., 0.), .4 + mx11);
}

float glass_sphere(vec3 p) {
    return sphere(p-vec3(0.,1.,.0), .5);
    //return box(p-vec3(0.,1.,0.), vec3(.5,.5,.5), 0.);
}

float glass_sphere_lil(vec3 p) {
    return sphere(p-vec3(0.,1.,-1), .3);
}

float red_cube(vec3 p) {
    return box(p-vec3(1.,.4,1.5), vec3(.7,.7,.7), 0.);
}

float blue_cube(vec3 p) {
    return box(p-vec3(3.,0.,-6.), vec3(.7,.7,.7), 0.);
}

float green_cube(vec3 p) {
    return box(p-vec3(-5.,-.1,4.), vec3(.4,.4,.4), 0.);
}

float mirror(vec3 p) {
    float height = .45;
    return box(p-vec3(.0,-1.*height,.0), vec3(1.7,.15,1.7), .0);
}

float plane(vec3 p) {
    float height = .5;
    float width = 15.;
    return box(p-vec3(.0,-1.*height,.0), vec3(width,.01,width), .0);
}

float scene(vec3 p) {
    return min(min( min(glass_sphere_lil(p), glass_sphere(p)), min(mirror(p), plane(p))), min(neon_sphere(p), min(min(blue_cube(p), green_cube(p)), red_cube(p))));
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p) {
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*scene( p + k.xyy*eps ) +
                      k.yyx*scene( p + k.yyx*eps ) +
                      k.yxy*scene( p + k.yxy*eps ) +
                      k.xxx*scene( p + k.xxx*eps ) );
}

march refraction(vec3 p, vec3 rd, float eta) {
    // Glass becomes a 'portal' through which refracted rays keep on raymarching
    march march;
    vec3 n = normal(p);
    march.ro = p;
    march.rd = normalize(refract(rd, n, eta));
    march.step = 0.0;
    for (int j=0; j < MAXSTEPS; j++) {
        p = march.ro + march.rd * march.step;
        float d = scene(p);
        march.step += max(abs(d),eps);
        if (d > eps) {
            // attenuate due to impurities, etc
            color *= 0.9;
            // second refraction
            march.ro = p;
            march.rd = normalize(refract(march.rd, -normal(p), 1./eta));
            break;
        }
    }
    return march;
}

march reflection(vec3 p, vec3 light_dir) {
    march march;
    vec3 n = normal(p);
    march.ro = p;
    march.rd = normalize(reflect(-light_dir, n));
    march.step = 0.0025;
    for (int j=0; j < MAXSTEPS; j++) {
        p = march.ro + march.rd * march.step;
        float d = scene(p);
        march.step += max(abs(d),eps);
        if (d > eps) {
            break;
        }
    }
    return march;
}

shading get_shading(material m, light l, vec3 p, vec3 eye) {
    shading s;
    float step;
    vec3 light_dir = l.position - p; // direction from p to light
    vec3 normal = normal(p); // surface normal at point p

    // Diffuse shading using Lambertian reflectance
    // https://en.wikipedia.org/wiki/Lambertian_reflectance
    s.diffuse = clamp(dot(normal,normalize(light_dir)), 0., 1.);

    // Specular reflections using phong
    // https://en.wikipedia.org/wiki/Phong_reflection_model
    vec3 reflection_dir = normalize(reflect(-light_dir, normal));
    s.specular = clamp(m.reflection_ratio * pow(dot(reflection_dir, normalize(eye-p)), m.shininess), 0., 1.);

    // Shadows - basically raymarch towards the light
    // If the path from p to each light source is obstructed add shadow
    // soft shadows taken from http://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
    s.shadow = 1.;
    step = 0.01;
    float ph = 1e10;
    for (int i = 0; i < MAXSTEPS/2; i++ ) {
        vec3 shadow_p = p + normalize(light_dir) * step;
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
        s.shadow = min(s.shadow, 32.*d/max(.0,step-y));
        ph = shadow_d;
        step += shadow_d;
    }
    s.shadow = clamp(s.shadow, 0., 1.);

    /*/ broken, but could be salvaged for fake caustics
    s.shadow = 1.;
    step = 0.001;
    //float ph = 1e10;
    for (int i = 0; i < MAXSTEPS; i++ ) {
        vec3 shadow_p = p + normalize(light_dir) * step;
        float shadow_d = scene(shadow_p);
        if ((shadow_d < eps) || (step > MAXDIST)) {
            break;
        }
        float y = pow(shadow_d, 1)/(2.0*ph);
        float d = sqrt(shadow_d*shadow_d-y*y);
        s.shadow = min( s.shadow, 64.*shadow_d/max(.0,step-y));
        ph = shadow_d;
        step += shadow_d;
    }
    s.shadow = clamp(s.shadow, 0., 1.);
    */

    // Ambient occlusion
    // Raymarch along normal a few steps, if hit a surface add shadow
    s.aoc = 1.;
    step = .0;
    for (int k =0; k < 5; k++) {
        step += .04;
        vec3 aoc_p = p + normal * step;
        float diff = 1.5*abs(length(aoc_p-p) - scene(aoc_p));
        s.aoc = min(s.aoc, s.aoc-diff);
    }

    // Ambient light
    s.amb = clamp(0.5+0.4*normal.y, 0., 1.);

    return s;
}

vec3 draw(in vec2 fragcoord) {
    vec3 eye = vec3(cos(time/4.)*3., 2., sin(time/4.)*3.);
    //eye = vec3(0.8,2., -2.);
    vec3 lookat = vec3(0.,0.,0.);
    vec3 fwd = normalize(lookat-eye);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x ));
    vec3 up = normalize(cross(fwd, right));
    //float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    //float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;
    float u = fragcoord.x * 2.0 / resolution.x - 1.0;
    float v = fragcoord.y * 2.0 / resolution.y - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 ro = eye;
    vec3 rd = normalize(1.4*fwd + u*right*aspect + v*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    //l1.position = vec3(1.2, 2.0, -1.9);
    color = vec3(1.);

    // raymarch a scene
    float step = .0;
    vec3 p, p_refr, p_refl;
    float eta;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = ro + rd * step;
        float d = scene(p);
        if (d > MAXDIST) {
            if (refracted) {
                mt1.color = vec3(0.,0.,1.);
                mt1.reflection_ratio = 10.9;
                mt1.shininess = 100.9;
                s1 = get_shading(mt1, l1, p_refr, eye);
                color += s1.specular*10.;
            }
            if (reflected) {
                // attenuate sky color
                color *= 0.8;
            }
            break;
        }
        if (d < eps) {
            if ((glass_sphere(p) == scene(p)) || (glass_sphere_lil(p) == scene(p)) || mirror(p) == scene(p)) {
                if (mirror(p) == scene(p)) {
                    // Compute reflection
                    p_refl = p;
                    reflected = true;
                    vec3 ld = normalize(eye - p);
                    march march = reflection(p, ld);
                    step = march.step;
                    rd = march.rd;
                    ro = march.ro;
                } else {
                    // Compute refraction
                    p_refr = p;
                    refracted = true;
                    eta = 1./2.22;
                    march march = refraction(p, rd, eta);
                    step = march.step;
                    rd = march.rd;
                    ro = march.ro;
                }
            } else {
                if (neon_sphere(p) == scene(p)) {
                    mt1.color = vec3(0.,1.,.0);
                    mt1.reflection_ratio = .9;
                    mt1.shininess = 15.9;
                    s1 = get_shading(mt1, l1, p, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    color = mt1.color * s1.diffuse * s1.shadow;
                    color += s1.specular;
                    color *= s1.aoc;
                    color += mt1.color * s1.amb *.09;
                    color = pow( color, vec3(1.0/2.2) );
                }
                if (red_cube(p) == scene(p)) {
                    mt1.color = vec3(1.,0.,0.);
                    mt1.reflection_ratio = .2;
                    mt1.shininess = 10.9;
                    s1 = get_shading(mt1, l1, p, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    color = mt1.color * s1.diffuse * s1.shadow;
                    color += s1.specular;
                    color *= s1.aoc;
                    color += mt1.color * s1.amb *.09;
                }
                if (blue_cube(p) == scene(p)) {
                    mt1.color = vec3(0.,0.07,1.);
                    mt1.reflection_ratio = .2;
                    mt1.shininess = 10.9;
                    s1 = get_shading(mt1, l1, p, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    color = mt1.color * s1.diffuse * s1.shadow * 1.2;
                    color += s1.specular;
                    color *= s1.aoc;
                    color += mt1.color * s1.amb *.3;
                }
                if (green_cube(p) == scene(p)) {
                    mt1.color = vec3(0.,1.,0.);
                    mt1.reflection_ratio = .2;
                    mt1.shininess = 10.9;
                    s1 = get_shading(mt1, l1, p, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    color = mt1.color * s1.diffuse * s1.shadow;
                    color += s1.specular;
                    color *= s1.aoc;
                    color += mt1.color * s1.amb *.09;
                }
                if (plane(p) == scene(p)) {
                    float checker = mod(floor(p.x)+floor(p.z)-.5, 2.0);
                    mt1.color = vec3(checker);
                    mt1.reflection_ratio = 0.01;
                    mt1.shininess = 3.;
                    s1 = get_shading(mt1, l1, p, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    color = mt1.color * s1.diffuse * s1.shadow;
                    color += s1.specular;
                    color *= s1.aoc;
                    color += mt1.color * s1.amb *.09;
                }
                if (refracted) {
                    mt1.color = vec3(0.,0.,1.);
                    mt1.reflection_ratio = 10.9;
                    mt1.shininess = 100.9;
                    s1 = get_shading(mt1, l1, p_refr, eye);
                    color += s1.specular*10.;
                }
                if (reflected) {
                    // attenuate reflected color
                    color *= 0.8;
                }
                break;
            }
        }
        step += d;
    }

    return color;
}

void main() {
    vec3 col;
    #ifdef AA
        // do some antialiasing
        float AA_size = 2.0;
        float count = 0.0;
        for (float aaY = 0.0; aaY < AA_size; aaY++) {
            for (float aaX = 0.0; aaX < AA_size; aaX++) {
                col += draw(gl_FragCoord.xy + vec2(aaX, aaY) / AA_size);
                count += 1.0;
            }
        }
        col /= count;
    #else
        col = draw(gl_FragCoord.xy);
    #endif

    // gamma correction
    col = pow( col, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(col, 1.);

}
