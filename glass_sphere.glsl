#version 460
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
uniform sampler2D backbuffer;
#define MAXSTEPS 128
#define MAXDIST 20.0
#define eps 0.0001

struct march { float step; vec3 ro; vec3 rd; };
struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float intensity; float shadow; float aoc; float amb; };

bool inside1 = false;
bool inside2 = false;
bool reflect1 = false;

light l1;
material m1,m2;
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
    return sphere(p - vec3(0., 2., 0.), .4);
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

float mirror(vec3 p) {
    float height = .45;
    return box(p-vec3(.0,-1.*height,.0), vec3(1.7,.01,1.7), .0);
}

float plane(vec3 p) {
    float height = .5;
    return box(p-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
}

float scene(vec3 p) {
    return min(min( min(glass_sphere_lil(p), glass_sphere(p)), min(mirror(p), plane(p))), min(neon_sphere(p), red_cube(p)));
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
            // attenuate sky color to 80%
            color *= 0.8;
            break;
        }
    }
    return march;
}

float fresnel(vec3 V, vec3 N, float R0) {
    float cosAngle = 1.0-max(dot(V, N), 0.0);
    float result = cosAngle * cosAngle;
    result = result * result;
    result = result * cosAngle;
    result = clamp(result * (1.0 - R0) + R0, 0.0, 1.0);
    return result;
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

    // Determine intensity of light at p
    float light_dist = scene(l.position);
    s.intensity = light_dist/length(light_dir);
    s.intensity = pow(s.intensity, 1.);

    // Shadows - basically raymarch towards the light
    // If the path from p to each light source is obstructed add shadow
    s.shadow = 1.;
    step = 0.0025;
    for (int i = 0; i < MAXSTEPS/3; i++ ) {
        vec3 shadow_p = p + normalize(light_dir) * step;
        float shadow_d = scene(shadow_p);
        if (shadow_d > length(light_dir)) {
            break;
        }
        if ((shadow_d < eps) || (step > MAXDIST*2)) {
            s.shadow = 0.;
            break;
        }
        // check http://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
        s.shadow = min( s.shadow, 32.*shadow_d/step );
        step += shadow_d;
    }
    s.shadow = clamp(s.shadow, 0., 1.);

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

void main() {
    vec3 eye = vec3(cos(time/3.)*3., 2., sin(time/3.)*3.);
    vec3 lookat = vec3(0.,0.,0.);
    vec3 fwd = normalize(lookat-eye);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x ));
    vec3 up = normalize(cross(fwd, right));
    float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 ro = eye;
    vec3 rd = normalize(1.4*fwd + u*right*aspect + v*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.),2.,cos(time/4.));
    color = vec3(1.);

    // raymarch a scene
    float step = .0;
    vec3 p_tmp, pp, ppp, pppp;
    float eta;
    for (int i = 0; i < MAXSTEPS; ++i) {
        vec3 p = ro + rd * step;
        float d = scene(p);
        if (d > MAXDIST) {
            break;
        }
        if (d < eps) {
            if ((glass_sphere(p) == scene(p)) || (glass_sphere_lil(p) == scene(p)) || mirror(p) == scene(p)) {
                if (mirror(p) == scene(p)) {
                    // Compute reflection
                    pp = p;
                    reflect1 = true;
                    vec3 ld = normalize(eye - p);
                    march march = reflection(p, ld);
                    step = march.step;
                    rd = march.rd;
                    ro = march.ro;
                }
                if (glass_sphere_lil(p) == scene(p)) {
                    // Glass becomes a 'portal' through which refracted rays keep on raymarching
                    ppp = p;
                    inside2 = true;
                    // eta here is a ratio of source material refraction index with the dest material refraction index
                    // eg. air: 1, water 1.5 => eta = 1./1.5
                    // eta = 1.0/1.5; // air / glass
                    eta = 1.0/2.22; // air / diamond
                    march march = refraction(p, rd, eta);
                    step = march.step;
                    rd = march.rd;
                    ro = march.ro;
                }
                if (glass_sphere(p) == scene(p)) {
                    // Glass becomes a 'portal' through which refracted rays keep on raymarching
                    pp = p;
                    inside1 = true;
                    eta = 1.0/1.33; // air / water
                    march march = refraction(p, rd, eta);
                    step = march.step;
                    rd = march.rd;
                    ro = march.ro;
                }
            } else {
                if (neon_sphere(p) == scene(p)) {
                    m1.color = vec3(1.,1.,.0);
                    m1.reflection_ratio = .9;
                    m1.shininess = 15.9;
                    s1 = get_shading(m1, l1, p, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    color = m1.color * s1.diffuse * s1.intensity * s1.shadow;
                    color += s1.specular;
                    color += m1.color * s1.aoc * s1.amb * 0.5;
                }
                if (red_cube(p) == scene(p)) {
                    m1.color = vec3(1.,0.,0.);
                    m1.reflection_ratio = .1;
                    m1.shininess = 0.9;
                    s1 = get_shading(m1, l1, p, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    color = m1.color * s1.diffuse * s1.intensity * s1.shadow;
                    color += s1.specular;
                    color += m1.color * s1.aoc * s1.amb * 0.6;
                }
                if (plane(p) == scene(p)) {
                    m1.color = vec3(.1,.1,1.);
                    m1.reflection_ratio = .1;
                    m1.shininess = 0.1;
                    s1 = get_shading(m1, l1, p, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    color = m1.color * s1.diffuse * s1.intensity * s1.shadow * 5.;
                    //color += s1.specular;
                    color *= s1.aoc;
                    color += s1.amb *.05;
                }
                if (reflect1) {
                    // attenuate reflected color
                    color *= 0.8;
                }
                if (inside1) {
                    // now add specular & reflection to the glass sphere
                    m1.color = vec3(0.,0.,1.);
                    m1.reflection_ratio = 10.9;
                    m1.shininess = 100.9;
                    s1 = get_shading(m1, l1, pp, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    float fr = fresnel(normalize(eye-pp), normal(pp), 0.5);
                    vec3 reflection = fr * texture(backbuffer, reflect(-1.*normalize(eye-pp).xy, normal(pp).xy)).xyz;
                    color += reflection*0.05;
                    color += s1.specular*10.;
                }
                if (inside2) {
                    // now add specular & reflection to the glass sphere
                    m1.color = vec3(0.,0.,1.);
                    m1.reflection_ratio = 10.9;
                    m1.shininess = 100.9;
                    s1 = get_shading(m1, l1, ppp, eye);
                    // check https://www.shadertoy.com/view/lsKcDD
                    float fr = fresnel(normalize(eye-pp), normal(pp), 0.5);
                    vec3 reflection = fr * texture(backbuffer, reflect(-1.*normalize(eye-pp).xy, normal(pp).xy)).xyz;
                    //color += m1.color * 1.;
                    color += reflection*0.05;
                    color += s1.specular*10.;
                }
                break;
            }
        }
        step += d;
    }

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    PixelColor = vec4(color, 0.);
}
