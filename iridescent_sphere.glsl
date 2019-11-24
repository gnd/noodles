#version 460
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
uniform sampler3D RGB2HSI;
uniform sampler3D HSI2RGB;
uniform sampler2D backbuffer;
#define MAXSTEPS 128
#define MAXDIST 20.0
#define eps 0.0001

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material m1,m2,m3;
shading s1,s2,s3;

vec3 sky, color; // 'sky color'

/* takes ~40-50ms on geforce 940m @ 1280x720 */

// Sinusoid bumps
float sinbumps(in vec3 p) {
    float frq = 10.;
    return sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31)*sin(time)*.07;
}

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

float iridescent_sphere1(vec3 p) {
    return sinbumps(p) + sphere(p-vec3(1.2,1.,1.), .7);
}

float iridescent_sphere2(vec3 p) {
    return sinbumps(p) + sphere(p-vec3(-1.2,1.,1.), .7);
}

float iridescent_sphere3(vec3 p) {
    return sinbumps(p) + sphere(p-vec3(0.,1.,-1.2), .7);
}

float plane(vec3 p) {
    float height = .5;
    float width = 15.;
    return box(p-vec3(.0,-1.*height,.0), vec3(width,.01,width), .0);
}

float scene(vec3 p) {
    return min(min(min(iridescent_sphere3(p), iridescent_sphere2(p)), iridescent_sphere1(p)), plane(p));
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

// from https://www.shadertoy.com/view/llcXWM
vec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ) {
    return a + b*cos( 6.28318*(c*t+d) );
}

vec3 spectrum(float n) {
    return pal( n, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );
}

shading get_shading(material m, light l, vec3 p, vec3 n, vec3 eye) {
    shading s;
    float step;
    vec3 light_dir = normalize(l.position - p); // direction from p to light

    // Diffuse shading using Lambertian reflectance
    // https://en.wikipedia.org/wiki/Lambertian_reflectance
    s.diffuse = clamp(dot(n,normalize(light_dir)), 0., 1.);

    // Specular reflections using phong
    // https://en.wikipedia.org/wiki/Phong_reflection_model
    vec3 reflection_dir = normalize(reflect(-light_dir, n));
    s.specular = clamp(m.reflection_ratio * pow(dot(reflection_dir, normalize(eye-p)), m.shininess), 0., 1.);

    // Shadows - basically raymarch towards the light
    // If the path from p to each light source is obstructed add shadow
    // soft shadows taken from http://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
    s.shadow = 1.;
    step = 0.01;
    float ph = 1e10;
    for (int i = 0; i < MAXSTEPS/2; i++ ) {
        vec3 shadow_p = p + light_dir * step;
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

    // Ambient occlusion
    // Raymarch along normal a few steps, if hit a surface add shadow
    s.aoc = 1.;
    step = .0;
    for (int k =0; k < 5; k++) {
        step += .04;
        vec3 aoc_p = p + n * step;
        float diff = 1.5*abs(length(aoc_p-p) - scene(aoc_p));
        s.aoc = min(s.aoc, s.aoc-diff);
    }

    // Ambient light
    s.amb = clamp(0.5+0.4*n.y, 0., 1.);

    return s;
}

vec3 draw(in vec2 fragcoord) {
    vec3 eye = vec3(cos(time/4.)*3., 2., sin(time/4.)*3.);
    //eye = vec3(-.2,2., -3.);
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
    sky = vec3(.0);
    color = sky;

    // raymarch a scene
    float step = .0;
    vec3 p;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = ro + rd * step;
        float d = scene(p);
        if (d > MAXDIST) {
            break;
        }
        if (d < eps) {
            vec3 n = normal(p);
            if (iridescent_sphere1(p) == scene(p)) {
                vec3 pp = mod(p*100., 2.) - .5;
                vec3 i = normalize(eye - pp);
                // first try at iridescnce -> color changes according to the angle of view
                vec3 ir = refract(i, n, 1./1.);
                vec3 ir2 = refract(i, -n, 1./2.);
                vec3 ir_color = vec3(cross(ir, n));
                vec3 ir_color2 = vec3(cross(ir2, -n));
                m1.color = ir_color2;
                m1.reflection_ratio = .9;
                m1.shininess = 0.9;
                s1 = get_shading(m1, l1, p, n, eye);
                m2.color = vec3(1.);
                m2.reflection_ratio = 2.9;
                m2.shininess = 1.9;
                s2 = get_shading(m2, l1, p, n, eye);
                m3.color = vec3(1.);
                m3.reflection_ratio = 0.5;
                m3.shininess = 0.1;
                s3 = get_shading(m3, l1, p, n, eye);
                // mix it
                color = ir_color2 * s1.diffuse * s1.shadow * 1.;
                color += ir_color2*1.3 + s2.specular * .1;
                color += s3.specular * 0.1;
                color += ir_color * s1.specular * 2.5;
                color += m1.color * s1.amb * 1.1;
            }
            if (iridescent_sphere2(p) == scene(p)) {
                float nv = dot(n, -rd);
                vec3 col = vec3(0.);
                // https://www.shadertoy.com/view/tdBXzR tweaked
                col += sin(nv * vec3(0.0, 1.0, 0.0) * 10.0 * 1.5) * 0.5 + 0.5;
                col += sin(nv * vec3(1.0, 0.0, 0.0) * 20.0 * 1.5) * 0.5 + 0.5;
                col += sin(nv * vec3(0.0, 0.0, 1.0) * 5.0 * 1.5) * 0.5 + 0.5;
                col = 1.1 - col;
                m1.color = clamp(normalize(col), 0.0, 1.0);
                m1.reflection_ratio = 3.5;
                m1.shininess = 70.;
                s1 = get_shading(m1, l1, p, n, eye);
                color = brightcon(m1.color, .1, 1.5) * s1.diffuse * 1.;
                color += s1.specular;
                //color *= s1.aoc;
                color += m1.color * s1.amb * 1.1;
                // a slight blue glov
                color += vec3(0.,0.,float(i)*(1./float(MAXSTEPS/3)))*.05;
            }
            // the tick was to remove the diffuse element & lower the abient and make the color almost black
            // the colors show up only as part of the specular component
            if (iridescent_sphere3(p) == scene(p)) {
                //https://www.shadertoy.com/view/llcXWM
                vec3 perturb = sin(p * 10.);
                m1.color = spectrum( dot(n + perturb * .05, eye) * 2.);
                m1.reflection_ratio = .5;
                m1.shininess = .1;
                s1 = get_shading(m1, l1, p, n, eye);
                color = 0.01 - m1.color;
                color += s1.specular;
                color += m1.color * s1.amb * 0.5;
            }
            if (plane(p) == scene(p)) {
                float checker = mod(floor(p.x)+floor(p.z)-.5, 2.0);
                m1.color = vec3(checker);
                m1.reflection_ratio = 0.01;
                m1.shininess = 3.;
                s1 = get_shading(m1, l1, p, n, eye);
                color = m1.color * s1.diffuse * s1.shadow/2.;
                color += s1.specular;
                color *= s1.aoc;
                color += m1.color * s1.amb *.2;
            }
            break;
        }
        step += d;
    }

    // Exponential distance fog
    color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * step * step));

    return color;
}

void main() {
    vec3 col;
    col = draw(gl_FragCoord.xy);

    // gamma correction
    col = pow( col, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(col, 1.);
}
