#version 460
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
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

float plane(vec3 p) {
    float height = .5;
    float width = 15.;
    return box(p-vec3(.0,-1.*height,.0), vec3(width,.01,width), .0);
}

float scene(vec3 p) {
    return plane(p);
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
    color = vec3(.0);
    sky = color;

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
