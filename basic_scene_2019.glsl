//taken from https://9bitscience.blogspot.com/2013/07/raymarching-distance-fields_14.html
uniform vec2 resolution;
uniform float time;
#define MAXSTEPS 128
#define MAXDIST 20.0
// should render sinbumps ?
//#define SINBUMPS
//#define REPEAT

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float intensity; float shadow; float soft_shadow; float aoc; };

light l1;
material m1,m2;
shading s1,s2;

// Sinusoid bumps
float sinbumps(in vec3 p){
    float frq = 1.7;
    return sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31)*0.4;
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

float scene_plane(vec3 p) {
    float height = .5;
    #ifdef SINBUMPS
        float plane_field = box(vec3(p.x,p.y+sinbumps(p),p.z)-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
    #else
        float plane_field = box(p-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
    #endif
    return plane_field;
}

float scene_object(vec3 p) {
    return min(box(p, vec3(.1,.1,2.), .0), min(box(p, vec3(1.,0.1,0.1),0.),sphere(p,.333)));
}

float scene(vec3 p) {
    #ifdef REPEAT
        // if doing domain repeat it is vital to avoid problems with zero-centered objects
        // push the domain by half of modulo to the side after mod
        p.x = mod(p.x,1.) - .5;
        p.z = mod(p.z,2.) - 1.;
    #endif
    return min(scene_plane(p), scene_object(p));
}

// aproximate surface normal by actually doing a difference between:
//  A' = A + dx; A'' = A - dx;
// dist(A') - dist(A) - (dist(A'') - dist(A))
// = dist(A') - dist(A'')
// the resulting vector is pointing/growing in the direction of the normal
vec3 normal(vec3 p, float dx) {
    vec3 n = vec3(0.);
    vec3 dxx = vec3(1.,.0,.0)*dx;
    vec3 dxy = vec3(0.,1.,0.)*dx;
    vec3 dxz = vec3(0.,0.,1.)*dx;
    n.x = (scene(p + dxx) - scene(p - dxx));
    n.y = (scene(p + dxy) - scene(p - dxy));
    n.z = (scene(p + dxz) - scene(p - dxz));
    return normalize(n);
}

shading get_shading(material m, light l, vec3 p, vec3 eye) {
    shading s;
    float step;
    vec3 light_dir = l.position-p; // direction from p to light
    vec3 normal = normal(p, 0.001); // surface normal at point p

    // Diffuse shading using Lambertian reflectance
    // https://en.wikipedia.org/wiki/Lambertian_reflectance
    s.diffuse = clamp(dot(normal,normalize(light_dir)), 0., .4);

    // Specular reflections using phong
    // https://en.wikipedia.org/wiki/Phong_reflection_model
    vec3 reflection_dir = 2.*dot(normalize(light_dir), normal)*normal-normalize(light_dir);
    s.specular = clamp(m.reflection_ratio * pow(dot(normalize(reflection_dir), normalize(eye-p)), m.shininess), 0., 1.);

    // Determine intensity of light at p
    float light_dist = scene(l.position);
    s.intensity = light_dist/length(light_dir);
    s.intensity = pow(s.intensity, 3.);

    // Shadows - basically raymarch towards the light
    // If the path from p to each light source is obstructed add shadow
    s.shadow = 1.;
    s.soft_shadow = 1.;
    step = 0.1;
    // we dont have to use MAXSTEPS here, 20 seems to be enough
    for (int j = 0; j < 20; j++) {
        vec3 shadow_p = p + normalize(light_dir) * step;
        float shadow_d = scene(shadow_p);
        if (shadow_d > length(light_dir)) {
            break;
        }
        if (shadow_d < 0.01) {
            s.shadow = 0.;
            break;
        }
        // check http://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
        s.soft_shadow = min(s.soft_shadow, shadow_d);
        step += shadow_d;
    }

    // Ambient occlusion
    // Raymarch along normal a few steps, ig hit a surface add shadow
    s.aoc = 1.;
    step = .0;
    for (int k =0; k < 5; k++) {
        step += .04;
        vec3 aoc_p = p + normal * step;
        float diff = .9*abs(length(aoc_p-p) - scene(aoc_p));
        s.aoc = min(s.aoc, s.aoc-diff);
    }

    return s;
}

void main() {
    vec3 eye = vec3(cos(time/3.), sin(time/3.)+1., -3.); // camera position
    vec3 up = vec3(0, 1, 0);
    vec3 right = vec3(1, 0, 0);
    vec3 fwd = normalize(vec3(0.)-eye);
    l1.position = vec3(cos(time),2.,sin(time)); // light position
    l1.color = vec3(1.,1.,1.); // light color
    float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;
    vec3 rd = normalize(fwd + right*u + up*v/1.77); // 'ray distance'
    vec3 ro = eye; // 'ray origin'
    vec3 color = vec3(0.0); // 'sky color'

    // raymarch a scene
    float step = .0;
    for (int i = 0; i < MAXSTEPS; ++i) {
        vec3 p = ro + rd * step;
        float d = scene(p);
        // if ray doesnt hit any surface, kill it after being longer than MAXDIST
        if (d > MAXDIST) {
            break;
        }
        // if ray very close to a surface, find out how to color that surface AKA what color is the pixel
        if (d < 0.001) {
            // red haze - basically just the number of iterations to reach the surface
            color = vec3(float(i)*(1./float(MAXDIST)), 0.,0.)*.7;

            // shading
            m1.color = vec3(1.,1.,1.);
            m1.reflection_ratio = 0.8;
            m1.shininess = 10.9;
            s1 = get_shading(m1, l1, p, eye);
            color += m1.color * l1.color * s1.diffuse * s1.intensity * s1.shadow * s1.soft_shadow*10. * s1.aoc;
            color += s1.specular;
            break;
        }
        step += d;
        // red haze makes it lil bit interesting
        color = vec3(float(i)*(1./32.), 0.,0.);
    }

    gl_FragColor = vec4(color, 0.);
}
