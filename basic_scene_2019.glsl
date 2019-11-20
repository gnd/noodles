//taken from https://9bitscience.blogspot.com/2013/07/raymarching-distance-fields_14.html
uniform vec2 resolution;
uniform float time;
#define MAXSTEPS 128
#define MAXDIST 20.0
// should render sinbumps ?
//#define SINBUMPS
//#define REPEAT

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

float scene(vec3 p) {
    float height = .5;
    #ifdef SINBUMPS
        float plane_field = box(vec3(p.x,p.y+sinbumps(p),p.z)-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
    #else
        float plane_field = box(p-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
    #endif
    #ifdef REPEAT
        // if doing domain repeat it is vital to avoid problems with zero-centered objects
        // push the domain by half of modulo to the side after mod
        p.x = mod(p.x,1.) - .5;
        p.z = mod(p.z,2.) - 1.;
    #endif
    return min(plane_field ,min(box(p, vec3(.1,.1,2.), .0), min(box(p, vec3(1.,0.1,0.1),0.),sphere(p,.333))));
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

void main() {
    vec3 eye = vec3(cos(time/3.), sin(time/3.)+1., -3.); // camera position
    vec3 up = vec3(0, 1, 0);
    vec3 right = vec3(1, 0, 0);
    vec3 fwd = normalize(vec3(0.)-eye);
    vec3 light = vec3(cos(time),2.,sin(time)); // light position (point light)
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

            // lighting follows
            vec3 light_col = vec3(1.,.8,1.)*1.; // light color
            vec3 material_col = vec3(1.,1.,1.); // material color
            vec3 light_dir = light-p; // direction from p to light
            vec3 normal = normal(p, 0.001); // surface normal at point p
            float reflection_ratio = 0.8; // for specular
            float shininess = 10.9; // for specular

            // Diffuse shading using Lambertian reflectance
            // https://en.wikipedia.org/wiki/Lambertian_reflectance
            float diffuse = clamp(dot(normal,normalize(light_dir)), 0., .4);

            // Determine intensity of light at p
            float light_dist = scene(light);
            float light_intensity = light_dist/length(light_dir);

            // Specular reflections using phong
            // https://en.wikipedia.org/wiki/Phong_reflection_model
            vec3 reflection_dir = 2.*dot(normalize(light_dir), normal)*normal-normalize(light_dir);
            float specular = clamp(reflection_ratio * pow(dot(normalize(reflection_dir), normalize(eye-p)), shininess), 0., 1.);

            // Shadows - basically raymarch towards the light
            // If the path from p to each light source is obstructed add shadow
            float shadow_col = 1.;
            float shadow_t = 0.1;
            float soft_shadow = 1.;
            // we dont have to use MAXSTEPS here, 20 seems to be enough
            for (int j = 0; j < 20; j++) {
                vec3 shadow_p = p + normalize(light_dir) * shadow_t;
                float shadow_d = scene(shadow_p);
                if (shadow_d > length(light_dir)) {
                    break;
                }
                if (shadow_d < 0.01) {
                    shadow_col = 0.;
                    break;
                }
                // check http://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
                soft_shadow = min(soft_shadow, shadow_d);
                shadow_t += shadow_d;
            }

            // Ambient occlusion
            // Raymarch along normal a few steps, ig hit a surface add shadow
            float aoc = 1.;
            float aoc_t = .0;
            for (int k =0; k < 5; k++) {
                aoc_t += .04;
                vec3 aoc_p = p + normal * aoc_t;
                float diff = .9*abs(length(aoc_p-p) - scene(aoc_p));
                aoc = min(aoc, aoc-diff);
            }

            // Put everything together
            color += material_col * light_col * pow(light_intensity, 3.) * diffuse * aoc * shadow_col * soft_shadow*10.;
            color += specular;
            break;
        }
        step += d;
        // red haze makes it lil bit interesting
        color = vec3(float(i)*(1./32.), 0.,0.);
    }

    gl_FragColor = vec4(color, 0.);
}
