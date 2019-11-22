#version 460
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
uniform sampler2D backbuffer;
#define MAXSTEPS 128
#define MAXDIST 20.0
#define eps 0.0001
#define SINBUMPS
//#define REPEAT
//#define REPEATPLANE
//#define FLY

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float intensity; float shadow; float aoc; float amb; };

light l1;
material m1,m2;
shading s1,s2;

// smooth min
float smin( float a, float b, float k ) {
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*h*k*(1.0/6.0);
}

// Sinusoid bumps
float sinbumps(in vec3 p){
    float frq = 1.7;
    return sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31)*sin(time)*2.;
}

float simple_bumps(in vec3 p) {
    float frq = 1.7;
    return (sin(p.x*frq)+cos(p.z*frq))/1.;
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
    #ifdef REPEATPLANE
        p.x = mod(p.x,10.) - 5.;
        p.z = mod(p.z,10.) - 5.;
    #endif
    #ifdef SINBUMPS
        float plane_field = box(vec3(p.x,p.y+sinbumps(p),p.z)-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
        //float plane_field = box(vec3(p.x,p.y+simple_bumps(p),p.z)-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
    #else
        float plane_field = box(p-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
    #endif
    return plane_field;
}

float scene_object(vec3 p) {
    #ifdef REPEAT
        p.x = mod(p.x,1.) - .5;
        p.z = mod(p.z,2.) - 1.;
    #endif
    return smin(box(p, vec3(.1,.1,2.), .0), smin(box(p, vec3(1.,0.1,0.1),0.),sphere(p,.333), .4), .4);
}

float glass_sphere(vec3 p) {
    return sphere(p-vec3(0.,1.,.0), .5);
}

float scene(vec3 p) {
    return min(min(scene_plane(p), scene_object(p)), glass_sphere(p));
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p) {
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*scene( p + k.xyy*eps ) +
                      k.yyx*scene( p + k.yyx*eps ) +
                      k.yxy*scene( p + k.yxy*eps ) +
                      k.xxx*scene( p + k.xxx*eps ) );
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
    vec3 light_dir = l.position-p; // direction from p to light
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
    #ifdef FLY
        vec3 eye = vec3(cos(time/3.), sin(time/3.)*3., time);
        vec3 lookat  = vec3(0.,0.,7.+time);
        l1.position = vec3(cos(time)*10.,2.,sin(time)*10.+time);
    #else
        vec3 eye = vec3(cos(time/3.)*3., sin(time/3.)*0. + 2., sin(time/3.)*3.);
        vec3 lookat = vec3(0.,0.,0.);
        l1.position = vec3(cos(time)*3.,2.,sin(time)*3.); // light position
    #endif

    vec3 fwd = normalize(lookat-eye);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x ));
    vec3 up = normalize(cross(fwd, right));
    float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 ro = eye;
    vec3 rd = normalize(1.1*fwd + u*right*aspect + v*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    //l1.position = vec3(cos(time)*3.,2.,sin(time)*3.); // light position
    vec3 color = vec3(0.,0.,0.); // 'sky color'

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
        if (d < eps) {
            // shading the plane
            if (scene_plane(p) == scene(p)) {
                m2.color = vec3(0.,0.,1.)*10.;
                m2.reflection_ratio = 0.9;
                m2.shininess = 100.9;
                s2 = get_shading(m2, l1, p, eye);
                // check https://www.shadertoy.com/view/lsKcDD
                color += m2.color * s2.diffuse * s2.intensity * s2.shadow * 2.;
                color += s2.specular;
                color += m2.color * s2.aoc * s2.amb * 0.2;
            }
            // shading the object
            if (scene_object(p) == scene(p)) {
                m1.color = vec3(.0,0.,0.);
                m1.reflection_ratio = .9;
                m1.shininess = 10.9;
                s1 = get_shading(m1, l1, p, eye);
                // check https://www.shadertoy.com/view/lsKcDD
                color = color * 0.9 + m1.color * s1.diffuse;// * s1.intensity *s1.shadow * .1;
                color += s1.specular;
                color += m1.color * s1.aoc * s1.amb * 0.2;
            }
            // shading the glass sphere
            if (glass_sphere(p) == scene(p)) {
                m1.color = vec3(0.0);
                m1.reflection_ratio = 10.9;
                m1.shininess = 100.9;
                s1 = get_shading(m1, l1, p, eye);
                // check https://www.shadertoy.com/view/lsKcDD
                float fr = fresnel(normalize(eye-p), normal(p), 0.5);
                vec3 reflection = fr * texture(backbuffer, reflect(-1.*normalize(eye-p), normal(p)).xy).xyz;
                vec3 refraction = (1.0-fr) * texture(backbuffer, refract(-1.*normalize(eye-p), normal(p), 0.8).xy).xyz;
                color *= m1.color;
    	        color += (reflection + refraction) * .5;
                //color += m1.color * s1.diffuse * s1.intensity * 0.01;// *s1.shadow * 1.;
                color += s1.specular * 2;
                //color += m1.color * s1.aoc * s1.amb * 0.2;
            }
            break;
        }
        step += d;
        // red haze makes it lil bit interesting
        color += vec3(float(i)*(1./float(MAXSTEPS/3)), 0.,0.)*.01;
    }

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    PixelColor = vec4(color, 0.);
}
