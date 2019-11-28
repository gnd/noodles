#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define MAXSTEPS 64
#define MAXDIST 20.0
#define eps 0.01
#define SINBUMPS
//#define REPEATPLANE
//#define REPEAT
//#define FLY

struct march { float step; vec3 ro; vec3 rd; };
struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material m1,m2;
shading s1,s2;

vec3 color = vec3(0.); //sky color
bool refracted = false;

// smooth min
float smin( float a, float b, float k ) {
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*h*k*(1.0/6.0);
}

// Sinusoid bumps
float sinbumps(in vec3 p){
    float frq = 1.7;
    return sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31)*sin(time)*1.7;
}

float simple_bumps(in vec3 p) {
    float frq = 1.7;
    return (sin(p.x*frq)+cos(p.z*frq))/1.;
}

float sphere( vec3 p, float r ) {
    return length(p) - r;
}

float vmax(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

// Cheap Box: distance to corners is overestimated
// from http://mercury.sexy/hg_sdf/
float cbox(vec3 p, vec3 b) { //cheap box
	return vmax(abs(p) - b);
}

float scene_plane(vec3 p) {
    float height = .5;
    #ifdef REPEATPLANE
        p.x = mod(p.x,10.) - 5.;
        p.z = mod(p.z,10.) - 5.;
    #endif
    #ifdef SINBUMPS
        float plane_field = cbox(vec3(p.x,p.y+sinbumps(p),p.z)-vec3(.0,-1.*height,.0), vec3(20.,.01,20.));
        //float plane_field = box(vec3(p.x,p.y+simple_bumps(p),p.z)-vec3(.0,-1.*height,.0), vec3(20.,.01,20.), .0);
    #else
        float plane_field = cbox(p-vec3(.0,-1.*height,.0), vec3(20.,.01,20.));
    #endif
    return plane_field;
}

float scene_object(vec3 p) {
    #ifdef REPEAT
        p.x = mod(p.x,1.) - .5;
        p.z = mod(p.z,2.) - 1.;
    #endif
    return smin(cbox(p, vec3(.1,.1,2.)), smin(cbox(p, vec3(1.,0.1,0.1)),sphere(p,.333), .4), .4);
}

float glass_sphere(vec3 p) {
    return sphere(p-vec3(0.,1.5,.0), .8);
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

// from https://www.shadertoy.com/view/llcXWM
vec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ) {
    return a + b*cos( 6.28318*(c*t+d) );
}

vec3 spectrum(float n) {
    return pal( n, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );
}

march refraction(vec3 p, vec3 n, vec3 rd, float eta) {
    // Glass becomes a 'portal' through which refracted rays keep on raymarching
    march march;
    march.ro = p;
    march.rd = normalize(refract(rd, n, eta));
    march.step = 0.0;
    for (int j=0; j < MAXSTEPS; j++) {
        p = march.ro + march.rd * march.step;
        float d = scene(p);
        march.step += max(abs(d),eps);
        if (d > eps) {
            // attenuate due to impurities, etc
            color *= 1.;
            // second refraction
            march.ro = p;
            march.rd = normalize(refract(march.rd, -normal(p), 1./eta));
            break;
        }
    }
    return march;
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
    #ifdef FLY
        vec3 eye = vec3(cos(time/3.), sin(time/3.)*3., time);
        vec3 lookat  = vec3(0.,0.,7.+time);
        l1.position = vec3(cos(time)*10.,2.,sin(time)*10.+time);
    #else
        vec3 eye = vec3(cos(time/3.)*3., sin(time/3.)*0. + 2., sin(time/3.)*3.);
        vec3 lookat = vec3(0.,0.,0.);
        l1.position = vec3(cos(time)*5.,4.,sin(time)*5.); // light position
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
    color = vec3(0.,0.,0.); // 'sky color

    // raymarch a scene
    float step = .0;
    vec3 p, p_refr;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = rd * step + ro;
        float d = scene(p);
        // if ray doesnt hit any surface, kill it after being longer than MAXDIST
        if (d > MAXDIST) {
            break;
        }
        // if ray very close to a surface, find out how to color that surface AKA what color is the pixel
        if (d < eps) {
            // precompute stuff
            vec3 n = normal(p);
            vec3 ld = normalize(l1.position-p);
            vec3 ed = normalize(ro-p);
            if (glass_sphere(p) < eps) {
                // Compute refraction
                p_refr = p;
                refracted = true;
                march march = refraction(p, n, rd, 1./2.22);
                step = march.step;
                rd = march.rd;
                ro = march.ro;
            } else {
                // shading the plane
                if (scene_plane(p) < eps) {
                    m2.color = vec3(0.,0.,1.)*10.;
                    m2.reflection_ratio = 0.9;
                    m2.shininess = 100.9;
                    s2 = get_shading(m2, l1, p, n, ld, ed);
                    color += m2.color * s2.diffuse * s2.shadow * 2.;
                    color += s2.specular;
                    color += m2.color * s2.amb * 0.05;
                }
                // shading the object
                if (scene_object(p) < eps) {
                    vec3 perturb = sin(p * 10.);
                    vec3 icolor = spectrum( dot(normal(p) + perturb * .05, eye) * 2.);
                    icolor = .7 -icolor;
                    m1.color = vec3(.0,0.,0.);
                    s1 = get_shading(m1, l1, p, n, ld, ed);
                    color += icolor * 0.04;
                    color += color * 0.1 + m1.color * clamp(dot(n,ld), 0., 1.);
                    color += clamp(1.5 * pow(dot(normalize(reflect(-ld, n)), ed), 10.9), 0., 1.) * .3;
                    color += m1.color * clamp(0.4*n.y+0.5, 0., 1.) * 0.01;
                }
                if (refracted) {
                    color += clamp(10.9 * pow(dot(normalize(reflect(-ld, n)), ed), 100.9), 0., 1.) * 10.;
                }
                break;
            }
        }
        step += d;
        // red haze makes it lil bit interesting
        color += vec3(float(i)*(1./float(MAXSTEPS/3)), 0.,0.)*.01;
    }

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    PixelColor = vec4(color, 0.);
}
