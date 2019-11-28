#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define MAXSTEPS 128
#define MAXDIST 8
#define eps 0.001

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material m1,m2,m3;
shading s1,s2,s3;

vec3 sky, color; // 'sky color'

/* un-optimized takes ~40ms on geforce 940m @ 1280x720
   - after optimization managed to get down to 14-20ms
   - further optimization got it down to 10-11ms */

// Sinusoid bumps
float sinbumps(in vec3 p, in float frq) {
    return sin(p.x*frq) * sin(p.y*frq+time) * sin(p.z*frq) * .11;
}

float vmax(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

// Cheap Box: distance to corners is overestimated
// from http://mercury.sexy/hg_sdf/
float cbox(vec3 p, vec3 b) { //cheap box
	return vmax(abs(p) - b);
}

float iridescent_sphere1(vec3 p, float bumps) {
    return length(p-vec3(1.2,1.,1.)) + bumps -.7;
}

float iridescent_sphere2(vec3 p, float bumps) {
    return length(p-vec3(-1.2,1.,1.)) + bumps -.7 ;
}

float iridescent_sphere3(vec3 p, float bumps) {
    return length(p-vec3(0.,1.,-1.2)) + bumps -.7;
}

float plane(vec3 p) {
    return cbox(p+vec3(0.,.5,.0), vec3(15.,.01,15.));
}

float scene(vec3 p) {
    float bumps = sinbumps(p, 10.);
    return min(min(min(iridescent_sphere3(p, bumps), iridescent_sphere2(p, bumps)), iridescent_sphere1(p, bumps)), plane(p));
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p) {
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*scene( k.xyy*eps + p ) +
                      k.yyx*scene( k.yyx*eps + p ) +
                      k.yxy*scene( k.yxy*eps + p ) +
                      k.xxx*scene( k.xxx*eps + p ) );
}

/* same speed as standard notation above ^^
vec3 normal(vec3 p) {
    // inspired by klems - a way to prevent the compiler from inlining map() 4 times
    vec3 n = vec3(0.0);
    for (int i = 0; i < 4; i++) {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*scene(p+0.0005*e);
    }
    return normalize(n);
}*/

// from https://alaingalvan.tumblr.com/post/79864187609/glsl-color-correction-shaders
vec3 brightcon(vec3 c, float brightness, float contrast) {
    return (c - 0.5) * contrast + 0.5 + brightness;
}

// from https://www.shadertoy.com/view/llcXWM
vec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ) {
    return b*cos( 6.28318*(c*t+d) ) + a;
}

vec3 spectrum(float n) {
    return pal( n, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );
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
    vec3 ro = vec3(cos(time/4.)*3., 2., sin(time/4.)*3.);
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
    sky = vec3(.0);
    color = sky;

    // raymarch a scene
    float step = .0;
    vec3 p;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = rd * step + ro;
        float d = scene(p);
        if (d > MAXDIST) {
            break;
        }
        if (d < eps) {
            // precompute stuff
            vec3 n = normal(p);
            vec3 ld = normalize(l1.position-p);
            vec3 ed = normalize(ro-p);
            float bumps = sinbumps(p, 10.);
            if (iridescent_sphere1(p, bumps) < eps) {
                vec3 pp = mod(p*100., 2.) - .5;
                vec3 i = normalize(ro - pp);
                vec3 ir = refract(i, n, 1./1.);
                vec3 ir2 = refract(i, -n, 1./2.);
                vec3 ir_color2 = vec3(cross(ir2, -n));
                // mix it
                color = ir_color2 * clamp(dot(n,ld), 0., 1.);
                color += ir_color2*1.3 + clamp(2.9 * pow(dot(normalize(reflect(-ld, n)), ed), 1.9), 0., 1.) * .1;
                color += clamp(.5 * pow(dot(normalize(reflect(-ld, n)), ed), .1), 0., 1.) * 0.1;
                color += vec3(cross(ir, n)) * clamp(.9 * pow(dot(normalize(reflect(-ld, n)), ed), .9), 0., 1.) * 2.5;
                color += ir_color2 * clamp(0.4*n.y+0.5, 0., 1.) * 1.1;
            }
            if (iridescent_sphere2(p, bumps) < eps) {
                float nv = dot(n, -rd);
                vec3 col = vec3(0.);
                col += sin(nv * vec3(0.0, 1.0, 0.0) * 10.0 * 1.5) * 0.5 + 0.5;
                col += sin(nv * vec3(1.0, 0.0, 0.0) * 20.0 * 1.5) * 0.5 + 0.5;
                col += sin(nv * vec3(0.0, 0.0, 1.0) * 5.0 * 1.5) * 0.5 + 0.5;
                col = 1.1 - col;
                m1.color = clamp(normalize(col), 0.0, 1.0);
                color = brightcon(m1.color, .1, 1.5) * clamp(dot(n,ld), 0., 1.) * 1.;
                color += clamp(3.5 * pow(dot(normalize(reflect(-ld, n)), ed), 70.), 0., 1.);;
                color += m1.color * clamp(0.4*n.y+0.5, 0., 1.) * 1.1;
                color += vec3(0.,0.,float(i)*(1./float(MAXSTEPS/3)))*.05; //glow
            }
            if (iridescent_sphere3(p, bumps) < eps) {
                vec3 perturb = sin(p * 10.);
                m1.color = spectrum( dot(perturb * .05 + n, ro) * 2.);
                //s1 = get_shading(m1, l1, p, n, ld, ed);
                color = 0.01 - m1.color;
                color += clamp(.5 * pow(dot(normalize(reflect(-ld, n)), ed), .1), 0., 1.);
                color += m1.color * clamp(0.4*n.y+0.5, 0., 1.) * 0.5;
            }
            if (plane(p) < eps) {
                m1.color = vec3(mod(floor(p.x)+floor(p.z)-.5, 2.0));
                m1.reflection_ratio = 0.01;
                m1.shininess = 3.;
                s1 = get_shading(m1, l1, p, n, ld, ed);
                color = m1.color * s1.diffuse * s1.shadow/2.;
                color += s1.specular;
                color += m1.color * s1.amb *.2;
            }
            break;
        }
        step += d;
    }

    // Exponential distance fog
    color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * step * step));

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(color, 1.);
}
