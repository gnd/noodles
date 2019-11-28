#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
uniform sampler2D sony;
#define MAXSTEPS 128
#define MAXDIST 20.0
#define eps 0.0001
#define FLY

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material m1;
shading s1;

vec3 sky, color; // 'sky color'
float speed_z = 4.;
float speed_x = speed_z * .5;
vec2 xz = vec2(cos(time/3.)*speed_x,time*speed_z);
vec3 light_position = vec3(xz.x + sin(time/4.)*0.,2.,xz.y + 8. + cos(time/4.)*0.);

// Random number generator
// implementation found at: lumina.sourceforge.net/Tutorials/Noise.html
float rnd(vec2 co){
        float k = sin(dot(co.xy ,vec2(12.9898,78.233)));
        return fract(k + k);
}

//https://www.shadertoy.com/view/XsXfRH
float hash(vec3 p) {
    p  = 50.0*fract( p*0.3183099 + vec3(0.71,0.113,0.419));
    return -1.0+2.0*fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
}

// http://www.iquilezles.org/www/articles/morenoise/morenoise.htm
float noise2f( in vec2 p ) {
        vec2 ip = vec2(floor(p));
        vec2 u = fract(p);
        u = u*u*u*((6.0*u-15.0)*u+10.0);
        float res = mix( mix(rnd(ip), rnd(ip+vec2(1.0,0.0)), u.x),
                         mix(rnd(ip+vec2(0.0,1.0)), rnd(ip+vec2(1.0,1.0)), u.x),
                        u.y);
        return res - 0.25;
}

float sphere(vec3 p, float r ) {
    return length(p) - r;
}

float maxcomp(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

float box(vec3 p, vec3 b, float r) {
    vec3 d = abs(p) - b;
    return min(maxcomp(d),0.0) - r + length(max(d,0.0));
}

float hm1(vec2 p) {
    return sin(p.r)*cos(p.g);
}

vec3 mask(vec3 p) {
    return (((p.x > 9) || (p.x < -11.)) || ((p.z > 4.) || (p.z < -6.))) ? vec3(0.) : vec3(1.);
}

vec3 texcol(vec3 p) {
    vec3 image = texture2D(sony, (vec2(9.-p.x,4.-p.z))*vec2(.05,.1)).xyz;
    return image;
}

float tex(vec3 p) {
    return dot(texcol(p), vec3(0.299, 0.587, 0.114))*2.;
}

float hm2(vec3 p) {
    return tex(p);
}

float plane(vec3 p) {
	return p.y + noise2f(p.xz/10.)*10. + noise2f(p.xz) + noise2f(p.xz*5.)*.1;
}

float sky_plane(vec3 p) {
	return -p.y + 10. - noise2f(p.xz/10.)*10. + noise2f(p.xz) + noise2f(p.xz*5.)*.1;
}

float lightsphere(vec3 p) {
    return sphere(p, 2.);
}

float scene(vec3 p) {
    return min(plane(p), sky_plane(p));
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
    return b*cos( 6.28318*(c*t+d) ) + a;
}

vec3 spectrum(float n) {
    return pal( n, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );
}

shading get_simple(material m, vec3 n, vec3 ld, vec3 ed) {
    shading s;
    s.diffuse = clamp(dot(n,ld), 0., 1.);
    s.specular = clamp(m.reflection_ratio * pow(dot(normalize(reflect(-ld, n)), ed), m.shininess), 0., 1.);
    s.amb = clamp(0.4*n.y+0.5, 0., 1.);
    return s;
}

shading get_shading(material m, light l, vec3 p, vec3 n, vec3 ld, vec3 ed) {
    shading s;
    float step;
    s.diffuse = clamp(dot(n,ld), 0., 1.);
    s.specular = clamp(m.reflection_ratio * pow(dot(normalize(reflect(-ld, n)), ed), m.shininess), 0., 1.);

    s.shadow = 1.;
    step = 0.01;
    float ph = 1e10;
    /*
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
    }*/
    s.shadow = clamp(s.shadow, 0., 1.);
    s.aoc = 1.;
    s.amb = clamp(0.4*n.y+0.5, 0., 1.);

    return s;
}

void main() {
    speed_z = 4.;
    speed_x = speed_z * .5;
    #ifdef FLY
        xz = vec2(cos(time/3.)*speed_x,time*speed_z);
        vec2 xz_ahead = vec2(cos(time+2./3.)*speed_x,(time+2.)*speed_z);
        vec3 ro = vec3(xz.x,2.-noise2f(xz.xy/10.)*10.,xz.y);
        vec3 lookat = vec3(xz_ahead.x,2.-noise2f(xz_ahead.xy/10.)*10.,xz_ahead.y);
        l1.position = light_position;
    #else
        vec3 ro = vec3(cos(time/4.)*5., 5., sin(time/4.)*5.); // ray origin, here also known as 'eye'
        vec3 lookat = vec3(0.,0.,0.);
        l1.position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    #endif
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    color = vec3(0.);
    sky = color;

    // raymarch a scene
    float step = .0;
    vec3 p;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = ro + (rd * step);
        float d = scene(p);
        if (d > MAXDIST) {
            break;
        }
        if (scene(p) < eps) {
            // precompute stuff
            vec3 n = normal(p);
            vec3 ld = normalize(l1.position-p);
            vec3 ed = normalize(ro-p);
            if (plane(p) < eps) {
                /*
                vec3 perturb = sin(p * .001);
                m1.color = spectrum( dot(perturb * .05 + n, ro) * 2.);
                color = 0.02 - m1.color;
                color += clamp(.5 * pow(dot(normalize(reflect(-ld, n)), ed), .1), 0., 1.) * 2.;
                color += m1.color * clamp(0.4*n.y+0.5, 0., 1.) * .5;
                */
                /*
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
                */
            }
            if (sky_plane(p) < eps) {
                /*
                vec3 perturb = sin(p * .001);
                m1.color = spectrum( dot(perturb * .05 + n, ro) * 2.);
                color = 0.07 - m1.color;
                color += clamp(.5 * pow(dot(normalize(reflect(-ld, n)), ed), .1), 0., 1.) * 2.;
                color += m1.color * clamp(0.4*n.y+0.5, 0., 1.) * 1.5;
                */
                /*
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
                */
            }
            break;
        }
        color += vec3(0.,0.,float(i)*(1./float(MAXSTEPS)))*0.2; //glow
        step += d;
    }

    // Exponential distance fog
    //color = mix(color, 0.8 * sky, .7 - exp2(-0.01 * step * step));

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(color, 1.);
}
