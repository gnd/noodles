#version 330
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
material m1;
shading s1;

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
    return p.y;
}

#define R(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)

// c.xy denote position:
// x is x, y is y, z is infinite
// c.z actually encodes cylinder diameter
float cyl(vec3 p, vec3 c) {
  return length(p.xy-c.xy)-c.z;
}

float cylinder(vec3 p) {
    float diameter = 4.;
    float c = cyl(p, vec3(0.,2.,diameter));
    float d = cyl(p, vec3(0.,2.,diameter));
    return max(c, -d);
}

float scene(vec3 p) {
    return min(plane(p + vec3(1.0,0.,0.)), cylinder(p));
    //return cylinder(p);
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
    vec3 ro = vec3(0., 2., 0.); // ray origin, here also known as 'eye'
    vec3 lookat = vec3(0.,0.,5.);
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(0.,.5,0.);
    //l1.position = vec3(1.2, 2.0, -1.9);
    color = vec3(.0);
    sky = color;

    // raymarch a scene
    float step = 1.;
    float d = 1.;
    vec3 p;
    for (int i = 0; i < MAXSTEPS; i++) {
        if ((d<eps) || (step > MAXDIST)) break;
        p = ro+rd*step;
        d = scene( p );
        step += d;
    }
    if( step>MAXDIST ) step=-1.;

    // apply materials
    if (step > 0.) {
        // precompute stuff
        vec3 n = normal(p);
        vec3 ld = normalize(l1.position-p);
        vec3 ed = normalize(ro-p);

        if (plane(p) < eps) {
            float checker = mod(floor(p.x)+floor(p.z)-.5, 2.0);
            m1.color = vec3(checker);
            m1.reflection_ratio = 0.01;
            m1.shininess = 3.;
            s1 = get_shading(m1, l1, p, n, ld, ed);
            color = m1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += m1.color * s1.amb *.2;
        }

        if (cylinder(p) < eps) {
            m1.color = vec3(1.,1.,0.);
            m1.reflection_ratio = 0.01;
            m1.shininess = 30.;
            s1 = get_shading(m1, l1, p, n, ld, ed);
            color = m1.color * s1.diffuse * s1.shadow/2.;
            color += s1.specular;
            color *= s1.aoc;
            color += m1.color * s1.amb *.2;
        }
    }

    // Exponential distance fog
    color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * step * step));

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(color, 1.);
}
