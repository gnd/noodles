#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform sampler2D sony;
uniform float time;
#define MAXSTEPS 128
#define MAXDIST 20.0
#define eps 0.0001
#define HEIGHT 10.

struct march { float step; vec3 ro; vec3 rd; };
struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material mat1;
shading s1;

bool refracted;
vec3 sky, color; // 'sky color'
vec3 p, p_refr;
vec2 sp = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
float zpos = time*3.;
vec3 pos = vec3( cos(time/4.)*4.-zpos, 0., sin(time/4.)*4. - zpos);

// sony height mapping
vec3 texcol(vec3 p) {
    p.x = mod(p.x, 18.) - 9.;
    p.z = mod(p.z, 8.) - 4.;
    vec3 pixel = texture2D(sony, (vec2(9.-p.x,4.-p.z))*vec2(.05,.1)).rgb;
    return vec3(dot(pixel, vec3(0.299, 0.587, 0.114)));
}

// sony height mapping + aliasing
vec3 texcol_aa(vec3 p) {
    p.x = mod(p.x, 18.) - 9.;
    p.z = mod(p.z, 8.) - 4.;
    vec3 col1 = texture2D(sony, (vec2(9.-p.x,4.-p.z))*vec2(.05,.1)).rgb;
    vec3 col2 = texture2D(sony, (vec2(9.-p.x+1.,4.-p.z))*vec2(.05,.1)).rgb;
    vec3 col3 = texture2D(sony, (vec2(9.-p.x,4.-p.z+1.))*vec2(.05,.1)).rgb;
    vec3 col4 = texture2D(sony, (vec2(9.-p.x+1.,4.-p.z+1.))*vec2(.05,.1)).rgb;
    vec3 pixel= mix(mix(mix(col3, col4,.5), col3, .5), col1, .5);
    return vec3(dot(pixel, vec3(0.299, 0.587, 0.114)));
}

vec3 texcol_real(vec3 p) {
    return texture2D(sony, sp).rgb;
}

vec3 txc(vec3 p) {
    //return texcol(p);
    //return texcol_aa(p);
    return texcol_real(p);
}

float sphere(vec3 p, float r) {
    return length(p) - r;
}

float ref(vec3 p) {
    return sphere(p-vec3(pos.x-10.,10.,pos.z-10.), 6.);
}

float plane(vec3 p) {
    return p.y - texcol(p).r * HEIGHT - txc(p).r * .2;
}

float scene(vec3 p) {
    return min(ref(p),plane(p));
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

void main() {
    vec3 ro = vec3(cos(time/4.)*4.-zpos, sin(time) + 10., sin(time/4.)*4. - zpos); // ray origin, here also known as 'eye'
    vec3 lookat = vec3(0.-zpos-20.,0.,-zpos - 20.);
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.)*4.-zpos,10.5,cos(time/4.)*3.-zpos + 2.);
    //l1.position = vec3(1.2, 2.0, -1.9);
    color = vec3(.0);
    sky = color;

    // raymarch a scene
    float step = .0;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = ro + (rd * step);
        float d = scene(p);
        if (d > MAXDIST) {
            break;
        }
        if (d < eps) {
            // precompute stuff
            vec3 n = normal(p);
            vec3 ld = normalize(l1.position-p);
            vec3 ed = normalize(ro-p);
            if (ref(p) < eps) {
                // Compute refraction
                p_refr = p;
                refracted = true;
                march march = refraction(p, n, rd, 1./2.22);
                step = march.step;
                rd = march.rd;
                ro = march.ro;
            } else {
                if (plane(p) < eps) {
                    mat1.color = texcol(p);
                    mat1.reflection_ratio = 0.01;
                    mat1.shininess = 3.;
                    s1 = get_shading(mat1, l1, p, n, ld, ed);
                    color = mat1.color * s1.diffuse * s1.shadow *10.;
                    color += mat1.color * s1.amb * .05;
                }
                if (p.y > 0.) {
                    color = vec3(.1*p.y,0.,0.);
                    color *= brightcon(texcol(p), 2., 3.);
                }
                if (p.y > HEIGHT/3.) {
                    color = vec3(0.,.08*p.y,0.);
                    color *= brightcon(texcol(p), 2., 3.);
                }
                if (p.y > (HEIGHT/3.)*1.5) {
                    color = vec3(0.,0.,.2*p.y);
                    color *= brightcon(texcol(p), 2., 3.);
                }
                color += vec3(float(i)*1./float(MAXSTEPS/10),0.,0.);
                break;
            }
        }
        step += d;
    }

    // Exponential distance fog
    color = mix(color, sky, 1.0 - exp2(-0.05 * step));

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(color, 1.);
}
