#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
uniform sampler2D backbuffer;

// PARAMS
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;

/////////////////// MIDIMIX //////////////////////
uniform float m5,m6,m7,m8,m9,m10;
uniform float m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30;
uniform float m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50;
uniform float m51,m52,m53,m54,m55,m56,m57,m58,m59;
float mx11 = m5/127.;
float mx12 = m6/127.;
float mx13 = m7/127.;
float mx14 = m8/127.;
float mx15 = m9/127.;
float mx16 = m10/127.;
float mx21 = m11/12.7;
float mx22 = m12/127.;
float mx23 = m13/127.;
float mx24 = m14/127.;
float mx25 = m15/127.;
float mx26 = m16/127.;
float mx31 = m17/127.;
float mx32 = m18/127.;
float mx33 = m19/127.;
float mx34 = m20/127.;
float mx35 = m21/127.;
float mx36 = m22/127.;
float mx41 = m23/127.;
float mx42 = m24/127.;
float mx43 = m25/127.;
float mx44 = m26/127.;
float mx45 = m27/127.;
float mx46 = m28/127.;
float mx51 = m29/127.;
float mx52 = m30/127.;
float mx53 = m31/127.;
float mx54 = m32/127.;
float mx55 = m33/127.;
float mx56 = m34/127.;
float mx61 = m35/127.;
float mx62 = m36/127.;
float mx63 = m37/127.;
float mx64 = m38/127.;
float mx65 = m39/127.;
float mx66 = m40/127.;
float mx71 = m41/127.;
float mx72 = m42/127.;
float mx73 = m43/127.;
float mx74 = m44/127.;
float mx75 = m45/127.;
float mx76 = m46/127.;
float mx81 = m47/127.;
float mx82 = m48/127.;
float mx83 = m49/127.;
float mx84 = m50/127.;
float mx85 = m51/127.;
float mx86 = m52/127.;
float mx91 = m53/127.;
float mx92 = m54/127.;
float mx93 = m55/127.;
float mx94 = m56/127.;
///////////////// END OF MIDIMIX ////////////////////


#define MAXSTEPS 96
#define MAXDIST 40.0
#define eps 0.001   // best results with 0.01 - lower seems to be too granullar

struct march { float step; vec3 ro; vec3 rd; };
struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material mt1,mt2,mt3;
shading s1,s2,s3;
vec2 p2;

bool reflected = false;
bool refracted = false;

// bumpsize is for the iridescent sphere
// bumpsize2 is for the object
// bumpsize3 is for the reflector
float bumpsize = mx16*10.;
float bumpsize2 = mx26*10. + m0/10.;
float bumpsize3 = mx36 + m0/10.; 
float miditime = time * 10. * mx91;

vec3 sky, color;
vec3 red = vec3(1.,0.,0.);
vec3 green = vec3(0.,1.,0.);
vec3 blue = vec3(0.,0.,1.);
vec3 yellow = vec3(1.,1.,0.);
vec3 purple = vec3(1.,0.,1.);
vec3 white = vec3(1.);
vec3 grey = vec3(.7);
vec3 black = vec3(0.);

// Counter 0 to 1
float cnt(in int m)
{
    float divider = 1000.0 / float(m)*10.0;
    return mod(time*divider, 10.0) / 10.0;
}

// Sinusoid bumps
float sinbumps(in vec3 p, in float frq) {
    return sin(p.x*frq) * sin(p.y*frq+miditime) * sin(p.z*frq) * .5;
}

// Sinusoid bumps - just y axis
vec3 ybumps(in vec3 p, in float frq) {
    return vec3(0., sin(p.x*frq) * sin(p.z*frq+miditime) * 0.5, 0.);
}

float smin( float a, float b, float k ) {
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

float vmax(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

// Cheap Box: distance to corners is overestimated
// from http://mercury.sexy/hg_sdf/
float cbox(vec3 p, vec3 b) { //cheap box
	return vmax(abs(p) - b);
}

float box( vec3 p, vec3 b ) {
    return length(max(abs(p)-b,0.0));
}

float iterbox_(vec3 p) {
    float res = 100.;
    for (int m=0; m<3; m++) {
        vec3 shift = vec3(0.41,.0,.0) * m;
        vec3 size = vec3(0.2);
        res = min (
            res,
            box(p - shift, size)
        );
    }
    return res;
}

float iterbox(vec3 p) {
    return cbox(p, vec3(0.2));
}

float plane(vec3 p) {
    return cbox(p+vec3(0.,.5,.0), vec3(15.,.01,15.));
}

float box(vec3 p) {
    return cbox(p, vec3(1.));
}

float sdCross( in vec3 p ) {
  float da = cbox(p.xyz,vec3(10.,1.0,1.0));
  float db = cbox(p.yzx,vec3(1.0,10.,1.0));
  float dc = cbox(p.zxy,vec3(1.0,1.0,10.));
  return min(da,min(db,dc));
}

float sdCross_smooth( in vec3 p ) {
  float da = cbox(p.xyz,vec3(2.,0.01,0.01));
  float db = cbox(p.yzx,vec3(0.01,2.,0.01));
  float dc = cbox(p.zxy,vec3(0.01,0.01,2.));
  float smoothness = mx81*10.;
  return smin(da,smin(db,dc,smoothness),smoothness);
}

float ssub( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h); 
}

#define R(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)

// Sphere field
float sphere(in vec3 p, in vec3 centerPos, float radius) {
        return length(p - centerPos) - radius;
}

float hollow_sphere(vec3 p, float bumps) {
    float size = mx11 * 20. + m1/10.;
    float d = sphere(p, vec3(0.), size) + bumps;
    float c = sphere(p, vec3(0.), size - 0.5) + bumps;
    return max(d,-c);
}

// Cube-oid refractive object
// TODO - rotation bumped by sound
float refractor(vec3 p, float bumps) {
    // rotate along the Y axis
    vec2 xz = R(p.xz, time);
    p = vec3(xz.x, p.y, xz.y);

    // thickness from -.1 ==> 1
    float thick = mx24 * 10.1 - .1;
    float d = cbox(p, vec3(1.,thick,1.)) + bumps;
    float e = cbox(p, vec3(thick,1.,1.)) + bumps;
    float f = cbox(p, vec3(1.,1.,thick)) + bumps;

    return min(min(d,e),f);
}

float reflector(vec3 p, vec3 bumps) {
    // plane
    p+=bumps;
    float size = mx31 * 20.;
    float position = (1. - mx35) * 10.;
    float upper = cbox(p + vec3(0.,-position,.0), vec3(size, mx34, size));
    float lower = cbox(p + vec3(0.,position,.0), vec3(size, mx34, size));
    return min(upper, lower);
}

float scene(vec3 p) {
    float bumps = sinbumps(p, bumpsize);
    float bumps2 = sinbumps(p, bumpsize2);
    vec3 bumps3 = ybumps(p, bumpsize3);

    float d = refractor(p*mx21, bumps2)/mx21;
    float c = reflector(p, bumps3);
    return min(c, min(hollow_sphere(p, bumps), d));
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p) {
    const vec2 k = vec2(1.,-1.);
    return normalize( k.xyy*scene( k.xyy*eps + p ) +
                      k.yyx*scene( k.yyx*eps + p ) +
                      k.yxy*scene( k.yxy*eps + p ) +
                      k.xxx*scene( k.xxx*eps + p ) );
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

march refraction(vec3 p, vec3 rd, float eta) {
    // Glass becomes a 'portal' through which refracted rays keep on raymarching
    march march;
    vec3 n = normal(p);
    march.ro = p;
    march.rd = normalize(refract(rd, n, eta));
    march.step = 0.0;
    for (int j=0; j < MAXSTEPS; j++) {
        p = march.ro + march.rd * march.step;
        float d = scene(p);
        march.step += max(abs(d),eps);
        if (d > eps) {
            // attenuate due to impurities, etc
            color *= .1;
            // second refraction
            march.ro = p;
            // lower means more pure
            float glass_purity = 0.47;
            march.rd = normalize(refract(march.rd, -normal(p), 1./eta*glass_purity));
            color = vec3(0.,0.,1.);
            break;
        }
    }
    return march;
}

march reflection(vec3 p, vec3 light_dir) {
    march march;
    vec3 n = normal(p);
    march.ro = p;
    march.rd = normalize(reflect(-light_dir, n));
    march.step = 0.0025;
    for (int j=0; j < MAXSTEPS; j++) {
        p = march.ro + march.rd * march.step;
        float d = scene(p);
        march.step += max(abs(d),eps);
        if (d > eps) {
            break;
        }
    }
    return march;
}

shading get_shading(material m, light l, vec3 p, vec3 n, vec3 ld, vec3 ed) {
    shading s;
    float step;
    vec3 light_dir = l.position - p; // direction from p to light
    vec3 normal = normal(p); // surface normal at point p

    // Diffuse shading using Lambertian reflectance
    // https://en.wikipedia.org/wiki/Lambertian_reflectance
    s.diffuse = clamp(dot(normal,normalize(light_dir)), 0., 1.);

    // Specular reflections using phong
    // https://en.wikipedia.org/wiki/Phong_reflection_model
    vec3 reflection_dir = normalize(reflect(-light_dir, normal));
    s.specular = clamp(m.reflection_ratio * pow(dot(reflection_dir, normalize(ed-p)), m.shininess), 0., 1.);

    // Shadows - basically raymarch towards the light
    // If the path from p to each light source is obstructed add shadow
    // soft shadows taken from http://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
    s.shadow = 1.;
    step = 0.01;
    float ph = 1e10;
    for (int i = 0; i < MAXSTEPS/2; i++ ) {
        vec3 shadow_p = p + normalize(light_dir) * step;
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
        s.shadow = min(s.shadow, 32.*d/max(.0,step-y));
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
        vec3 aoc_p = p + normal * step;
        float diff = 1.5*abs(length(aoc_p-p) - scene(aoc_p));
        s.aoc = min(s.aoc, s.aoc-diff);
    }

    // Ambient light
    s.amb = clamp(0.5+0.4*normal.y, 0., 1.);

    return s;
}

vec3 str(in float xp, in float yp, in float x, in float y, in vec3 c, in vec3 col) {
        vec3 cl = ((p2.x > xp-x*.51) && (p2.x < x*.51+xp) && (p2.y > yp-y*.51) && (p2.y < y*.51+yp)) ? -c+col : vec3(0.0);
        return cl;
}

// Feedback
vec3 feedb_sqr(in float xpos, in float ypos, in float xsiz, in float ysiz, in float bsiz, in vec3 c) {
    vec3 ccc = texture2D(backbuffer, (vec2(p2.x,p2.y)-0.51)*bsiz*1.5+0.51).xyz;
    c+=str(xpos,ypos,xsiz,ysiz,c,ccc);
    return c;
}

void main() {
    // TODO - bump camera y according to sound
    vec3 ro = vec3(0.,m0*10., 5.);
    //ro = vec3(cos(t/4.)*20.*mx81, 1.5, sin(t/4.)*20.*mx81);
    //ro = vec3(0.,0.,0.);
    vec3 lookat = vec3(0.,0.,0.);
    //vec3 lookat = vec3(cos(time/4.)*2.,sin(time/10.),sin(time/4.)*2.);
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    sky = black;
    color = sky;

    // raymarch a scene
    float step = .0;
    vec3 p, p_refr, p_refl;
    float eta;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = rd * step + ro;
        float d = scene(p);
        if (d > MAXDIST) {
            if (refracted) {
                mt1.color = vec3(1.,1.,1.);
                mt1.reflection_ratio = 10.9;
                mt1.shininess = 100.9;
                vec3 n = normal(p);
                s1 = get_shading(mt1, l1, p, n, p_refr, ro);
                color += s1.specular*10.;
            }
            if (reflected) {
                // attenuate sky color
                color *= 0.1;
            }
            break;
        }
        if (d < eps) {
            // precompute stuff
            vec3 n = normal(p);
            vec3 ld = normalize(l1.position-p);
            vec3 ed = normalize(ro-p);
            float bumps = sinbumps(p, bumpsize);
            float bumps2 = sinbumps(p, bumpsize2);
            vec3 bumps3 = ybumps(p, bumpsize3);
            
            if (refractor(p*mx21, bumps2)/mx21 < eps) {
                // Compute refraction
                    p_refr = p;
                    refracted = true;
                    eta = 1./2.22;
                    march march = refraction(p, rd, eta);
                    step = march.step;
                    rd = march.rd;
                    ro = march.ro;
                    color = vec3(0.,0.,0.);
            } else if (reflector(p, bumps3) < eps) {
                // Compute reflection
                    p_refl = p;
                    reflected = true;
                    vec3 ld = normalize(ro - p);
                    march march = reflection(p, ld);
                    step = march.step;
                    rd = march.rd;
                    ro = march.ro;
            } else  {
                if ((hollow_sphere(p, bumps) < eps) && (mx11 > 0.)) {
                    if (mx12 == 0) {
                        mt1.color = vec3(1.,1.,.0);
                        mt1.reflection_ratio = .1;
                        mt1.shininess = 10.9;
                        s1 = get_shading(mt1, l1, p, n, ld, ro);
                        // check https://www.shadertoy.com/view/lsKcDD
                        color = mt1.color * s1.diffuse * s1.shadow;
                        color += s1.specular;
                        color *= s1.aoc;
                        color += mt1.color * s1.amb *.09;
                    }
                    else if (mx12 == 1) {
                        float gamma = 0.7 - (mx11 * 0.6);
                        float shininess = 0.7 - (mx11 * 0.5);
                        vec3 perturb = cos(p * 30.);
                        mt1.color = spectrum( dot(perturb * .05 + n, ro) * 2.);
                        s1 = get_shading(mt1, l1, p, n, ld, ed);
                        color = gamma - mt1.color;
                        color += clamp(shininess * pow(dot(normalize(reflect(-ld, n)), ed), 7.), .0, 1.);
                        color = pow( color, vec3(1.0/2.2) );
                    }
                }

                if (refracted) {
                    mt1.color = vec3(10.,10.,1.);
                    mt1.reflection_ratio = 10.9;
                    mt1.shininess = 100.9;
                    s1 = get_shading(mt1, l1, p, n, p_refr, ro);
                    color += s1.specular*10.;
                }
                if (reflected) {
                    // attenuate reflected color
                    color *= 0.8;
                }
                 break;
            }
        }
        step += d;
    }

    // Exponential distance fog
    //color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * step * step));

    // gamma correction
    //color = pow( color, vec3(1.0/2.2) );

    // Strobe
    //PixelColor += cnt(100);

    // set coords for feedback
    float asp = resolution.x / resolution.y;
    p2 = gl_FragCoord.xy / resolution;

    // feedback normal
    if (mx72 > 0.) {
        vec3 col = feedb_sqr(.5, .5, 1., 1., mx74, color.xyz)*mx71*1.1;
        color += col;
    }

    PixelColor = vec4(color, 1.0);
}