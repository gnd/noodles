#version 450
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time; // number of seconds since start as float: 1.234, 1.235, etc
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
uniform float m51,m52,m53,m54,m55,m56;
float mx11 = m5/127.;
float mx12 = m6/127.;
float mx13 = m7/127.;
float mx14 = m8/127.;
float mx15 = m9/127.;
float mx16 = m10/127.;
float mx21 = m11/127.;
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

/// Speeds
uniform float m57,m58,m59,m60,m61,m62,m63,m64,m65,m66;
float sp0 = m57;
float sp1 = m58;
float sp2 = m59;
float sp3 = m60;
float sp4 = m61;
float sp5 = m62;
float sp6 = m63;
float sp7 = m64;
float sp8 = m65;
float sp9 = m66;

#define MAXSTEPS 128
#define MAXDIST 192
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
// bumpsize2 is for the refractor
// bumpsize3 is for the reflector
float bumpsize  = m1 * 1.5 + (mx16*10.);
float bumpsize2 = m2 * 0.5 + (mx26*10.);
float bumpsize3 = m0 * 0.1 + mx36; 
float miditime  = (time * 10.) * mx61;

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
// time is the number of sec.ms from beginning: 1.234, 1.235, .. 
float cnt(in int m) {
    float divider = 10000.0 / float(m);
    return mod(time*divider, 10.0) / 10.0;
}

// Random number generator
// lumina.sourceforge.net/Tutorials/Noise.html
float rnd(vec2 co){
    float k = sin(dot(co.xy ,vec2(12.9898,78.233)));
    return fract(k + k);
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

//==============================================================
// otaviogood's noise from https://www.shadertoy.com/view/ld2SzK
//--------------------------------------------------------------
// This spiral noise works by successively adding and rotating sin waves while increasing frequency.
// It should work the same on all computers since it's not based on a hash function like some other noises.
// It can be much faster than other noise functions if you're ok with some repetition.
const float nudge = 4.; // size of perpendicular vector
float normalizer = 1.0 / sqrt(1.0 + nudge*nudge);   // pythagorean theorem on that perpendicular to maintain scale
float SpiralNoiseC(vec3 p) {
    //float n = -mod(time * 0.2,-2.); // noise amount
    float n = mx64 * -7. + 5.; // noise amount
    float iter = 2.0;
    for (int i = 0; i < 8; i++)
    {
        // add sin and cos scaled inverse with the frequency
        n += -abs(sin(p.y*iter) + cos(p.x*iter)) / iter;    // abs for a ridged look
        // rotate by adding perpendicular and scaling down
        p.xy += vec2(p.y, -p.x) * nudge + time;
        p.xy *= normalizer;
        // rotate on other axis
        p.xz += vec2(p.z, -p.x) * nudge;
        p.xz *= normalizer;
        // increase the frequency
        iter *= 1.733733;
    }
    return n;
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

// rotation
// dont use - this modifies all of p
#define R(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)

// rotation as a function
vec2 rot2d(in vec2 p, in float angle) {
    return cos(radians(angle))*p + sin(radians(angle))*vec2(p.y, -p.x);
}

// precompute stuff for the tunnel 
// its actually two planes with a noise heightmap
float speed_z = mx41 * 100 + abs(m1 * 0.01);
float speed_x = speed_z * .1 + 5.;
vec2 path2d = vec2(cos(time * 0.33) * speed_x, speed_z * time);
vec3 path = vec3( path2d.x, 2.-noise2f(path2d.xy/10.)*10., path2d.y );

float low_plane(vec3 p) {
    // bend
    float snd = 0;
    if (mx43 > 0.) {
        snd = m4 * 0.3;
    }
    float k = 0.09 * mx44 + 0.002 + snd;
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2 m = mat2(c,-s,s,c);
    vec3 q = vec3(m*p.xy,p.z);
    return noise2f(q.xz*0.02)*(20.+m5) + noise2f(q.xz*0.1)*10. + noise2f(q.xz*0.5)*2. + q.y + 15.;
}

float sky_plane(vec3 p) {
    // bend
    float snd = 0;
    if (mx43 > 0.) {
        snd = m4 * 0.3;
    }
    float k = -0.09 * mx44 - 0.002 - snd;
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2 m = mat2(c,-s,s,c);
    vec3 q = vec3(m*p.xy,p.z);
    return noise2f(q.zx*0.02)*(20.+m5) + noise2f(q.zx*0.1)*10. + noise2f(q.zx*0.5)*2. - q.y + 15.;
}

float back_plane(vec3 p) {
    p = vec3(p.x,p.y+sp1*mx41,p.z*.7);
    return 20. + p.z + noise2f(p.xy) + noise2f(p.xy/2.)*2. + noise2f(p.xy/3.)*3. + noise2f(p.xy/1.)*2.;
}

float planes(vec3 p) {
    //vec3 q = p - path;
    //return min(low_plane(q), sky_plane(q));
    return back_plane(p);
}

// Cheap Box: distance to corners is overestimated
// from http://mercury.sexy/hg_sdf/
float cbox(vec3 p, vec3 b) { //cheap box
	return vmax(abs(p) - b);
}

float plane(vec3 p) {
    return cbox(p+vec3(0.,.5,.0), vec3(15.,.01,15.));
}

// Sphere field
float sphere(in vec3 p, in vec3 centerPos, float radius) {
        return length(p - centerPos) - radius;
}

float hollow_sphere(vec3 p, float bumps) {
    float size = mx11 * 50. + m1 * 0.2 + 0.01;
    float d = sphere(p, vec3(0.), size) + bumps;
    float c = sphere(p, vec3(0.), size - 1.) + bumps;
    return max(d,-c);
}

// Cube-oid refractive object
float refractor(vec3 p, float bumps) {
    // rotate along the Y axis
    vec3 q1 = p;
    vec2 xz = rot2d(q1.xz, sp0*360.);
    q1 = vec3(xz.x, q1.y, xz.y);

    // rotate along the X axis
    vec3 q2 = p;
    vec2 yz = rot2d(q2.yz, sp0*360.+90.);
    q2 = vec3(q2.x, yz.x, yz.y);

    // rotate along the Z axis
    vec3 q3 = p;
    vec2 xy = rot2d(q3.xy, sp0*360.0-90.);
    q3 = vec3(xy.x, xy.y, q3.z);

    // thickness from -.1 ==> 1
    float thick = mx24 * 8.1 - .1;
    float size = mx21 + sp2;
    float d = cbox(q2/size*0.2, vec3(1.,thick,1.))*size*5. + bumps;
    float e = cbox(q3/size*0.2, vec3(thick,1.,1.))*size*5. + bumps;
    float f = cbox(q1/size*0.2, vec3(1.,1.,thick))*size*5. + bumps;

    return smin(smin(d,e,2.),f,1.0);
}

float reflector(vec3 p, vec3 bumps) {
    // plane
    p+=bumps;
    float size = mx31 * 50.;
    float position = (1. - mx35) * 20.;
    float upper = cbox(p + vec3(0.,-position,.0), vec3(size, mx34, size));
    float lower = cbox(p + vec3(0.,position,.0), vec3(size, mx34, size));
    return min(upper, lower);
}

float scene(vec3 p) {
    // pre-compute bumps
    float bumps = sinbumps(p, bumpsize);
    float bumps2 = sinbumps(p, bumpsize2);
    vec3 bumps3 = ybumps(p, bumpsize3);

    // max distance field to loose it in min() when off
    float hollow_d = 100.;
    float refractor_d = 100.;
    float reflector_d = 100.;
    float planes_d = 100.;

    if (mx12 > 0.) {
        hollow_d = hollow_sphere(p, bumps);
    }

    if (mx22 > 0.) {
        refractor_d = refractor(p, bumps2);
    }

    if (mx32 > 0.) {
        reflector_d = reflector(p, bumps3);
    }
    
    if (mx42 > 0.) {
        planes_d = planes(p);
    }

    return min(min(min(hollow_d, refractor_d), reflector_d), planes_d);
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

float fogFactorExp2(const float dist, const float density) {
  const float LOG2 = -1.442695;
  float d = density * dist;
  return 1. - clamp(exp2(d * d * LOG2), 0.0, 1.0);
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
            float glass_purity = .47;
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
    for (int i = 0; i < MAXSTEPS; i++ ) {
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
    vec3 ro;
    if (mx62 > 0.) {
        ro = vec3(5.*sin(m2),1., 5.*cos(m2)+5.);    
    } else {
        ro = vec3(0.,0.,5.);
    }
    
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
    l1.position = vec3(sin(time/4.)*4.,2.5,10.);
    //l1.position = vec3(0.);
    sky = black;
    color = sky;

    // raymarch a scene
    float sdf_step = .0;
    vec3 p, p_refr, p_refl;
    float eta;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p = rd * sdf_step + ro;
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
            
            if ((refractor(p, bumps2) < eps) && (mx22 > 0.)) {
                // Compute refraction
                    p_refr = p;
                    refracted = true;
                    eta = 1./2.22;
                    march march = refraction(p, rd, eta);
                    sdf_step = march.step;
                    rd = march.rd;
                    ro = march.ro;
                    color = vec3(0.,0.,0.);
            } else if ((reflector(p, bumps3) < eps) && (mx32 > 0.)) {
                // Compute reflection
                    p_refl = p;
                    reflected = true;
                    vec3 ld = normalize(ro - p);
                    march march = reflection(p, ld);
                    sdf_step = march.step;
                    rd = march.rd;
                    ro = march.ro;
            } else  {
                if ((hollow_sphere(p, bumps) < eps) && (mx12 > 0.)) {
                    if (mx13 == 0) {
                        mt1.color = vec3(1.,1.,.0);
                        mt1.reflection_ratio = .1;
                        mt1.shininess = 10.9;
                        s1 = get_shading(mt1, l1, p, n, ld, ro);
                        // check https://www.shadertoy.com/view/lsKcDD
                        color = mt1.color * s1.diffuse * s1.shadow;
                        color += s1.specular;
                        color *= s1.aoc;
                        color += mt1.color * s1.amb *.09;
                        color = pow( color, vec3(1.0/4.2) );
                    }
                    else if (mx13 == 1) {
                        float gamma = 0.5 - (mx11 * 0.6);
                        float shininess = 0.7 - (mx11 * 0.5);
                        vec3 perturb = cos(p * 30.);
                        mt1.color = spectrum( dot(perturb * .05 + n, ro) * 2.);
                        s1 = get_shading(mt1, l1, p, n, ld, ed);
                        color = gamma - mt1.color;
                        color += clamp(shininess * pow(dot(normalize(reflect(-ld, n)), ed), 7.), .0, 1.);
                        color = pow( color, vec3(1.0/2.2) );
                    }
                }
                if ((planes(p) < eps) && (mx42 == 1)) {
                    color += vec3(float(i)*(1./float(MAXSTEPS/2.)),0.,0.);
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
                    if (mx33 > 0) {
                        color *= 1.9;
                    } else {
                        color *= .6;
                    }
                    
                }
                 break;
            }
        }
        sdf_step += d;
    }

    // Exponential distance fog
    color = mix(color, color * 0.05, fogFactorExp2(sdf_step, 0.005));

    // gamma correction
    //color = pow( color, vec3(1.0/2.2) );

    // set coords for feedback
    float asp = resolution.x / resolution.y;
    p2 = gl_FragCoord.xy / resolution;

    // feedback normal
    if (mx72 > 0.) {
        vec3 col = feedb_sqr(.5, .5, 1., 1., mx74, color.xyz)*mx71*1.1;
        color += col;
    }

    // channel mix
    color = vec3(color.r*mx86,color.g*mx85,color.b*mx84);

    // strobe
    color += step(0.75,vec3(cnt(50))) * mx83*.9;

    // permanent strobe
    color += step(0.98,cnt(100)) * mx82;

    PixelColor = vec4(color, 1.0) * mx81;
}