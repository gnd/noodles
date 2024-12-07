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

/// Control data
uniform float m57,m58,m59,m60,m61,m62,m63,m64,m65,m66;
uniform float m67,m68,m69,m70,m71,m72,m73,m74,m75,m76;
float cx0 = m57;
float cx1 = m58;
float cx2 = m59;
float cx3 = m60;
float cx4 = m61;
float cx5 = m62;
float cx6 = m63;
float cx7 = m64;
float cx8 = m65;
float cx9 = m66;
float cx10 = m57;
float cx11 = m58;
float cx12 = m59;
float cx13 = m60;
float cx14 = m61;
float cx15 = m62;
float cx16 = m63;
float cx17 = m64;
float cx18 = m65;
float cx19 = m66;

#define MAXSTEPS 64
#define MAXDIST 128
#define eps 0.001   // best results with 0.01 - lower seems to be too granullar
#define SPHERES 18

struct march { float step; vec3 ro; vec3 rd; };
struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material mt1,mt2,mt3;
shading s1,s2,s3;
vec4 spheres[SPHERES];

bool reflected = false;
bool refracted = false;

float bumpsize  = m0 + mx16*10.;
float bumpsize2 = m1 * 0.5;
float bumpsize3 = m2 * 0.1; 
float miditime  = time * 10.;

vec3 sky, color;
vec3 red = vec3(1.,0.,0.);
vec3 green = vec3(0.,1.,0.);
vec3 blue = vec3(0.,0.,1.);
vec3 yellow = vec3(1.,1.,0.);
vec3 purple = vec3(1.,0.,1.);
vec3 white = vec3(1.);
vec3 grey = vec3(.7);
vec3 black = vec3(0.);

// set coords for feedback
vec2 p2 = gl_FragCoord.xy / resolution;

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

// Sinusoid bumps set freq
float sinbumps_set(in vec3 p){
    float frq = 1.7;
    return sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31)*sin(time)*1.7;
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

vec3 blend(in vec3 d1, in vec3 d2, in float k) {
    return k * d1 + (1.-k) * d2;
}

float blend(in float d1, in float d2, in float k) {
    return k * d1 + (1.-k) * d2;
}

// rotation as a function
vec2 rot2d(in vec2 p, in float angle) {
    return cos(radians(angle))*p + sin(radians(angle))*vec2(p.y, -p.x);
}

float back_plane(vec3 p, vec3 bumps) {
    p = vec3(p.x,p.y+cx1*mx41,p.z*.7) + bumps;
    return 20. + p.z + noise2f(p.xy) + noise2f(p.xy/2.)*2. + noise2f(p.xy/3.)*3. + noise2f(p.xy/1.)*2.*m0;
}

float planes(vec3 p, vec3 bumps) {
    return back_plane(p, bumps);
}

// Cheap Box: distance to corners is overestimated
// from http://mercury.sexy/hg_sdf/
float cbox(vec3 p, vec3 b) { //cheap box
	return vmax(abs(p) - b);
}

float plane(vec3 p) {
    return cbox(p+vec3(0.,.5,.0), vec3(15.,.01,15.));
}

float grid_plane(vec3 p, float bumps) {
    return p.y + 1. + bumps * .1;
}

float sphere(in vec3 p, in vec3 centerPos, float radius) {
    return length(p - centerPos) - radius;
}

float sphere_simple(vec3 p, float r) {
    return length(p) - r;
}

float lunch_sphere(in vec3 p, float bumps) {
    vec3 center = vec3(0., .5, 1.);
    return sphere(p, center, 1.*10.*m1) + bumps;
}

float lunch_sphere_hollow(in vec3 p, float bumps) {
    vec3 center = vec3(0., .5, 2.);
    float outer = sphere(p, center, 0.75 + mx11 * 50.) + bumps;
    float inner = sphere(p, center, mx11 * 50.) + bumps;
    return max(outer, -inner);
}

float low_plane(vec3 p) {
    // bend
    float k = 0.07; // or some other amount
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2 m = mat2(c,-s,s,c);
    p = vec3(m*p.xy,p.z);
    return p.y + 20.+ noise2f(p.xz/10.)*10. + noise2f(p.xz) + noise2f(p.xz*5.);
}

float sky_plane(vec3 p) {
    // bend
    float k = -0.09; // or some other amount
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2 m = mat2(c,-s,s,c);
    p = vec3(m*p.xy,p.z);
    return -p.y + 20. - noise2f(p.xz/10.)*10. + noise2f(p.xz) + noise2f(p.xz*5.);
}

// taken from https://www.shadertoy.com/view/4ttGDH
float sinfield (vec3 p) {
    // camera orientation
    vec3 ro;
    if (mx52 > 0.) {
        ro = vec3(5.*sin(m2),1., 5.*cos(m2)+5.);    
    } else {
        ro = vec3(0.,0.,5.);
    }

    p.xy -= ro.xy;

    // precompute stuff
    float speed = pow(1.01,100.);
    float speed_z = 0.1+(100.*mx15);
    float speed_x = speed_z * .5;

    float camera_rot = 0.;
    vec2 path2d = vec2(cos(time/3.)*(speed_x+camera_rot), time*speed_z);
    vec2 path2d_next = vec2(cos(time+2./3.)*speed_x+camera_rot/2., (time+2.)*speed_z);

    // s1 / heightmap
    vec3 path = vec3( path2d.x*0.1, 2.-noise2f(path2d.xy/10.)*10., path2d.y );
    p -= path;


    p = cos(p*.315*1.25*(1.-mx56) + sin(p.zxy*.875*1.25*mx55*10.*m0));
    float bumps = abs(sinbumps_set(p))*mx54;
    float n = length(p);
    return (n - 1.025)*1.33;//-bumps;
}

float tunel (vec3 p) {
    // precompute stuff
    float speed = pow(1.01,100.);
    float speed_z = 0.1+(100.*mx15);
    float speed_x = speed_z * .5;

    float camera_rot = 0.;
    vec2 path2d = vec2(cos(time/3.)*(speed_x+camera_rot), time*speed_z);
    vec2 path2d_next = vec2(cos(time+2./3.)*speed_x+camera_rot/2., (time+2.)*speed_z);

    // s1 / heightmap
    vec3 path = vec3( path2d.x*0.1, 2.-noise2f(path2d.xy/10.)*10., path2d.y );
    p -= path;
    return min(low_plane(p), sky_plane(p));
}

// Cube-oid refractive object
float refractor(vec3 p, float bumps) {
    p += vec3(0.,-.5,0.);

    // rotate along the Y axis
    vec3 q1 = p;
    vec2 xz = rot2d(q1.xz, cx0*360);
    q1 = vec3(xz.x, q1.y, xz.y);

    // rotate along the X axis
    vec3 q2 = p;
    vec2 yz = rot2d(q2.yz, cx1*360.);
    q2 = vec3(q2.x, yz.x, yz.y);

    // rotate along the Z axis
    vec3 q3 = p;
    vec2 xy = rot2d(q3.xy, -90.);
    q3 = vec3(xy.x, xy.y, q3.z);

    float thick = 20.*mx65;
    float size = mx64*2.;
    float d = cbox(q2/size*0.2, vec3(1.,thick,1.))*size*5. + bumps;
    float e = cbox(q1/size*0.2, vec3(thick,1.,1.))*size*5. + bumps;
    float f = cbox(q3/size*0.2, vec3(1.,1.,thick))*size*5. + bumps;

    return smin(f,smin(d,e,2.1),2.1);
}

// from https://www.shadertoy.com/view/4t2GDG
float blob (vec3 p) {
    float res = 10000.;
    for( int i = 0; i < SPHERES; i++ ) {
        res = smin(res, sphere_simple(p-spheres[i].xyz, spheres[i].w), 8.);
   }
   return res;
}

// dummy
float reflector(vec3 p, vec3 bumps) {
    return 1.;
}

float scene(vec3 p) {
    // pre-compute bumps
    float bumps = sinbumps(p, bumpsize);
    float bumps2 = sinbumps(p, bumpsize2);
    vec3 plane_bumps = ybumps(p, 1.8);

    // max distance field to loose it in min() when off
    float lunch_sphere_d = MAXDIST;
    float grid_plane_d = MAXDIST;
    float planes_d = MAXDIST;
    float refractor_d = MAXDIST;
    float blob_d = MAXDIST;
    float tunel_d = MAXDIST;
    float sinfield_d = MAXDIST;


    if (mx12 > 0.) {
        lunch_sphere_d = lunch_sphere_hollow(p, bumps);
    }

    if (mx22 > 0.) {
        grid_plane_d = grid_plane(p, bumps2);
    }

    if (mx32 > 0.) {
        planes_d = planes(p, plane_bumps);
    }

    if (mx42 > 0.) {
        //refractor_d = refractor(p, bumps2);
        tunel_d = tunel(p);
    }

    if (mx52 > 0.) {
        blob_d = blob(p);
    }

    if (mx53 > 0.) {
        refractor_d = refractor(p, bumps2);
    }

    if (mx43 > 0.) {
        sinfield_d = sinfield(p);
    }

    return min(sinfield_d,min(refractor_d,min(tunel_d,min(blob_d, min(planes_d, min(lunch_sphere_d, grid_plane_d))))));
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

// Texture copies
vec3 texcopy(in float tx_posx, in float tx_posy, in float tx_zoom_x, in float tx_zoom_y) {
    float tx_size_x = 1. / tx_zoom_x;
    float tx_size_y = 1. / tx_zoom_y;
    vec2 coord = gl_FragCoord.xy / resolution;
    vec2 sample_coord = gl_FragCoord.xy / resolution - vec2(tx_posx-(1./tx_zoom_x*.5),tx_posy-(1./tx_zoom_y*.5));
    vec3 bg_buff = texture(backbuffer, vec2(sample_coord.x*tx_zoom_x,sample_coord.y*tx_zoom_y), 1.).xyz;
    return ((coord.x > tx_posx-tx_size_x*.5) && (coord.x < tx_posx+tx_size_x*.5) && (coord.y > tx_posy-tx_size_y*.5) && (coord.y < tx_posy+tx_size_y*.5)) ? bg_buff-color : vec3(0.0);
}

// Texture copies - physical modelling
vec3 texcopy_pm(in float tx_posx, in float tx_size_x, in float tx_posy, in float tx_size_y) {
    float tx_zoom_x = 1. / tx_size_x;
    float tx_zoom_y = 1. / tx_size_y;
    vec2 coord = gl_FragCoord.xy / resolution;
    vec2 sample_coord = gl_FragCoord.xy / resolution - vec2(tx_posx-(1./tx_zoom_x*.5),tx_posy-(1./tx_zoom_y*.5));
    vec3 bg_buff = texture(backbuffer, vec2(sample_coord.x*tx_zoom_x,sample_coord.y*tx_zoom_y), 1.).xyz;
    return ((coord.x > tx_posx-tx_size_x*.5) && (coord.x < tx_posx+tx_size_x*.5) && (coord.y > tx_posy-tx_size_y*.5) && (coord.y < tx_posy+tx_size_y*.5)) ? bg_buff-color : vec3(0.0);
}

void main() {
    // TODO - bump camera y according to sound
    vec3 ro;
    if (mx62 > 0.) {
        ro = vec3(5.*sin(cx1*m0),1., 5.*cos(cx1*m0)+5.);    
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


    // setup spheres
    if (mx52 > 0.) {
        vec4 tmp;
        float multi = m2+m1;
        for( int i = 0; i < SPHERES; i++ ) {
            tmp.x = cos(time *.13 * (float( i )+2.))*multi;
            tmp.y = sin(time * .075 * (float( i )+4.))*multi;
            tmp.z = sin(time * .1 * (float( i )+3.3))*m1;
            tmp.w = .1 * (sin(time * .1  *(float( i) +1.))+2.) + mx51*6.0;
            spheres[i] = tmp;
        }
    }

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
             vec3 plane_bumps = ybumps(p, 1.8);
            
            if ((refractor(p, bumps2) < eps) && (mx53 > 0.)) {
                // Compute refraction
                p_refr = p;
                refracted = true;
                eta = 1./2.22;
                march march = refraction(p, rd, eta);
                sdf_step = march.step;
                rd = march.rd;
                ro = march.ro;
                color = vec3(0.,0.,0.);
            } else if ((lunch_sphere(p, bumps) < eps) && (mx12 > 0.)) {
                // Compute refraction
                p_refr = p;
                refracted = true;
                eta = 1./2.22;
                march march = refraction(p, rd, eta);
                sdf_step = march.step;
                rd = march.rd;
                ro = march.ro;
                color = vec3(0.,0.,0.);
            } else if ((blob(p) < eps) && (mx52 > 0.)) {
                // Compute refraction
                p_refr = p;
                refracted = true;
                eta = 1./2.22;
                march march = refraction(p, rd, eta);
                sdf_step = march.step;
                rd = march.rd;
                ro = march.ro;
                color = vec3(0.,0.,0.);

            // reflector not beinh g used currently
            } else if ((reflector(p, bumps3) < eps) && (1. > 0.)) {
                // Compute reflection
                p_refl = p;
                reflected = true;
                vec3 ld = normalize(ro - p);
                march march = reflection(p, ld);
                sdf_step = march.step;
                rd = march.rd;
                ro = march.ro;
            } else  {
                if ((lunch_sphere(p, bumps) < eps) && (mx12 > 0.)) {
                    // rotate along the Y axis
                    vec3 q1 = p;
                    vec2 xz = rot2d(q1.xz, 10.*cx0*360.);
                    q1 = vec3(xz.x, q1.y, xz.y);

                    float checker = fract(floor(q1.y*10.)*.5+floor(q1.x*10.)*0.5)*2.+.5;
                    mt1.color = vec3(checker);
                    mt1.reflection_ratio = 0.01;
                    mt1.shininess = 3.;
                    s1 = get_shading(mt1, l1, p, n, ld, ed);
                    color = mt1.color * s1.diffuse * s1.shadow/20.;
                    color += s1.specular;
                    color *= s1.aoc;
                    color += mt1.color * s1.amb *0.6;
                }
                if ((grid_plane(p, bumps) < eps) && (mx22 > 0.)) {
                    float checker = mod(floor(p.x)+floor(p.z)-.5, 2.0);
                    mt1.color = vec3(checker);
                    mt1.reflection_ratio = 0.01;
                    mt1.shininess = 3.;
                    s1 = get_shading(mt1, l1, p, n, ld, ed);
                    color = mt1.color * s1.diffuse * s1.shadow/2.;
                    color += s1.specular;
                    color *= s1.aoc;
                    color += mt1.color * s1.amb *.2;
                }
                if ((blob(p) < eps) && (mx52 > 0.)) {
                        float nv = dot(n, -rd);
                        vec3 col = vec3(0.);
                        col += sin(nv * vec3(0.0, 1.0, 0.0) * 10.0 * 1.5) * 0.5 + 0.5;
                        col += sin(nv * vec3(1.0, 0.0, 0.0) * 20.0 * 1.5) * 0.5 + 0.5;
                        col += sin(nv * vec3(0.0, 0.0, 1.0) * 5.0 * 1.5) * 0.5 + 0.5;
                        col = 1.1 - col;
                        mt1.color = clamp(normalize(col), 0.0, 1.0);
                        float spec = clamp(3.5 * pow(dot(normalize(reflect(-ld, n)), ed), 70.), 0., 1.);
                        color += mt1.color * clamp(0.4*n.y+0.5, 0., 1.) * 3.;
                }
                if ((planes(p, plane_bumps) < eps) && (mx32 == 1)) {
                    color += vec3(float(i)*(1./float(MAXSTEPS/20.)),0.,0.)*m1;
                }
                if ((tunel(p) < eps) && (mx42 > 0.)) {
                    vec3 perturb = sin(p * mx14);
                    mt1.color = spectrum( dot(perturb * .05 + n, ro) * 2.);
                    //s1 = get_shading(m1, l1, p, n, ld, ed);
                    color = 0.01 - mt1.color;
                    color += clamp(.5 * pow(dot(normalize(reflect(-ld, n)), ed), .1), 0., 1.);
                    color += mt1.color * clamp(0.4*n.y+0.5, 0., 1.) * 10.5;
                }
                if ((sinfield(p) < eps) && (mx43 > 0.)) {
                    vec3 pp = mod(p*100., 2.) - .5;
                    vec3 i = normalize(ro - pp);
                    vec3 ir = refract(i, n, 1./1.);
                    vec3 ir2 = refract(i, -n, 1./2.);
                    vec3 ir_color2 = vec3(cross(ir2, -n));
                    // mix it
                    mt1.color = ir_color2 * clamp(dot(n,ld), 0., 1.);
                    mt1.color += ir_color2*1.3 + clamp(2.9 * pow(dot(normalize(reflect(-ld, n)), ed), 1.9), 0., 1.) * .1;
                    mt1.color += clamp(.5 * pow(dot(normalize(reflect(-ld, n)), ed), .1), 0., 1.) * 0.1;
                    mt1.color += vec3(cross(ir, n)) * clamp(.9 * pow(dot(normalize(reflect(-ld, n)), ed), .9), 0., 1.) * 2.5;
                    mt1.color += ir_color2 * clamp(0.4*n.y+0.5, 0., 1.) * 1.1;
                    color = mt1.color * smoothstep(0.,0.2, 1.);
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
    color = mix(color, color * 0.05, fogFactorExp2(sdf_step, 0.05));

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

    // physically modelled feedback layer
    if (mx73 > 0.) {
        color += texcopy_pm(cx5, cx6, cx7, cx8)*10.;
        //color += texcopy_pm(cx9, cx10, cx11, cx12);
    }

    // channel mix
    color = vec3(color.r*mx86,color.g*mx85,color.b*mx84);

    // strobe
    color += step(0.75,vec3(cnt(50))) * mx83*.9;

    // gamma correction
    color = pow( color, vec3(1.0/mx74) );

    // permanent strobe
    color += step(0.98,cnt(100)) * mx82;

    PixelColor = vec4(color, 1.0) * mx81;
}