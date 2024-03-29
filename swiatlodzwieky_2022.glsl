#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
uniform sampler2D backbuffer;
vec2 p2;

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


#define MAXSTEPS 256
#define MAXDIST 12
#define eps 0.00001

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material mt1,mt2,mt3;
shading s1,s2,s3;

vec3 sky, color; // 'sky color'
vec3 cr = vec3(1.,0.,0.);

// Sinusoid bumps
float sinbumps(in vec3 p, in float frq) {
    return sin(p.x*frq) * sin(p.y*frq+time*4.) * sin(p.z*frq) * .17;
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

float fract2(vec3 p) {
    
/*
    float kk = 10.0*mx21; // or some other amount
    float cc = cos(kk*p.y);
    float ss = sin(kk*p.y);
    mat2  mm = mat2(cc,-ss,ss,cc);
    p = vec3(mm*p.xz,p.y); */

    // i = 0
    //p = vec3(R(p.xy, mx11*2.), p.z);
    //p = vec3(p.x, R(p.yz, mx11*2.));
    vec2 xz = R(p.xz, mx11*2.);
    p = vec3(xz.x, p.y, xz.y);
    float d = sdCross_smooth(p);

    // i = 1
    float s = 3.0;
    float c = 0.5;
    vec3 l = vec3(4.,1.,4.);
    p = p-c*clamp(round(vec3(R(p.xy, mx11*2.), p.z)/c),-l,l);
    d = min(d, sdCross_smooth(p*s)/s);

    // i = 2
    s *= 3;
    c = 1.;
    l = vec3(4.,1.,4.);
    p = p-c*clamp(round(p/c),-l,l);
    d = min(d, sdCross_smooth(p*s)/s);

    float disp = mx31*10.;
    float pd = sin(disp*p.x)*sin(disp*p.y)*sin(disp*p.z);
    
    return d+pd;
}

float fract(vec3 p) {
    float d = box(p);
    float s = 1.;
    for( int m=0; m<4; m++ ) {
        vec3 a = mod( p*s, 2.0 )-1.0;
        s *= 3.0;
        vec3 r = 1.0 - 3.0*abs(a);

        float c = sdCross(r)/s;
        d = max(d,c);
    }
    return d;
}

float scene(vec3 p) {
    return fract2(p);

    //min (
      //          fract(p),
          //      plane(p) 
        //        );
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
    vec3 ro = vec3(3., 2.5, 3.);
    ro = vec3(cos(time/4.)*2., 1., sin(time/4.)*2.);
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
            
            if (fract(p) < eps) {
                mt1.color = vec3(0.,1.,0.);
                mt1.reflection_ratio = 5.;
                mt1.shininess = 8.;
                s1 = get_shading(mt1, l1, p, n, ld, ed);
                color = mt1.color * s1.diffuse * s1.shadow/2.;
                color += s1.specular;
                color += mt1.color * s1.amb *.2;
                if (n.x > 0.5) {
                    color+=vec3(1.,0.,p.z);
                } else {
                    color+=vec3(0.,p.y,1.);
                }
            }

            if (fract2(p) < eps) {
                mt1.color = vec3(0.,1.,0.);
                mt1.reflection_ratio = 5.;
                mt1.shininess = 8.;
                s1 = get_shading(mt1, l1, p, n, ld, ed);
                color = mt1.color * s1.diffuse * s1.shadow/2.;
                color += s1.specular;
                color += mt1.color * s1.amb *.2;
                if (n.x > 0.5) {
                    color+=vec3(1.,0.,p.z);
                } else {
                    color+=vec3(0.,p.y,1.);
                }
            }

                    
            /*if (plane(p) < eps) {
                m1.color = vec3(mod(floor(p.x)+floor(p.z)-.5, 2.0));
                m1.reflection_ratio = 0.01;
                m1.shininess = 3.;
                s1 = get_shading(m1, l1, p, n, ld, ed);
                color = m1.color * s1.diffuse * s1.shadow/2.;
                color += s1.specular;
                color += m1.color * s1.amb *.2;
            }*/
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

    // set coords for feedback
    float asp = resolution.x / resolution.y;
    p2 = gl_FragCoord.xy / resolution;

    // feedback normal
//    if (mx72 > 0.) {
        vec3 col = feedb_sqr(.5, .5, 1., 1., mx74, PixelColor.xyz)*mx71*1.1;
        PixelColor += vec4(col, 1.);
  //  }
}
