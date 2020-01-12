#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform sampler2D sony;
uniform float time;
#define MAXSTEPS 256
#define MAXDIST 10.0
#define eps 0.00001

/// PARAMS
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;
uniform float m5;
uniform float m6;
uniform float m7;

/////////////////// MIDIMIX //////////////////////
uniform float m8,m9,m10;
uniform float m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30;
uniform float m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50;
uniform float m51,m52,m53,m54,m55,m56,m57,m58,m59;
float mx11 = m8/127.;
float mx12 = m9/127.;
float mx13 = m10/127.;
float mx14 = m11/127.;
float mx15 = m12/127.;
float mx16 = m13/127.;
float mx21 = m14/127.;
float mx22 = m15/127.;
float mx23 = m16/127.;
float mx24 = m17/127.;
float mx25 = m18/127.;
float mx26 = m19/127.;
float mx31 = m20/127.;
float mx32 = m21/127.;
float mx33 = m22/127.;
float mx34 = m23/127.;
float mx35 = m24/127.;
float mx36 = m25/127.;
float mx41 = m26/127.;
float mx42 = m27/127.;
float mx43 = m28/127.;
float mx44 = m29/127.;
float mx45 = m30/127.;
float mx46 = m31/127.;
float mx51 = m32/127.;
float mx52 = m33/127.;
float mx53 = m34/127.;
float mx54 = m35/127.;
float mx55 = m36/127.;
float mx56 = m37/127.;
float mx61 = m38/127.;
float mx62 = m39/127.;
float mx63 = m40/127.;
float mx64 = m41/127.;
float mx65 = m42/127.;
float mx66 = m43/127.;
float mx71 = m44/127.;
float mx72 = m45/127.;
float mx73 = m46/127.;
float mx74 = m47/127.;
float mx75 = m48/127.;
float mx76 = m49/127.;
float mx81 = m50/127.;
float mx82 = m51/127.;
float mx83 = m52/127.;
float mx84 = m53/127.;
float mx85 = m54/127.;
float mx86 = m55/127.;
float mx91 = m56/127.;
float mx92 = m57/127.;
float mx93 = m58/127.;
float mx94 = m59/127.;
///////////////// END OF MIDIMIX ////////////////////

struct light { vec3 position; vec3 color; };
light l1;
vec3 color;

float vmax(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

// sony height mapping
vec3 texcol(vec3 p) {
    return texture2D(sony, vec2(p.x*.27+.5,p.z*.24+.5)).rgb;
}

float scene(vec3 p) {
    // initial cube
    float d = vmax(abs(p) - mx31*2.);
    float s = 1.;

    for (int m=0; m<5; m++) {
         vec3 a = mod(p*s, mx11*5.)-(mx11*2.5);
         s *= 3.0;
         vec3 r = mx21*2. - 3.*abs(a);
         float aa = max(r.y, r.z);
         float ab = max(r.x, r.z);
         float ac = max(r.x, r.y);
         float c = (min(aa, min(ab,ac))-.22)/s;
         d = max(d,-c);
    }
    return (1.-texcol(p).y)*.01 + d;
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p) {
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*scene( p + k.xyy*eps ) +
                      k.yyx*scene( p + k.yyx*eps ) +
                      k.yxy*scene( p + k.yxy*eps ) +
                      k.xxx*scene( p + k.xxx*eps ) );
}

void main() {
    float speed = mx51*1100.;
    float sintime = sin(time/speed);
    float costime = cos(time/speed);
    // original eye
    vec3 eye = vec3(sintime*10., sintime*2.+2., costime*10.);
    //eye = vec3(0,10.,0.);

    // original lookat
    vec3 lookat = vec3(costime*.5, 2.*mx61, sintime*.5);

    vec3 fwd = normalize(lookat-eye);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x ));
    vec3 up = normalize(cross(fwd, right));
    //float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    //float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 ro = eye;
    vec3 rd = normalize(fwd * (mx41*150.) + right * (uv.x*aspect) + up * uv.y);

    // enviroment
    l1.color = vec3(1.,1.,1.);
    l1.position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    //l1.position = vec3(1.2, 2.0, -1.9);
    color = vec3(.0);

    // raymarch a scene
    float step = .0;
    vec3 p, p_refr;
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
            vec3 col = vec3(.4);
            // no need to distinguish objects here
            color = col * clamp(dot(n,ld), 0., 1.);
            color += col * clamp(dot(normalize(reflect(-ld, n)), ed), 0., 1.2);
            color += col * clamp(.4*n.y+.5, 0., .01);
            // add step fog
            color += vec3(float(i)/float(MAXSTEPS),0.,0.)*2.;
            color /= texcol(p) * .9;
            break;
        }
        step += d;
    }

    // Exponential distance fog
    color = mix(color, 0.8 * vec3(0.), 1.0 - exp2(-.003 * step * step * step));

    // gamma correction
    color = pow( color, vec3(1.0/2.2) );

    // send to screen
    PixelColor = vec4(color, 1.);
}
