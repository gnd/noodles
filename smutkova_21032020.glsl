#version 330
out vec3 PixelColor;
uniform vec2 rand;
uniform vec2 resolution;
uniform sampler2D sony;
uniform sampler2D backbuffer;
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
float mx11 = m8/25.4;
float mx12 = m9/127.;
float mx13 = m10/127.;
float mx14 = m11/127.;
float mx15 = m12/127.;
float mx16 = m13/127.;
float mx21 = m14/63.5;
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
vec3 p;
vec3 color;
float asp = 1.;
//float vmax(in vec3 p ) {
//    return max(p.x,max(p.y,p.z));
//}

float scene(vec3 p, float ty) {
    // initial cube
    //float d = vmax(abs(p) - mx31);
    float d = max(abs(p.x),max(abs(p.y),abs(p.z))) - 1.;
    float s = 1.;

    for (int m=0; m<5; m++) {
         vec3 a = mod(p*s, mx11)-(mx11*.5);
         vec3 r = 3.*-abs(a) + mx21+m0*mx51;
         s *= 3.0;
         float aa = max(r.y, r.z);
         float ab = max(r.x, r.z);
         float ac = max(r.x, r.y);
         float c = (min(aa, min(ab,ac))-.22)/s;
         d = max(d,-c);
    }
    //return (1.-ty)*.01 + d;
    return ty * -.01 + .01 + d;
}

// taken from http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal(vec3 p, float ty) {
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*scene( k.xyy*eps + p, ty) +
                      k.yyx*scene( k.yyx*eps + p, ty) +
                      k.yxy*scene( k.yxy*eps + p, ty));// +
                      //k.xxx*scene( p + k.xxx*eps , ty) );
}

// Random number generator
float rnd(vec2 co){
        // implementation found at: lumina.sourceforge.net/Tutorials/Noise.html
        float k = sin(dot(co.xy ,vec2(12.9898,78.233)));
        return fract(k + k);
}

// Noise2f
float noise2f( in vec2 p )
{
        vec2 ip = vec2(floor(p));
        vec2 u = fract(p);
        // http://www.iquilezles.org/www/articles/morenoise/morenoise.htm
        //u = u*u;
        u = u*u*u*((6.0*u-15.0)*u+10.0);

        float res = mix(
                                        mix(rnd(ip), rnd(ip+vec2(1.0,0.0)), u.x),
                                        mix(rnd(ip+vec2(0.0,1.0)), rnd(ip+vec2(1.0,1.0)), u.x),
                                        u.y);

        return res - 0.25;
        //return 2.0* (res-0.5);
}

// Stripe
vec3 str(in float xp, in float yp, in float x, in float y, in vec3 c, in vec3 col)
{
        float ys = y * asp;
        vec3 cl = ((p.x > xp- x/2.) && (p.x < xp + x/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? -c+col : vec3(0.0);
        return cl;
}

vec3 feedb_sqr(in float xpos, in float ypos, in float xsiz, in float ysiz, in float bsiz, in float blnd, in float brg, in vec3 c)
{
	vec2 mm = (p.x > xpos-xsiz/2.) && (p.x < xpos+xsiz/2.) && (p.y > ypos-ysiz/2.) && (p.y < ypos+ysiz/2.) ? (vec2(p.x,p.y)-0.5)*bsiz+0.5 : vec2(0.);
    vec3 ccc = texture2D(backbuffer, mm).xyz;
    c+=str(xpos,ypos,xsiz,ysiz,c,ccc);
	return c;
}


void main() {
    float speed = time/((mx41+m2)*1100.);
    float sintime = sin(speed);
    float costime = cos(speed);
    // original eye
    vec3 eye = vec3(sintime*10., sintime*2.+2., costime*10.);
    // original lookat
    vec3 lookat = vec3(costime*.5, mx61*2.-1., sintime*.5);
    vec3 fwd = normalize(lookat-eye);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x ));
    //vec3 up = normalize(cross(fwd, right));
    vec3 up = vec3(0.,1.,0.);
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    //vec3 ro = eye;
    vec3 rd = normalize(fwd * (mx31*550.) + right * (uv.x*aspect) + up*uv.y);
    // enviroment
    vec3 lpos = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
    color = vec3(.0);

    vec3 pp = p;
    // raymarch a scene
    float step = .0;
    for (int i = 0; i < MAXSTEPS; i+=2) {
        p = rd * step + eye;
        vec3 texxcol = texture(sony, vec2(p.x*.27+.5,p.z*.24+.5)).rgb;
        float d = scene(p, texxcol.y);
        if (d < eps) {
            // precompute stuff
            vec3 n = normal(p, texxcol.y);
            vec3 ld = lpos-p;
            vec3 ed = eye-p;
            color = vec3(clamp(dot(n,ld), 0.,.4)
                    + clamp(dot(reflect(-ld, n), ed),0.,.48)
                    + n.y*.4+.01);
            // add step fog
            color.r += float(i)/float(MAXSTEPS)*2.;
            color /= texxcol;;
            break;
        }
        step += d;
    }

    // Exponential distance fog FIXME
    PixelColor = mix(color, vec3(0.), .6 - exp2(-.9 * step));

    p = gl_FragCoord.xyy / resolution.xyy;
	asp = resolution.x / resolution.y;

    PixelColor += vec3(noise2f(p.xy*rand*100.)*mx91*m2,0.,0.);
    PixelColor += feedb_sqr(mx86, mx85, mx76, mx75, mx81, mx71, mx61, PixelColor)*mx71;
}
