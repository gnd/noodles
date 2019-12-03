#version 330
uniform vec2 resolution;
uniform float time;
uniform sampler2D backbuffer;
out vec4 PixelColor;

/// PARAMS
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;
uniform float m5;
uniform float m6;
uniform float m7;

/// COLORS
vec3 cx;
vec3 ck = vec3(0.);
vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cg = vec3(0.,1.,0.);
vec3 cb = vec3(0.,0.,1.);
//////////////////// END OF GLOBALS //////////////

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
#define MAXSTEPS 84
#define MAXDIST 10.0
#define eps .001
#define SPHERES 10
#define FLY

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material mat1;
shading s1;

vec2 p;
vec3 p3;
float asp, spec;
vec3 sky, color; // 'sky color'
vec4 spheres[SPHERES];
float speed = pow(1.01,100.);

// set the camera trajectory
// todo precompute stuff
// make this controllable by midimix
#ifdef FLY
	vec3 ro = vec3(cos(time/3.)*4., sin(time/3.)*3., time*speed);
	vec3 lookat = vec3(cos(time+2./3.)*4.,0.,time*speed + 7.);
	vec3 light_position = vec3(cos(time)*4.,0.,time*speed + sin(time)*4);
#else
	vec3 ro = vec3(cos(time/4.)*3., 2., sin(time/4.)*3.); // ray origin, here also known as 'eye'
	vec3 lookat = vec3(0.,0.,0.);
	vec3 light_position = vec3(sin(time/4.)*4.,2.5,cos(time/4.)*3.);
#endif


// Counter 0 to 1
float cnt(int m) {
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

float smin( float a, float b, float k ) {
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

float blend(float d1, float d2, float k) {
    return smin(d1, d2, k);
}

float sphere(vec3 p, float r) {
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
    return box(p-vec3(.0,-1.*.5,.0), vec3(15.,.01,15.), .0);
}

// taken from https://www.shadertoy.com/view/4ttGDH
float sinfield (vec3 p) {
    //p.xy -= camPath(p.z).xy;
    p.xy -= ro.xy;
	p = cos(p*.315*1.25 + sin(p.zxy*.875*1.25*mx24*10.));
    float n = length(p);
    return (n - 1.025)*1.33;
}

// from https://www.shadertoy.com/view/4t2GDG
float blob (vec3 p) {
    float res = 10000.;
    for( int i = 0; i < SPHERES; i++ ) {
        res = blend(res, sphere(p-spheres[i].xyz, spheres[i].w), 8.);
   }
   return res;
}

float test(vec3 p) {
	return box(p-vec3(0.,0.,time*speed+4), vec3(5.5, 3., 0.1), 0.);
}

float scene(vec3 p) {
    return blend(blob(p), sinfield(p), 6.);
	//return test(p);
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

// Stripe
vec3 str(in float xp, in float yp, in float x, in float y, in vec3 c, in vec3 col) {
        //float ys = y * asp;
        float ys=y;
        vec3 cl = ((p.x > xp-x*.51) && (p.x < x*.51+xp) && (p.y > yp-ys*.51) && (p.y < ys*.51+yp)) ? -c+col : vec3(0.0);
        return cl;
}

// Feedback
vec3 feedb_sqr(in float xpos, in float ypos, in float xsiz, in float ysiz, in float bsiz, in vec3 c) {
    vec3 ccc = texture2D(backbuffer, (vec2(p.x,p.y)-0.51)*bsiz*1.5+0.51).xyz;
    c+=str(xpos,ypos,xsiz,ysiz,c,ccc);
	return c;
}

void main() {
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);

    // enviroment
    l1.color = vec3(1.,1.,1.);
	l1.position = light_position;
    //l1.position = vec3(1.2, 2.0, -1.9);
    color = vec3(.0);
    sky = color;

	// feedback - front
    asp = resolution.x / resolution.y;
    p = gl_FragCoord.xy / resolution;
    color += feedb_sqr(mx76, mx75, mx66, mx65, mx74+m1*.1, color)*mx71*1.05*mx73;

    // setup spheres
    vec4 tmp;
    for( int i = 0; i < SPHERES; i++ ) {
        tmp.x = cos(time *.13 * (float( i )+2.));
        tmp.y = sin(time * .075 * (float( i )+4.));
        tmp.z = sin(time * .1 * (float( i )+3.3)) + time*speed+4.;
        tmp.w = .1 * (sin(time * .1  *(float( i) +1.))+2.) + mx31*2.5;
        spheres[i] = tmp;
    }

    // raymarch a scene
    float step = .0;
    for (int i = 0; i < MAXSTEPS; ++i) {
        p3 = ro + (rd * step);
        float d = scene(p3);
        if (d > MAXDIST) {
            break;
        }
        if (d < eps) {
            // precompute stuff
            vec3 n = normal(p3);
            vec3 ld = normalize(l1.position-p3);
            vec3 ed = normalize(ro-p3);
            if (scene(p3) < eps) {
                /*vec3 pp = mod(p3*100., 2.) - .5;
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
				color = texture2D(sony, (vec2(p3.x+5.5,3.-p3.y))*vec2(.09,.16)).xyz;
				*/
            }
            if (scene(p3) < eps) {
				// needs eps >= .001
                float nv = dot(n, -rd);
                vec3 col = vec3(0.);
                col += sin(nv * vec3(0.0, 1.0, 0.0) * 10.0 * 1.5) * 0.5 + 0.5;
                col += sin(nv * vec3(1.0, 0.0, 0.0) * 20.0 * 1.5) * 0.5 + 0.5;
                col += sin(nv * vec3(0.0, 0.0, 1.0) * 5.0 * 1.5) * 0.5 + 0.5;
                col = 1.1 - col;
                mat1.color = clamp(normalize(col), 0.0, 1.0);
                spec = clamp(3.5 * pow(dot(normalize(reflect(-ld, n)), ed), 70.), 0., 1.);;
                color += mat1.color * clamp(0.4*n.y+0.5, 0., 1.) * 3.;
            }
            break;
        }
        step += d;
    }

    // Exponential distance fog
    color = mix(color, 0.8 * sky, 1.0 - exp2(-0.010 * step * step));

    // gamma correction
    color = mix(color, pow(color, vec3(1.0/2.2)), mx82);

	// feedback - back
	asp = resolution.x / resolution.y;
	p = gl_FragCoord.xy / resolution;
	color += feedb_sqr(mx76, mx75, mx66, mx65, mx74+m1*.1, color)*mx71*1.05*mx72;

    // channel mix
	color = vec3(color.r*mx86,color.g*mx85,color.b*mx84);

    // strobe
    color += vec3(cnt(50))*mx83*.9;

    // send to screen
    PixelColor = vec4(color+spec, 1.)*mx81;
}
