//////////////////// GLOBALS /////////////////
uniform float time;
uniform vec2 resolution;
uniform vec2 rand;
uniform vec2 mouse;
uniform sampler2D backbuffer;
vec2 p; //relative position
vec3 c; //color
float mw = 0.001; //minimal width
float asp; //aspect ratio

///////// PARAMS
uniform float m0;
uniform float m1;
uniform float m2;
uniform float m3;
uniform float m4;
uniform float m5;
uniform float m6;
uniform float m7;
uniform float m8;
uniform float m9;
uniform float m11;
uniform float m12;
uniform float m14;
uniform float m15;
uniform float m20;

float mx11 = m8 / 127.;
float mx12 = m9 / 127.;
float mx14 = m11 / 127.;
float mx15 = m12 / 127.;
float mx21 = m14 / 127.;
float mx22 = m15 / 127.;
float mx31 = m20 / 127.;

// counter 0 to 1 
float cnt(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

// horizontal line
vec3 liney(in float yp, in float width, in vec3 col)
{
        return ((p.y <= yp+width/2.0) && (p.y >= yp-width/2.0)) ? col : vec3(0.0);
	
}

// square
vec3 sqr(in float xp, in float yp, in float s, in vec3 col)
{
	float ys = s * asp;
	vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
	return cl;
}

// 2D rotation
vec2 rot(in float angle)
{
	float dst = length(vec2(p.x,p.y/asp));
	float gamma = (p.y > 0.) ? degrees(acos(p.x / dst)) : 360. - degrees(acos(p.x / dst));
	float alfa = mod(angle*360. + gamma, 360.);
	return dst * vec2(cos(radians(alfa))/asp, sin(radians(alfa))) * asp;
}


// glslsandbox
float Hash( vec2 p)
{
	vec3 p2 = vec3(p.xy,0.5);
	return fract(sin(dot(p2,vec3(372.1,621.7, 200.4)))*378.);
}

float noise(in vec2 p)
{
	vec2 i = floor(p);
	vec2 f = fract(p);
	f *= f * (3.0-2.0*f);
	return mix(mix(Hash(i + vec2(0.,0.)), Hash(i + vec2(1.,0.)),f.x),
               mix(Hash(i + vec2(0.,1.)), Hash(i + vec2(1.,1.)),f.x),
               f.y);
}

float fbm(vec2 p)
{
	float v = 0.0;
	v += noise(p*0.5)*1.5;
	v += noise(p*2.)*-1.25;
	v += noise(p*4.)*-1.125;
	v += noise(p*10.)*.125;

	return v * 91.0;
}

vec4 mod289(vec4 x)
{
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}
 
vec4 permute(vec4 x)
{
    return mod289(((x*34.0)+1.0)*x);
}
 
vec4 taylorInvSqrt(vec4 r)
{
    return 1.79284291400159 - 0.85373472095314 * r;
}
 
vec2 fade(vec2 t) {
    return t*t*t*(t*(t*6.0-15.0)+10.0);
}

// Classic Perlin noise
float cnoise(vec2 P)
{
    vec4 Pi = floor(P.xyxy) + vec4(0.0, 0.0, 1.0, 1.0);
    vec4 Pf = fract(P.xyxy) - vec4(0.0, 0.0, 1.0, 1.0);
    Pi = mod289(Pi); // To avoid truncation effects in permutation
    vec4 ix = Pi.xzxz;
    vec4 iy = Pi.yyww;
    vec4 fx = Pf.xzxz;
    vec4 fy = Pf.yyww;
     
    vec4 i = permute(permute(ix) + iy);
     
    vec4 gx = fract(i * (1.0 / 41.0)) * 2.0 - 1.0 ;
    vec4 gy = abs(gx) - 0.5 ;
    vec4 tx = floor(gx + 0.5);
    gx = gx - tx;
     
    vec2 g00 = vec2(gx.x,gy.x);
    vec2 g10 = vec2(gx.y,gy.y);
    vec2 g01 = vec2(gx.z,gy.z);
    vec2 g11 = vec2(gx.w,gy.w);
     
    vec4 norm = taylorInvSqrt(vec4(dot(g00, g00), dot(g01, g01), dot(g10, g10), dot(g11, g11)));
    g00 *= norm.x;  
    g01 *= norm.y;  
    g10 *= norm.z;  
    g11 *= norm.w;  
     
    float n00 = dot(g00, vec2(fx.x, fy.x));
    float n10 = dot(g10, vec2(fx.y, fy.y));
    float n01 = dot(g01, vec2(fx.z, fy.z));
    float n11 = dot(g11, vec2(fx.w, fy.w));
     
    vec2 fade_xy = fade(Pf.xy);
    vec2 n_x = mix(vec2(n00, n01), vec2(n10, n11), fade_xy.x);
    float n_xy = mix(n_x.x, n_x.y, fade_xy.y);
    return 2.3 * n_xy;
}

// Fractional Brownian Motion 
float fbm2(vec2 p, int octaves, float lacunarity, float gain)
{
    float sum = 0.0;
    float amp = 1.0;
    vec2 pp = p;
     
    int i;
     
    for(i = 0; i < octaves; i+=1)
    {
        amp *= gain; 
        sum += amp * cnoise(pp);
        pp *= lacunarity;
    }
    return sum;
}


// FLAKY
vec3 flacks(in int iter, in vec3 col) {
	vec3 fc = vec3(.0);
	for( int i=1; i < iter; ++i ) { 
	        float t = abs(((0.1*time)*0.1) /  ((p.x + fbm( p+ time/float(i)))*50.));
	        fc +=  t * col;
	}
	return fc; 
}


// FLAKY
vec3 flacks2(in int iter, in vec3 col) {
        vec3 fc = vec3(.0);                                                                                                                                                                                  
        for( int i=1; i < iter; ++i ) {                                                                                                                                                                      
                float t = abs(((0.1*time)*0.1) /  ((p.x + fbm2( p+ time/float(i), 1, mx21, mx31*10.)    )*50.));
                fc +=  t * col;                                                                                                                                                                              
        }                                                                                                                                                                                                    
        return fc*.051;
}      


void main( void ) {
////////////////////// INIT ////////////////////////
p = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
p = vec2(p.y, p.x);
asp = (resolution.y / resolution.x); ///1.3333;	//	1.333 is a fix for 16:9 monitors at 1024x768 monitors (lunchmeat related)

vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cg = vec3(0.,0.,1.);
vec3 cb = vec3(0.,1.,0.);

////////////////////// END OF INIT /////////////////
c=vec3(0.);

float cen = .5;
p -= vec2(cen);
p = rot(m0);
p += vec2(cen);
//c += sqr(0.5, 0.5, 0.1, cw);

// CHECK BACKBUFFER
if (mx22 == 1.) {
vec2 pp = vec2(p.y,1.-p.x); 		// why reverse the position ?
c += texture2D(backbuffer, (pp-0.5)*(1.+0.1)+0.5).xyz*mx14;
} 

// RGB LINES
//c+=liney(0.2,0.1,cr);
//c+=liney(0.4,0.1,cb);
//c+=liney(0.6,0.1,cg);

p*=mx15*20.;
// FLAKY
if (mx12 == 1.) {
c+=flacks(15, vec3(cnt(100)));
} else {
p-=0.5;
c+=flacks2(15, vec3(cnt(100)));
}

// DRAW
gl_FragColor = vec4(c, 1.0);

}
