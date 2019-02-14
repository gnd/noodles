/////////////////////////////////////// PRIMITIVES START
// to check for errors directly
//uniform float time;
//vec2 p = vec2(0.);
//float asp = 1.;
//uniform sampler2D backbuffer;
//uniform vec2 rand;

/////////////////////////////////////// SIMPLE FUNCTIONS

// Sin between 0 - 1
float zsin(in float a) {
	return (sin(a) + 1.0) / 2.0;
}

// Sin between 0 - 1 - bouncing
float a_sin(in float a) {
	return abs(sin(a));
}

// Cos between 0 - 1
float zcos(in float a) {
    return (cos(a) + 1.0) / 2.0;
}

// Cos between 0 -1 - bouncing
float a_cos(in float a) {
    return abs(cos(a));
}

// Modulo
float mdl(in float a, in float b)
{
        return a - b * floor(a/b + 0.001);
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

// Counter 0 to 1
float cnt (in int ms)
{
	float p = float(ms) / 1000.;	// we need a fraction of 1000ms
	return mod(time,p) / p; 	// we return normalized modulo of 1000ms
}


// Old Counter 0 to 1 - anyway WTF
float ocnt(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

/////////////////////////////////////////// NOISE

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
float fbm(vec2 p, int octaves, float lacunarity, float gain)
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

// simpler noise  + fbm (from glslsandbox)
float hash( vec2 p)
{
	vec3 p2 = vec3(p.xy,0.5);
	return fract(sin(dot(p2,vec3(372.1,621.7, 200.4)))*378.);
}

float noise(in vec2 p)
{
	vec2 i = floor(p);
	vec2 f = fract(p);
	f *= f * (3.0-2.0*f);
	return mix(mix(hash(i + vec2(0.,0.)), hash(i + vec2(1.,0.)),f.x),
		mix(hash(i + vec2(0.,1.)), hash(i + vec2(1.,1.)),f.x),
		f.y);
}

// simpler FBM (seems its also faster)
float fbm2(vec2 p)
{
	float v = 0.0;
	v += noise(p*0.5)*1.5;
	v += noise(p*2.)*-1.25;
	v += noise(p*4.)*-1.125;
	v += noise(p*10.)*.125;

	return v * 91.0;
}

float noise_ta(in float amp)
{
	amp += .199;
	float noise = 1.-noise2f(vec2(sin(cnoise(p)*time)*amp*10.));
        return pow(noise,5.) * step(0.2,amp);
}

float noise_pow(in float amp)
{
        vec2 _p = p;
        float tt = 30.+sin(time)*10.*amp;
        float cx = cos(time)/100.;
        _p+= noise2f(_p*tt*amp);
        cx+= noise2f(_p*tt);
        float cxx=cx*cx*10.4;
        return 1.-cxx;
}



//domain distort - perlin x
vec2 dx_x(in float amp)
{
	return vec2(p.x + cnoise(vec2(p.y*60., rand.x))*amp, p.y);
}


// domain distort - perlin #1
vec2 dx_p1(in float mod1, in int iter, in float gain)
{
        vec2 q = vec2( fbm( p + vec2(0.0,mod1), iter, 2., gain ),
                       fbm( p + vec2(5.2,1.3), iter, 2., gain ) );

        return (p+4.*q);
}


// domain distort - perlin #2
vec2 dx_p2(in float mod1, in int iter, in float gain)
{
        vec2 q = vec2( fbm( p + vec2(0.0,mod1), iter, 2., gain ),
                       fbm( p + vec2(5.2,1.3), iter, 2., gain ) );

        vec2 r = vec2( fbm( p + 4.0*q + vec2(1.7,9.2), iter, 2., gain ),
                       fbm( p + 4.0*q + vec2(8.3,2.8), iter, 2., gain ) );
	return (p+4.*r);
}


/////////////////////////////////////////// MACHINELETS

float strobe(in float freq)
{
	float fr = freq / 1000.; 					// normalize
	float cutoff = (fr > 1.) ? 0.98 : log(10.*fr)*0.1+0.8;		// modify cutoff according to strobe freq, devised via http://www.intmath.com/functions-and-graphs/graphs-using-jsxgraph.php
	return step(cutoff, cnt(int(freq)));
}


/////////////////////////////////////////// GEOMETRIC MODIFIERS

// 2D rotation
vec2 rot2d(in float angle)
{
	float dst = length(vec2(p.x,p.y/asp));
	float gamma = (p.y > 0.) ? degrees(acos(p.x / dst)) : 360. - degrees(acos(p.x / dst));
	float alfa = mod(angle*360. + gamma, 360.);
	return dst * vec2(cos(radians(alfa))/asp, sin(radians(alfa))) * asp;
}


mat4 inverse(mat4 m) {
  float
      a00 = m[0][0], a01 = m[0][1], a02 = m[0][2], a03 = m[0][3],
      a10 = m[1][0], a11 = m[1][1], a12 = m[1][2], a13 = m[1][3],
      a20 = m[2][0], a21 = m[2][1], a22 = m[2][2], a23 = m[2][3],
      a30 = m[3][0], a31 = m[3][1], a32 = m[3][2], a33 = m[3][3],

      b00 = a00 * a11 - a01 * a10,
      b01 = a00 * a12 - a02 * a10,
      b02 = a00 * a13 - a03 * a10,
      b03 = a01 * a12 - a02 * a11,
      b04 = a01 * a13 - a03 * a11,
      b05 = a02 * a13 - a03 * a12,
      b06 = a20 * a31 - a21 * a30,
      b07 = a20 * a32 - a22 * a30,
      b08 = a20 * a33 - a23 * a30,
      b09 = a21 * a32 - a22 * a31,
      b10 = a21 * a33 - a23 * a31,
      b11 = a22 * a33 - a23 * a32,

      det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

  return mat4(
      a11 * b11 - a12 * b10 + a13 * b09,
      a02 * b10 - a01 * b11 - a03 * b09,
      a31 * b05 - a32 * b04 + a33 * b03,
      a22 * b04 - a21 * b05 - a23 * b03,
      a12 * b08 - a10 * b11 - a13 * b07,
      a00 * b11 - a02 * b08 + a03 * b07,
      a32 * b02 - a30 * b05 - a33 * b01,
      a20 * b05 - a22 * b02 + a23 * b01,
      a10 * b10 - a11 * b08 + a13 * b06,
      a01 * b08 - a00 * b10 - a03 * b06,
      a30 * b04 - a31 * b02 + a33 * b00,
      a21 * b02 - a20 * b04 - a23 * b00,
      a11 * b07 - a10 * b09 - a12 * b06,
      a00 * b09 - a01 * b07 + a02 * b06,
      a31 * b01 - a30 * b03 - a32 * b00,
      a20 * b03 - a21 * b01 + a22 * b00) / det;
}


mat4 rmat(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat4(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s, oc * axis.z * axis.x + axis.y * s, 0.0,
                oc * axis.x * axis.y + axis.z * s, oc * axis.y * axis.y + c, oc * axis.y * axis.z - axis.x * s, 0.0,
                oc * axis.z * axis.x - axis.y * s, oc * axis.y * axis.z + axis.x * s, oc * axis.z * axis.z + c, 0.0,
                0.0, 0.0, 0.0, 1.0);
}

vec3 rot( vec3 p, mat4 m )
{
    vec4 q = inverse(m)*vec4(p,1.);
    return q.xyz;
}

float opS( float d1, float d2 )
{
    return max(-d2,d1);
}

vec2 opU( vec2 d1, vec2 d2 )
{
	return (d1.x<d2.x) ? d1 : d2;
}

// smooth min from: https://github.com/stackgl/glsl-smooth-min/blob/master/poly.glsl
float smin(float a, float b, float k) {
	float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
	return mix(b, a, h) - k * h * (1.0 - h);
}

///////////////////////////////////////////// RAYMARCHING PRIMITIVES

// Sphere field
float sphere(in vec3 p, in vec3 centerPos, float radius) {
        return length(p-centerPos) - radius;
}

// Box field
float box( vec3 p, vec3 b ) {
	return length(max(abs(p)-b,.5));
}

float hex( vec3 p, vec2 h )
{
	vec3 q = abs(p);
	return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
}

// Sinusoid bumps
float sinusoidBumps(in vec3 p, float distort){
        float frq = 10.;
        return distort * 3.*sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31);
}

////////////////////////////////////// SIMPLE PRIMITIVES

// Vertical line
vec3 linx(in float xp, in float width, in vec3 col)
{
        return ((p.x <= xp+width/2.) && (p.x >= xp-width/2.)) ? col : vec3(0.0);
}

// Horizontal line
vec3 liny(in float yp, in float width, in vec3 col)
{
        return ((p.y <= yp+width/2.) && (p.y >= yp-width/2.)) ? col : vec3(0.0);

}

// Square
vec3 sqr(in float xp, in float yp, in float s, in vec3 col)
{
        float ys = s * asp;
        vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
        return cl;
}

// Square - remove background
vec3 sqr_del(in float xp, in float yp, in float s, in vec3 c, in vec3 col)
{
        float ys = s * asp;
        vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? -c+col : vec3(0.0);
        return cl;
}

// Square ring
vec3 sqrr(in float xp, in float yp, in float s, in float w, in vec3 col)
{
        float r = s / 2.;
        float rs = s * asp / 2.;
        w /= 2.;
        vec3 cl = ((((p.x > xp-r-w) && (p.x < xp-r+ w)) || ((p.x < xp+r+w) && (p.x > xp+r-w))) && (p.y > yp-rs-w) && (p.y < yp+rs+w)) ||
                  ((((p.y > yp-rs-w*asp) && (p.y < yp-rs+w*asp)) || ((p.y < yp+rs+w*asp) && (p.y > yp+rs-w*asp))) && (p.x > xp-r-w) && (p.x < xp+r+w))
                  ? col : vec3(0.0);
        return cl;
}

// Rectangle
vec3 rect(in float xp, in float yp, in float sx, in float sy, in vec3 col)
{
        float syy = sy*asp;
        vec3 cl = ((p.x > xp - sx/2.) && (p.x < xp + sx/2.) && (p.y > yp-syy/2.) && (p.y < yp + syy/2.)) ? col : vec3(0.0);
        return cl;
}

// Rectangle by corner
vec3 rect_c(in float xp, in float yp, in float sx, in float sy, in vec3 col)
{
        float syy = sy;
        vec3 cl = ((p.x > xp) && (p.x < xp + sx) && (p.y > yp) && (p.y < yp + syy)) ? col : vec3(0.0);
        return cl;
}

// Stripe
vec3 str(in float xp, in float yp, in float x, in float y, in vec3 c, in vec3 col)
{
        float ys = y * asp;
        vec3 cl = ((p.x > xp- x/2.) && (p.x < xp + x/2.) && (p.y > yp-ys/2.) && (p.y < yp + ys/2.)) ? -c+col : vec3(0.0);
        return cl;
}

// Circle
vec3 circ(in float xp, in float yp, in float r, in vec3 col)
{
        // correct aspect ratio transfer
        xp = xp * asp;
        float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
        vec3 cl = dist < r ? col : vec3(0.);
        return cl;
}

// Circle - remove background
vec3 circ_del(in float xp, in float yp, in float r, in vec3 c, in vec3 col)
{
        // correct aspect ratio transfer
        xp = xp * asp;
        float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
        vec3 cl = dist < r ? -c+col : vec3(0.);
        return cl;
}

// Circle ring
vec3 circr(in float xp, in float yp, in float r, in float w, in vec3 col)
{
        // correct aspect ratio transfer
        xp = xp * asp;
        float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
        vec3 cl = (((dist < r) && (dist > r - w/2.)) || ((dist > r) && (dist < r + w/2.)))  ? col : vec3(0.);
        return cl;
}

// Oscillator on the x-axis
vec3 oscx(in float xp, in float w, in float freq, in float amp, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.x - xp) + sin(p.y * freq * 100.)*amp);
        return cl * c;
}

// Oscillator on the y-axis
vec3 oscy(in float yp, in float w, in float freq, in float amp, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.y - yp) + cos(p.x * freq * 100.)*amp);
        return cl * c;
}

// vertical Z-lines repeated into sin-bordel
float zrpty(in float xp, in float w, in float rep, in float amp, in float freq, in float mod_amp)
{
        w /= 10.;
        rep = 1. - rep*.999;
        mod_amp += .1;
        float rev_xp = 1. - xp;
        float sinmod = sin(p.y*freq);
        float diagonal_lines = w / abs((p.x-xp) + mod(p.y,rep) + sinmod*amp);
        float vertical_mask = step(0.99,rep / (1.-p.x-rev_xp+sinmod*-amp));
        float horizontal_lines = w / abs(mod(p.y,rep));

        return mod ( diagonal_lines + horizontal_lines*vertical_mask, 1./mod_amp ) * mod_amp;
}

// horizontal Z-lines repeated into sin-bordel
float zrptx(in float yp, in float w, in float rep, in float amp, in float freq, in float mod_amp)
{
        w /= 10.;
        rep = 1. - rep*.999;
        mod_amp += .1;
        float rev_yp = 1. - yp;
        float sinmod = sin(p.x*freq);
        float diagonal_lines = w / abs((p.y-yp) + mod(p.x,rep) + sinmod*amp);
        float horizontal_mask = step(0.99,rep / (1.0-p.y-rev_yp+sinmod*-amp));
        float vertical_lines = w / abs(mod(p.x,rep));

        return mod ( diagonal_lines + vertical_lines*horizontal_mask, 1./mod_amp ) * mod_amp;
}


// basically iterative fbm
vec3 flacks(in int iter, in vec3 col, in float distort ) {
        vec3 fc = vec3(.0);
        for( int i=1; i < iter; ++i ) {
                float t = abs(((0.1*time)*0.1) /  ((p.x + fbm( p+ time/float(i), 1, 0.5, distort)    )*50.));
                fc +=  t * col;
        }

        // this can be adjusted for +/- brightness
        return fc;
}


// basically iterative fbm
vec3 flacks2(in int iter, in vec3 col, in float distort ) {
        vec3 fc = vec3(.0);
        for( int i=1; i < iter; ++i ) {
                float t = abs(((0.1*time)*0.1) /  ((p.x + fbm2( p+ time/float(i) )    )*50.));
                fc +=  t * col;
        }

        // this can be adjusted for +/- brightness
        return fc*.051;
}

/////////////////////////////////////////// TEXTURES AND PAINTING
vec3 feedb_sqr(in float xpos, in float ypos, in float xsiz, in float ysiz, in float bsiz, in float blnd, in float brg, in vec3 c)
{
        vec2 mm = (p.x > xpos-xsiz/2.) && (p.x < xpos+xsiz/2.) && (p.y > ypos-ysiz/2.) && (p.y < ypos+ysiz/2.) ? (vec2(p.x,p.y)-0.5)*bsiz+0.5 : vec2(0.);
        vec3 ccc = texture2D(backbuffer, mm).xyz;
        c+=str(xpos,ypos,xsiz,ysiz,c,ccc);
        return c;
}

vec3 feedb_crc(in float xpos, in float ypos, in float siz, in float bsiz, in vec3 c)
{
      float dist = sqrt( (p.x - xpos) * (p.x - xpos) + (p.y - ypos) * (p.y - ypos) );
      vec2 mm = dist < xpos*asp ? (vec2(p.x,p.y)-0.5)*bsiz+0.5 : vec2(0.);
      vec3 ccc = texture2D(backbuffer, mm).xyz;
      c+=circ_del(xpos,ypos,siz,c,ccc);
      return c;
}




////////////////////////////////////////////// PRIMITIVES END

// to check for errors directl
//void main(void) {
//gl_FragColor = vec4(vec3(1.), 1.0);
//}
