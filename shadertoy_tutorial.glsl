uniform vec2 resolution;
uniform float time;
/*
https://www.shadertoy.com/view/XllGW4
HOWTO Get Started With Ray Marching
by Michael Pohoreski aka MysticReddit
Version 0.51, July 2017
*/

#define AUTO_ROTATE     // uncomment to stop auto camera rotation
//#define BACKGROUND_BLUE // uncomment for blue background, else cubemap background
//#define VIEW_ISOMETRIC  // Nice isometric camera angle
//#define LOW_Q // uncomment for low quality if your GPU is a potato

#ifdef LOW_Q
    #define MARCHSTEPS 25
#else
    #define MARCHSTEPS 50
    #define AMBIENT_OCCLUSION
    #define DOUBLE_SIDED_TRANSPARENCY
#endif
#define MAX_DIST 10.0
#define SPECULAR
#define REFLECTIONS
#define TRANSPARENCY
#define SHADOWS
#define FOG
#define DIRECTIONAL_LIGHT
#define DIRECTIONAL_LIGHT_FLARE
#define PI 3.141592654
#define kNt  -1.0 //no trans
#define kTt   1.0 //yes trans
#define kIt   0.0 //inverse trans

const float MATERIAL_1 = 1.0;
const float MATERIAL_2 = 2.0;
float gMaterial  = MATERIAL_1;

// TODO: Document these structure member fields!
// rd Ray Direction
// rl Ray Length
struct sRay   { vec3 ro ; vec3  rd ; float sd; float rl; };
struct sHit   { vec3 hp ; float hd ; vec3 oid; };
struct sSurf  { vec3 nor; vec3  ref; vec3 tra; };
struct sMat   { vec3 ctc; float frs; float smt; vec2 par; float trs; float fri; };
struct sShade { vec3 dfs; vec3  spc; };
struct sLight { vec3 rd ; vec3  col; };

// __ Matrix functions __ _____________________________________

    // Return 2x2 rotation matrix
    // With vector swizzle/mask can use as a 3x3 xform
    // For y, you need to invert
    // angle in radians
    // ========================================
    mat2 Rot2(float a ) {
        float c = cos( a );
        float s = sin( a );
        return mat2( c, -s, s, c );
    }

    // http://www.songho.ca/opengl/gl_anglestoaxes.html
    // Return 4x4 rotation X matrix
    // angle in radians
    // ========================================
    mat4 Rot4X(float a ) {
        float c = cos( a );
        float s = sin( a );
        return mat4( 1, 0, 0, 0,
                     0, c,-s, 0,
                     0, s, c, 0,
                     0, 0, 0, 1 );
    }

    // Return 4x4 rotation Y matrix
    // angle in radians
    // ========================================
    mat4 Rot4Y(float a ) {
        float c = cos( a );
        float s = sin( a );
        return mat4( c, 0, s, 0,
                     0, 1, 0, 0,
                    -s, 0, c, 0,
                     0, 0, 0, 1 );
    }

    // Return 4x4 rotation Z matrix
    // angle in radians
    // ========================================
    mat4 Rot4Z(float a ) {
        float c = cos( a );
        float s = sin( a );
        return mat4(
            c,-s, 0, 0,
            s, c, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
         );
    }

    // Translate is simply: p - d
    // opTx will do transpose(m)
    // p' = m*p
    //    = [m0 m1 m2 m3 ][ p.x ]
    //      [m4 m5 m6 m7 ][ p.y ]
    //      [m8 m9 mA mB ][ p.z ]
    //      [mC mD mE mF ][ 1.0 ]
    // ========================================
    mat4 Loc4( vec3 p ) {
        p *= -1.;
        return mat4(
            1,  0,  0,  p.x,
            0,  1,  0,  p.y,
            0,  0,  1,  p.z,
            0,  0,  0,  1
        );
    }


    // if no support for GLSL 1.2+
    //     #version 120
    // ========================================
    mat4 transposeM4(in mat4 m ) {
        vec4 r0 = m[0];
        vec4 r1 = m[1];
        vec4 r2 = m[2];
        vec4 r3 = m[3];

        mat4 t = mat4(
             vec4( r0.x, r1.x, r2.x, r3.x ),
             vec4( r0.y, r1.y, r2.y, r3.y ),
             vec4( r0.z, r1.z, r2.z, r3.z ),
             vec4( r0.w, r1.w, r2.w, r3.w )
        );
        return t;
    }


// __ Smoothing functions _____________________________________

    // Smooth Min
    // http://www.iquilezles.org/www/articles/smin/smin.htm

    // Min Polynomial
    // ========================================
    float sMinP( float a, float b, float k ) {
        float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
        return mix( b, a, h ) - k*h*(1.0-h);
    }

    // Min Exponential
    // ========================================
    float sMinE( float a, float b, float k) {
        float res = exp( -k*a ) + exp( -k*b );
        return -log( res )/k;
    }

    // Min Power
    // ========================================
    float sMin( float a, float b, float k ) {
        a = pow( a, k );
        b = pow( b, k );
        return pow( (a*b) / (a+b), 1.0/k );
    }

// __ Surface Primitives ____________________________

    // Return max component x, y, or z
    // ========================================
    float maxcomp(in vec3 p ) {
        return max(p.x,max(p.y,p.z));
    }

// Signed

    // b.x = Width
    // b.y = Height
    // b.z = Depth
    // Leave r=0 if radius not needed
    // ========================================
    float sdBox(vec3 p, vec3 b, float r) {
        vec3 d = abs(p) - b;
        return min(maxcomp(d),0.0) - r + length(max(d,0.0));
        // Inlined maxcomp
        //return min(max(d.x,max(d.y,d.z)),0.0) - r + length(max(d,0.0));
    }

    // ========================================
    float sdCappedCylinder( vec3 p, vec2 h ) {
        vec2 d = abs(vec2(length(p.xz),p.y)) - h;
        return min(max(d.x,d.y),0.0) + length(max(d,0.0));
    }

    // ========================================
    float sdCapsule( vec3 p, vec3 a, vec3 b, float r ) {
        vec3 pa = p - a, ba = b - a;
        float h = clamp( dot(pa,ba) / dot(ba,ba), 0.0, 1.0 );
        return length( pa - ba*h ) - r;
    }

    // c.x Width
    // c.y Base Radius
    // c.z Depth
    // Note: c must be normalized
    // ========================================
    float sdCone( vec3 p, vec3 c) // TODO: do we need to use 'in' for all primitives?
    {
        // c.x = length
        // c.y = base radius
        //float q = length( p.xy );
        //return dot( c, vec2( q, p.z ) ); // BUG in iq's docs -- laying on side

        float q = length( p.xz );
        return dot( c.xy, vec2( q, p.y ) );

        // Alt. cone formula given in: ???
        //vec2 q = vec2( length( p.xz ), p.y );
        //float d1 = -p.y - c.z;
        //float d2 = max( dot(q,c.xy), p.y );
        //return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.0);
    }

    // ========================================
    float sdCylinder( vec3 p, vec3 c ) {
        return length(p.xz - c.xy) - c.z;
    }

    // n.xyz = point on plane
    // n.w   = distance to plane
    // Note: N must be normalized!
    // ========================================
    float sdPlane( vec3 p, vec4 n ) {
        return dot( p, n.xyz ) + n.w;
    }

    // 4 sided pyramid
    // h.x = base X
    // h.y = height
    // h.z = base Z (usually same as h.x)
    // ========================================
    float sdPyramid4( vec3 p, vec3 h ) {
        p.xz = abs(p.xz);                   // Symmetrical about XY and ZY
        vec3 n = normalize(h);
        return sdPlane(p, vec4( n, 0.0 ) ); // cut off bottom
    }

    // ========================================
    float sdSphere( vec3 p, float r ) {
        return length(p) - r;
    }

    // ========================================
    float sdSphere2( vec3 p, float r ) {
        return abs(length(p) - r);
    }

    // ========================================
    float sdTorus( vec3 p, vec2 t ) {
        vec2 q = vec2(length(p.xy) - t.x, p.z);
        return length(q) - t.y;
    }

    // TODO: document/derive magic number 0.866025
    // ========================================
    float sdTriPrism( vec3 p, vec2 h ) {
        vec3 q = abs(p);
        return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
    }

// Unsigned

    // Box
    // ========================================
    float udBox( vec3 p, vec3 b ) {
        return length( max( abs(p) - b, 0.0 ) );
    }

    // Round Box
    // ========================================
    float udRoundBox(vec3 p, vec3 b, float r)
    {
        return length(max(abs(p) - b, 0.0))- r;
    }

// __ Distance Operations _____________________________________

// Basic
    // Op Union
    // ========================================
    float opU( float d1, float d2 ) {
        return min( d1, d2 );
    }

    // Op Union
    // ========================================
    vec4 opU2( vec4 d1, vec4 d2 ) {
        return min( d1, d2 );
    }

    // Op Union
    // ========================================
    vec4 opU( vec4 a, vec4 b ) {
        return mix(a, b, step(b.x, a.x));
    }

    // Op Subtraction
    // ========================================
    float opS( float a, float b ) {
        return max( -b, a ); // BUG in iq's docs: -a, b
    }
    // Op Subtraction
    // ========================================
    vec4 opS( vec4 a, vec4 b ) {
        return max( -b, a );
    }

    // Op Intersection
    // ========================================
    float opI( float a, float b ) {
        return max( a, b );
    }

    // Op Intersection
    // ========================================
    vec4 opI( vec4 a, vec4 b ) {
        return max( a, b );
    }

// Advanced
    // ========================================
    float opBlend( float a, float b, float k ) {
        return sMin( a, b, k );
    }

    // a angle
    // ========================================
    float displacement( vec3 p, float a ) {
        return sin(a*p.x)*sin(a*p.y)*sin(a*p.z); // NOTE: Replace with your own!
    }

    // ========================================
    float opDisplace( vec3 p, float d1, float d2 ) {
        return d1 + d2;
    }

    // Op Union Translated
    // ========================================
    vec4 opUt( vec4 a, vec4 b, float fts ){
        vec4 vScaled = vec4(b.x * (fts * 2.0 - 1.0), b.yzw);
        return mix(a, vScaled, step(vScaled.x, a.x) * step(0.0, fts));
    }


// __ Domain Operations _______________________________________
// Basic

    // Op Repetition
    // ========================================
    vec3 opRep( vec3 p, vec3 spacing ) {
        return mod(p,spacing) - 0.5*spacing;
    }

// Deformations

    // Op Twist X
    // ========================================
    vec3 opTwistX( vec3 p, float angle ) {
        mat2 m = Rot2( angle * p.x );
        return   vec3( m*p.yz, p.x );
    }

    // Op Twist Y
    // ========================================
    vec3 opTwistY( vec3 p, float angle ) {
#if 0 // original
        float c = cos( angle * p.y );
        float s = sin( angle * p.y );
        mat2  m = mat2( c, -s, s, c );
        vec3  q = vec3( m*p.xz, p.y );
        // return primitive(q); // BUG in iq's docs, should be: return q
        return q;
#else // cleaned up
        mat2 m = Rot2( angle * p.y );
        return   vec3( m*p.xz, p.y );
#endif
    }

    // Op Twist Z
    // ========================================
    vec3 opTwistZ( vec3 p, float angle ) {
        mat2 m = Rot2( angle * p.z );
        return   vec3( m*p.xy, p.z );
    }

    // iq's bend X
    // ========================================
    vec3 opCheapBend( vec3 p, float angle ) {
#if 0 // original // broken :-(
        float c = cos( angle * p.y );
        float s = sin( angle * p.y );
        mat2  m = mat2( c, -s, s, c );
        vec3  q = vec3( m*p.xy, p.z ); // BUG in iq's docs, should be: p.yx
#else
        mat2  m = Rot2( angle * p.y );
        vec3  q = vec3( m*p.yx, p.z );
#endif
        return q;
    }

    // Op Cheap Bend X
    // ========================================
    vec3 opBendX( vec3 p, float angle ) {
        mat2 m = Rot2( angle * p.y );
        return   vec3( m*p.yx, p.z );
    }

    // Op Cheap Bend Y
    // ========================================
    vec3 opBendY( vec3 p, float angle ) {
        mat2 m = Rot2( angle * p.z );
        return   vec3( m*p.zy, p.x );
    }

    // Op Cheap Bend Z
    // ========================================
    vec3 opBendZ( vec3 p, float angle ) {
        mat2 m = Rot2( angle * p.x );
        return   vec3( m*p.xz, p.y );
    }

    // d = distance to move
    // ========================================
    vec3 opTrans( vec3 p, vec3 d ) {
        return p - d;
    }

    // Note: m must already be inverted!
    // TODO: invert(m) transpose(m)
    // Op Rotation / Translation
    // ========================================
    vec3 opTx( vec3 p, mat4 m ) {   // BUG in iq's docs, should be q
        return (transposeM4(m)*vec4(p,1.0)).xyz;
    }

    // Op Scale
    // ========================================
    float opScale( vec3 p, float s ) {
        return sdBox( p/s, vec3(1.2,0.2,1.0), 0.01 ) * s; // TODO: FIXME: NOTE: replace with primative sd*()
    }

// The fun starts here!
#define BOX_W  1.50
#define BOX_H  0.25
#define BOX_D  0.75
#define RADIUS 0.666

float draw( vec3 p )
{
    vec3  q = p  ; // save original point
    float d = 0.0; // distance function; default to no intersection
    p = q;
    d = sdSphere( p, 0.5 ); // position, radius
    float b1 = sdBox(p+vec3(.0,.0,.0), vec3(.1,1.1,.7), 0.);
    float b2 = sdBox(opTx(p,Rot4X(radians(45.)))+vec3(.0,.0,.0), vec3(.7,.01,.9), 0.);
    return opU(opU(d, b1), b2);
}

// ========================================
vec4 DE( vec3 hp, float fts ) {
    vec4 vResult = vec4(MAX_DIST, -1.0, 0.0, 0.0);
    vec4 vDist = vec4( draw(hp), MATERIAL_1, hp.xz);
    vDist.y = gMaterial; // v0.42 draw may over-ride material
    return opUt(vResult, vDist, fts);
}


// ========================================
sMat getMaterial( sHit hitInfo ) {
    sMat mat;
    if(hitInfo.oid.x == MATERIAL_1) {
        mat.frs = 0.31;
        mat.smt = 1.0;
        mat.trs = 1.0;
        mat.fri = 0.75;
        const float fExtinctionScale = 2.0;
        vec3 tc = vec3(0.93,0.96,1.0);        //tex/col
        mat.ctc = (vec3(1.0) - tc) * fExtinctionScale;
    } else
    if(hitInfo.oid.x == MATERIAL_2) {
        mat.frs = 0.0;
        mat.smt = 1.0;
        mat.trs = 0.0;
        mat.fri = 0.0;
        mat.ctc = vec3(0.25,0.5,0.75); // Beautiful Baby Blue
    }
    return mat;
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


// Circle
vec3 circ(in float xp, in float yp, in float r, in vec3 col, in vec2 p)
{
    float asp = 1.77;
	// correct aspect ratio transfer
	xp = xp * asp;
	float dist = sqrt( (p.x*asp - xp) * (p.x*asp - xp) + (p.y - yp) * (p.y - yp) );
	vec3 cl = dist < r ? col : vec3(0.);
	return cl;
}

// ========================================
vec3 getBackground( vec3 rd, vec2 p ) {
    vec3  tc = vec3(noise2f(tan(rd.xy*sin(time/10.)*70.)),0.,0.)*2.;
    tc += circ(mod(time/10.,1.),sin(time),0.4,vec3(0.,1.,0.),rd.xy)*10. * mod(time*10.,1.);
    return tc;
}

// ========================================
sLight getDirLight() {
    sLight result;
    result.rd  = normalize(vec3(-0.2, -0.3, 0.5));
    result.col = vec3(8.0, 7.5, 7.0);
    return result;
}

// ========================================
vec3 getAmbient( vec3 nor, vec2 p) {
    return getBackground(nor, p);
}

// ========================================
vec3 normal( vec3 p, float fts ) {
    vec3 e = vec3(0.01,-0.01,0.0);
    return normalize( vec3(
        e.xyy*DE(p+e.xyy,fts).x +
        e.yyx*DE(p+e.yyx,fts).x +
        e.yxy*DE(p+e.yxy,fts).x +
        e.xxx*DE(p+e.xxx,fts).x)
    );
}

// ========================================
void march( sRay ray, out sHit res, int maxIter, float fts ) {
    res.hd = ray.sd;
    res.oid.x = 0.0;

    for( int i=0;i<=MARCHSTEPS;i++ ) {
        res.hp = ray.ro + ray.rd * res.hd;
        vec4 r = DE( res.hp, fts );
        res.oid = r.yzw;
        if((abs(r.x) <= 0.01) || (res.hd >= ray.rl) || (i > maxIter))
            break;
        res.hd = res.hd + r.x;
    }
    if(res.hd >= ray.rl) {
        res.hd = MAX_DIST;
        res.hp = ray.ro + ray.rd * res.hd;
        res.oid.x = 0.0;
    }
}

// ========================================
float getShadow( vec3 hp, vec3 nor, vec3 lrd, float d ) {
#ifdef SHADOWS
    sRay ray;
    ray.rd = lrd;
    ray.ro = hp;
    ray.sd = 0.05 / abs(dot(lrd, nor));
    ray.rl = d - ray.sd;
    sHit si;
    march(ray, si, 32, kNt);
    float s = step(0.0, si.hd) * step(d, si.hd );
    return s;
#else
    return 1.0;
#endif
}

// ========================================
float getAmbientOcclusion( sHit hi, sSurf s ) {
#ifdef AMBIENT_OCCLUSION
    vec3 hp = hi.hp;
    vec3 nor = s.nor;
    float ao = 1.0;

    float d = 0.0;
    for( int i=0; i<=5; i++ ) {
        d += 0.1;
        vec4 r = DE(hp + nor * d, kNt);
        ao *= 1.0 - max(0.0, (d - r.x) * 0.2 / d );
    }
    return ao;
#else
    return 1.0;
#endif
}

// ========================================
vec3 getFog( vec3 color, sRay ray, sHit hi, vec2 p ) {
#ifdef FOG
    float a = exp(hi.hd * - 0.05);
    vec3 fog = getBackground(ray.rd, p);

    #ifdef DIRECTIONAL_LIGHT_FLARE
        sLight lig = getDirLight();
        float f = clamp(dot(-lig.rd, ray.rd), 0.0, 1.0);
        fog += lig.col * pow(f, 10.0);
    #endif

    color = mix(fog, color, a);
#endif

    return color;
}

// More complex empirically motivated phase functions are efficiently approximated by the Schluck function [BLS93].
// ========================================
float getSchlick(vec3 nor, vec3 v, float frs, float sf) {
    float f = dot(nor, -v);
    f = clamp((1.0 - f), 0.0, 1.0);
    float fDotPow = pow(f, 5.0);
    return frs + (1.0 - frs) * fDotPow * sf;
}

// http://en.wikipedia.org/wiki/Fresnel_equations
// ========================================
vec3 getFresnel( vec3 dif, vec3 spe, vec3 nor, vec3 v, sMat m ) {
    float f = getSchlick(nor, v, m.frs, m.smt * 0.9 + 0.1);
    return mix(dif, spe, f);
}

// ========================================
float getPhong( vec3 ird, vec3 lrd, vec3 nor, float smt ) {
    vec3  v  = normalize(lrd - ird);
    float f  = max(0.0, dot(v, nor));
    float sp = exp2(4.0 + 6.0 * smt);
    float si = (sp + 2.0) * 0.125;
    return pow(f, sp) * si;
}

// ========================================
sShade setDirLight( sLight l, vec3 p, vec3 d, vec3 nor, sMat m ) {
    sShade s;
    vec3 lrd = -l.rd;
    float sf = getShadow( p, nor, lrd, 8.0 );
    vec3 il = l.col * sf * max(0.0, dot(lrd, nor));
    s.dfs = il;
    s.spc = getPhong( d, lrd, nor, m.smt ) * il;
    return s;
}

// ========================================
vec3 setColor( sRay ray, sHit hi, sSurf sc, sMat m, vec2 p) {
    vec3 color;
    sShade s;
    s.dfs = vec3(0.0);
    s.spc = vec3(0.0);
    float ao = getAmbientOcclusion(hi, sc);
    vec3 al = getAmbient(sc.nor, p) * ao;
    s.dfs += al;
    s.spc += sc.ref;

#ifdef DIRECTIONAL_LIGHT
    sLight dl = getDirLight();
    sShade sh = setDirLight(dl, hi.hp, ray.rd, sc.nor, m);
    s.dfs += sh.dfs;
    s.spc += sh.spc;
#endif
    vec3 dr = s.dfs * m.ctc;
    dr = mix(dr, sc.tra, m.trs);
#ifdef SPECULAR
    color = getFresnel(dr , s.spc, sc.nor, ray.rd, m);
#else
    color = dr;
#endif
    return color;
}

// ========================================
vec3 getColor( sRay ray, vec2 p ) {
    sHit hi;
    march(ray, hi, 32, kNt);
    vec3 color;

    if(hi.oid.x < 0.5) {
        color = getBackground(ray.rd, p);
    } else {
        sSurf s;
        s.nor  = normal(hi.hp, kNt);
        sMat m = getMaterial( hi );
        s.ref  = getBackground(reflect(ray.rd, s.nor), p);
        m.trs  = 0.0;
        color  = setColor(ray, hi, s, m, p);
    }

    color = getFog(color, ray, hi, p);
    return color;
}

// ========================================
vec3 getReflection( sRay ray, sHit hitInfo, sSurf s, vec2 p ) {
#ifdef REFLECTIONS
    sRay rRay;
    rRay.rd = reflect(ray.rd, s.nor);
    rRay.ro = hitInfo.hp;
    rRay.rl = 16.0;
    rRay.sd = 0.1 / abs(dot(rRay.rd, s.nor));
    return getColor(rRay, p);
#else
    return getBackground(reflect(ray.rd, s.nor), p);
#endif
}

// ========================================
vec3 getTransparency( sRay ray, sHit hit, sSurf s, sMat m, vec2 p ) {
#ifdef TRANSPARENCY
    sRay rRay;
    rRay.rd = refract(ray.rd, s.nor, m.fri);
    rRay.ro = hit.hp;
    rRay.rl = 16.0;
    rRay.sd = 0.05 / abs(dot(rRay.rd, s.nor));

    sHit hit2;
    march(rRay, hit2, 32, kIt);
    vec3 nor = normal(hit2.hp, kIt);
    sRay rRay2;
    rRay2.rd = refract(rRay.rd, nor, 1.0 / m.fri);
    rRay2.ro = hit2.hp;
    rRay2.rl = 16.0;
    rRay2.sd = 0.0;
    float ed = hit2.hd;
    vec3 color = getColor( rRay2, p );

    return color * clamp(exp(-(m.ctc * ed)),0.0,1.0);
#else
    return getBackground(reflect(ray.rd, s.nor), p);
#endif
}

// ========================================
vec3 getRayColor( sRay ray, vec2 p ) {
    sHit i;
    march(ray, i, MARCHSTEPS, kTt); //256

    vec3 color;
    if(i.oid.x < 0.5) {
        color = getBackground(ray.rd, p);
    } else  {
        sSurf s;
        s.nor  = normal(i.hp, kTt);
        sMat m = getMaterial( i );
        s.ref  = getReflection(ray, i, s, p);
        if(m.trs > 0.0) s.tra = getTransparency(ray, i, s, m, p);
        color  = setColor(ray, i, s, m, p);
    }

    getFog(color, ray, i, p); // BUG? Is this intentional that color is not updated??
    return color;
}

// ========================================
sRay setCameraRay( vec3 hp, vec3 i , vec2 fragCoord) {
    float fRatio = resolution.x / resolution.y; // Aspect Ratio

    vec3 f   = normalize(i - hp);
    vec3 vUp = vec3(0.0, 1.0, 0.0);
    vec2 vvc = 2.*fragCoord.xy/resolution.xy-1.;
    vvc.y /= fRatio;

    sRay ray;
    ray.ro = hp;
    vec3 r = normalize(cross(f, vUp));
    vUp    = cross(r, f);
    ray.rd = normalize( r * vvc.x + vUp * vvc.y + f);
    ray.sd = 0.0;
    ray.rl = MAX_DIST;
    return ray;
}

void main(void ) {
    vec2 m = vec2(0.0); // Default OpenGL camera: Look down -z axis
    vec2 p = vec2( gl_FragCoord.x / resolution.x, gl_FragCoord.y / resolution.y);

#ifdef VIEW_ISOMETRIC
    m = vec2( 3.5, 1.0 ) / PI; // fake isoemetric
#else
    m += 2. * 1. / resolution.xy;
    m.x += 1.;
#endif
    float nRotate = time *0.05; // slow rotation
    float h  = PI * (m.x - nRotate);
    float e  = mix(0.0, 2.5, m.y); // eye
    float d  = mix(2.5, 2.5 + (1.0 > 0.0 ? 4.0 : 2.0), m.y); // eye distance

    // ro RayOrigin
    vec3 ro  = vec3(sin(h) * cos(e), sin(e), cos(h) * cos(e)) * d;
    vec3 ta  = vec3(0.0, 0.0, 0.0);

    sRay ray = setCameraRay( ta + ro, ta, gl_FragCoord.xy);
    vec3 col = getRayColor( ray, p);
    gl_FragColor = vec4( col, 1.0 );
}
