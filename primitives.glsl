/////////////////////////////////////// PRIMITIVES START
vec2 p = vec2(0.);
float asp = 1.;
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
float cnt(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

////////////////////////////////////// SIMPLE PRIMITIVES

// Vertical line
vec3 linex(in float xp, in float width, in vec3 col)
{
        return ((p.x <= xp+width/2.) && (p.x >= xp-width/2.)) ? col : vec3(0.0);
}

// Horizontal line
vec3 liney(in float yp, in float width, in vec3 col)
{
        return ((p.y <= yp+width/2.0) && (p.y >= yp-width/2.0)) ? col : vec3(0.0);

}

// Square
vec3 sqr(in float xp, in float yp, in float s, in vec3 col)
{
	float ys = s * asp;
	vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-s/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
	return cl;
}

// Rectangle
vec3 rect(in float xp, in float yp, in float sx, in float sy, in vec3 col)
{
	float syy = sy*asp;
	vec3 cl = ((p.x > xp - sx/2.) && (p.x < xp + sx/2.) && (p.y > yp-syy/2.) && (p.y < yp + syy/2.)) ? col : vec3(0.0);
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
vec3 oscx(in float xp, in float w, in float sinm, in float sinm2, in vec3 c)
{
	float cl = 0.5*w / abs(4. * (p.x - xp) + sin(p.y * sinm * 100.)*sinm2);
	return cl * c;
}

// Oscillator ??
vec3 oscxt(in float xp, in float w, in float sinm, in float sinm2, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.x - xp) + sin(p.y * 10. * sinm) + tan(p.y * sinm * 150.)*sinm2);
        return cl * c;
}

// Oscillator on the y-axis
vec3 oscy(in float yp, in float w, in float sinm, in float sinm2, in vec3 c)
{
        float cl = 0.5*w / abs(4. * (p.y - yp) + cos(p.x * sinm * 100.)*sinm2);
        return cl * c;
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


/////////////////////////////////////////// GEOMETRIC MODIFIERS
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
    return vec3(q.x,q.y,q.z);
}

float opS( float d1, float d2 )
{
    return max(-d2,d1);
}

vec2 opU( vec2 d1, vec2 d2 )
{
	return (d1.x<d2.x) ? d1 : d2;
}

///////////////////////////////////////////// RAYMARCHING PRIMITIVES

// Sphere field
float sphere(in vec3 p, in vec3 centerPos, float radius) {
        return length(p-centerPos) - radius;
}

// Box field
float box( vec3 p, vec3 b ) {
	return length(max(abs(p)-b,0.0));
}

float hex( vec3 p, vec2 h )
{
	vec3 q = abs(p);
	return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
}

// Sinusoid bumps
float sinusoidBumps(in vec3 p){
        float frq = 10.;
        return 3.*sin(p.x*frq+time*0.57)*atan(p.y*frq+time*2.17)*sin(p.z*frq-time*1.31);
}

// The whole scene
float scene(in vec3 p, float type, float distort) {
	if (type == 0.) {
	return box(p, vec3(0.2,0.2,0.2)) + distort*sinusoidBumps(p);
	}
	if (type == 1.) {
	return sphere(p, vec3(0.,0.,0.), 0.3) + distort*sinusoidBumps(p);
        }
	if (type == 2.) {
	return hex(p, vec2(0.2,0.4)) + distort*sinusoidBumps(p);
        }
}

// Gets the surface normal for p
vec3 getNormal(in vec3 original_p, float type, float rotation, float distort) {
	mat4 rotmat = rmat(vec3(0.,1.,0.), radians(rotation));
        vec3 p = rot(original_p, rotmat);
        return normalize(vec3(
                scene(vec3(p.x+eps,p.y,p.z),type, distort)-scene(vec3(p.x-eps,p.y,p.z), type, distort),
                scene(vec3(p.x,p.y+eps,p.z),type, distort)-scene(vec3(p.x,p.y-eps,p.z), type, distort),
                scene(vec3(p.x,p.y,p.z+eps),type, distort)-scene(vec3(p.x,p.y,p.z-eps), type, distort)
        ));
}

// Raymarches
float rayMarching( vec3 origin, vec3 dir, float start, float end, float type, float rotation, float distort) {
        float sceneDist = 1e4;
        float rayDepth = start;
        for ( int i = 0; i < maxIterations; i++ ) {
		mat4 rotmat = rmat(vec3(0.,1.,0.), radians(rotation));
  		vec3 rotated_p = rot(origin + dir * rayDepth, rotmat);
                sceneDist = scene( rotated_p, type, distort );
                if (( sceneDist < stopThreshold ) || (rayDepth >= end)) {
                        break;
                }
                rayDepth += sceneDist * stepScale;
        }
        if ( sceneDist >= stopThreshold ) rayDepth = end;
        else rayDepth += sceneDist;
        return rayDepth;
}

// lighting
vec3 lights(float type, float rotation, float distort, float colormod, float aspect_in) {
        vec2 aspect = vec2(asp, aspect_in*70.);
        vec2 screenCoords = (2.0*gl_FragCoord.xy/resolution.xy - 1.0)*aspect*.5;
        vec3 lookAt = vec3(0.,0.,0.);
        //rotate by hand
        mx11*=10.;
        vec3 camPos = vec3(sin(mx11), sin(mx11), cos(mx11));
        //static
        //vec3 camPos = vec3(0., 0., -1.);
        vec3 forward = normalize(lookAt-camPos);
        vec3 right = normalize(vec3(forward.z, 0., -forward.x ));
        vec3 up = normalize(cross(forward,right));
        float FOV = .9;
        vec3 ro = camPos;
        vec3 rd = normalize(forward + FOV*screenCoords.x*right + FOV*screenCoords.y*up);
        vec3 bgcolor = vec3(1.,0.97,0.92)*0.15;
        float dist = rayMarching(ro, rd, clipNear, clipFar, type, rotation, distort );
        vec3 BACK = vec3(1.1);
        if ( dist >= clipFar ) {
            //c += vec3(bgcolor);
            //gl_FragColor = vec4(bgcolor, .3);
            //BACK = vec3(bgcolor);
            c=ck;
            //return;
        }
        vec3 sp = ro + rd*dist*colormod*10.;
        vec3 surfNormal = getNormal(sp*sp, type, rotation, distort);
        //vec3 lp = vec3(1.5*sin(time*.5), 0.75+0.25*cos(time*0.5), -1.0);
        vec3 lp = vec3(noise2f(p),noise2f(p*p),noise2f(p*p*p)*m0);
        vec3 ld = lp-sp;
        vec3 lcolor = vec3(1.,1.,1.);
        float len = length( ld );
        ld /= len;
        float lightAtten = min( 1.0 / ( 0.25*len*len ), 1.0 ); // Keeps things between 0 and 1.
        vec3 ref = reflect(-ld, surfNormal);
        vec3 sceneColor = vec3(0.0);
        vec3 objColor = vec3(1., .0, 0.);
        float bumps =  sinusoidBumps(sp);
        objColor += clamp(objColor-vec3(0.3, 0.4, 0.)*bumps, 0.0, 1.0);
        float ambient = 1.3;
        float specularPower = 100.0;
        float diffuse = max( 0.0, dot(surfNormal, ld) );
        diffuse = pow(diffuse, 1000.);
        float specular = max( 0.0, dot( ref, normalize(camPos-sp)) );
        specular = pow(specular, specularPower);
        sceneColor = vec3(0.);
        sceneColor += (objColor*(diffuse*0.8+ambient)+specular*0.5)*lcolor*lightAtten;
	return sceneColor;
}

////////////////////////////////////////////// PRIMITIVES END
