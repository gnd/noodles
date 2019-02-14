//precision mediump float;

//////////////////// GLOBALS /////////////////
uniform float time;
uniform vec2 resolution;
uniform sampler2D backbuffer;
uniform vec2 rand;
vec2 p; //relative position
vec3 c; //color
float mw = 0.001; //minimal width
float asp; //aspect ratio

/// RAYMARCH
#define PI 3.1415926535898
const float eps = 0.005;
const int maxIterations = 256;
const float stepScale = 0.5;
const float stopThreshold = 0.005;
const float clipNear = 0.0;
const float clipFar = 4.0;

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
vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cg = vec3(0.,1.,0.);
vec3 cb = vec3(0.,0.,1.);
vec3 ck = vec3(0.);
vec3 ca, cc, cd, ce, cf;
//////////////////// END OF GLOBALS //////////////

/////////////////// MIDIMIX - DUAL//////////////////////
uniform float m8,m9,m10,m11;
uniform float m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31;
uniform float m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51;
uniform float m52,m53,m54,m55,m56,m57,m58,m59,m60,m61,m62,m63;
float mx11 = m12/127.;
float mx12 = m13/127.;
float mx13 = m14/127.;
float mx14 = m15/127.;
float mx15 = m16/127.;
float mx16 = m17/127.;
float mx21 = m18/127.;
float mx22 = m19/127.;
float mx23 = m20/127.;
float mx24 = m21/127.;
float mx25 = m22/127.;
float mx26 = m23/127.;
float mx31 = m24/127.;
float mx32 = m25/127.;
float mx33 = m26/127.;
float mx34 = m27/127.;
float mx35 = m28/127.;
float mx36 = m29/127.;
float mx41 = m30/127.;
float mx42 = m31/127.;
float mx43 = m32/127.;
float mx44 = m33/127.;
float mx45 = m34/127.;
float mx46 = m35/127.;
float mx51 = m36/127.;
float mx52 = m37/127.;
float mx53 = m38/127.;
float mx54 = m39/127.;
float mx55 = m40/127.;
float mx56 = m41/127.;
float mx61 = m42/127.;
float mx62 = m43/127.;
float mx63 = m44/127.;
float mx64 = m45/127.;
float mx65 = m46/127.;
float mx66 = m47/127.;
float mx71 = m48/127.;
float mx72 = m49/127.;
float mx73 = m50/127.;
float mx74 = m51/127.;
float mx75 = m52/127.;
float mx76 = m53/127.;
float mx81 = m54/127.;
float mx82 = m55/127.;
float mx83 = m56/127.;
float mx84 = m57/127.;
float mx85 = m58/127.;
float mx86 = m59/127.;
float mx91 = m60/127.;
float mx92 = m61/127.;
float mx93 = m62/127.;
float mx94 = m63/127.;
///////////////////////// END OF MIDIMIX - DUAL /////////////////

//#include primitives_sandbox.glsl

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

vec3 rd_gen(vec3 camPos, float aspect_in) {
        vec2 aspect = vec2(asp, aspect_in*70.);
        vec2 screenCoords = (2.0*gl_FragCoord.xy/resolution.xy - 1.0)*aspect*.5;
        vec3 lookAt = vec3(0.,0.,0.); 
        vec3 forward = normalize(lookAt-camPos);
        vec3 right = normalize(vec3(forward.z, 0., -forward.x ));
        vec3 up = normalize(cross(forward,right));
        float FOV = .9;
        return normalize(forward + FOV*screenCoords.x*right + FOV*screenCoords.y*up);
}

vec3 lights(vec3 camPos, vec3 sp, vec3 surfNormal) {
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
        return sceneColor*1.;
}


void main(void) {
	/// INIT 
	p = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
	vec2 pp = p;	//save some clean p 
	asp = resolution.x / resolution.y;

	//c+=circ(0.5, 0.5, 0.2, cr);	

	float type = mx12; 		// type denotes the scene that is rendered
	float rotation = mx81*360.; 	// rotation in degrees
	float distort = mx16; 		// height of sinx bumps
	float colormod = mx15; 		// some kind of color setting
	float aspect = mx14; 		// aspect ratio (maybe useless)

	mx11*=10.; 			// camera position
	vec3 camPos = vec3(sin(mx11), sin(mx11), cos(mx11)); //static vec3(0., 0., -1.); //rotate vec3(sin(mx11), sin(mx11), cos(mx11)); 
	vec3 rd = rd_gen(camPos, aspect);
	float dist = rayMarching(camPos, rd, clipNear, clipFar, type, rotation, distort );
	float outside = 0.;
        if ( dist >= clipFar ) {
	    outside = 1.;
        }
	vec3 sp = camPos + rd*dist*3.5;
	vec3 surfNormal = getNormal(sp*sp, type, rotation, distort);
	if (outside == 0.) {
	//	c += lights(camPos, sp, surfNormal);
	}

//	c+=circ_del(0.8, 0.5, mx24, c, cb);	

	
	float _xp = 0.5;
	float _w = 0.1;
	float _freq = mx24;
	float _amp = mx25;
	float cl = .5*_w / (4. * (p.x-_xp) + sin(p.y*_freq*100.)*_amp);
        //c+=cl;

	//rnd lines, keep then small
	//c += 0.001 / abs(p.x-0.5 + rnd(vec2(p.y+p.x,0.1))*mx24);
	//c += 0.001 / abs(p.x-0.5 + rnd(vec2(p.y-p.x,0.1))*mx24);

	//horizon
	//c += mod(0.001 / abs(ppp.y - mx26),0.001)*1000.;
	//c += mod(abs(ppp.y - mx26),.1)*10.;	// bigger, not 1/x

	// maybe nice lines
	//float amod = mx24;
        //c += 0.001 / abs( (p.x-0.5) + mod(p.y,amod));
        //c += 0.001 / abs(mod(p.y,amod)) * step(0.9,0.1 / abs(p.x -0.5 ));

	// Z !!
	//float amod = mx24;
        //c += 0.001 / abs( (p.x-0.5) + mod(p.y,amod));
        //c += 0.001 / abs(mod(p.y,amod)) * step(0.99,amod / (1.0-p.x -0.5 ));

	// zrpty - first
	//float vert_rep = mx24;
	//_w = mx21;
	//_amp = mx26;
	//_freq = mx25*100.;
	//float mod_count = mx34*10.;
	//c += mod(_w/10. / abs( (p.x-0.5) + mod(p.y,amod) + sin(p.y*_freq)*_amp) + _w/10. / abs(mod(p.y,amod)) * step(0.99,amod / (1.0-p.x -0.5 +sin(p.y*_freq)*-_amp)),1./mod_count)*mod_count; 

	// retty cool
	p = dx_p1(-m0+time,int(mx35*10.),mx36);
	c += zrpty(0.5, mx21, mx24, mx25, mx26*100., mx31*100.);
	//c -= zrpty(0.9, mx21, mx24, mx25, mx26*100., mx31*100.); // this is a perfect second layer

	//c -= zrptx(0.5, mx21, mx24, mx25, mx26*100., mx31*100.);

	// standard perlin noise
	//c+=cnoise(vec2(p.x,rand.x)*time)*0.5;


	// perlin line
	//c+= 0.001 / abs((p.x-0.5) + cnoise(vec2(p.y*60.,rand.y))*0.1 +cnoise(vec2(p.y*6.,rand.x))*0.2 );

	// more perlin lines
	//float sep = mx36/10.;
	//for (int i = 0; i<10; i++) {
	//	c+= 0.001 / abs((p.x-float(i)*sep - 0.3) + cnoise(vec2(p.y*60.,rand.y+float(i)))*0.1*m0 +cnoise(vec2(p.y*6.,rand.x-float(i)))*0.2*m1 );
	//}

	// not bad - as a overlay, psych filter
	//p-=0.5;
	//c+=cnoise(vec2(p*cnt(mx36*1000.)*mx35*100.));

	// perlin colour
	//p-=mx35;		//nice to move this around, maybe in pd wit rand + line
	//float pamp = mx36*1000.;
	//c+= vec3 ( cnoise(p*pamp* (rand.x+rand.y)), cnoise(p*pamp*rand.x), cnoise(p*pamp*rand.y));

	// domain distortion
	//c+=fbm(p, 10, 2., 1.);

	int it = int(mx35*10.);
	vec2 q = vec2( fbm( p + vec2(0.0,time), it, 2., mx36 ),
                     fbm( p + vec2(5.2,1.3), it, 2., mx36 ) );
	vec2 r = vec2( fbm( p + 4.0*q + vec2(1.7,9.2), it, 2., mx36 ),
                     fbm( p + 4.0*q + vec2(8.3,2.8), it, 2., mx36 ) );
      	//c += fbm( p + 4.0*r , it, 2., mx36);

	p = dx_p1(time,int(mx35*10.),mx36);
	c += circ(0.5, 0.5, mx11, cb);
	p = pp;
	
	// DRAW
	gl_FragColor = vec4(c, 1.0);
}
