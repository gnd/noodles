//precision mediump float;

//////////////////// GLOBALS /////////////////
uniform float time;
uniform vec2 resolution;
uniform sampler2D backbuffer;
uniform vec2 rand;
//vec2 p; //relative position
vec3 c; //color
float mw = 0.001; //minimal width
//float asp; //aspect ratio

/// RAYMARCH
#define PI 3.1415926535898
const float eps = 0.005;
const int maxIterations = 256;
const float stepScale = 0.5;
const float stopThreshold = 0.005;
const float clipNear = 0.0;
const float clipFar = 8.0;

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
vec3 ca, cc, cd, ce, cf, cx;

/// INIT 
vec2 p = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
vec2 pp = p;    //save some clean p 
float asp = resolution.x / resolution.y;
/////////////////// MIDIMIX //////////////////////
uniform float m8,m9,m10;
uniform float m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30;
uniform float m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50;
uniform float m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,m61,m62;
uniform float m63,m64,m65;
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
float mx76 = m49/126.;
float mx81 = m50/127.;
float mx82 = m51/127.;
float mx83 = m52/127.;
float mx84a = m53/127.; 
float mx85a = m54/127.;
float mx86a = m55/127.;
float mx84b = m56/127.;
float mx85b = m57/127.;
float mx86b = m58/127.;
float mxs1 = m59/127.;
float mxs2 = m60/127.;
float mxs3 = m61/127.;
float mx91 = m62/127.;
float mx92 = m63/127.;
float mx93 = m64/127.;
float mx94 = m65/127.;
///////////////// END OF MIDIMIX ////////////////////

//#include primitives_sandbox.glsl

// sinbump2
float sinbumps2(in vec3 p, float distort){
        float frq = 100.;
	float force = PI * distort;
	distort *= 0.05;
        return distort * cos(p.x/frq + force)*atan(p.y*frq + force)*tan(p.z*frq - time*100.31);
}

// bumps (for col + geom)
float bumps(in vec3 p, float type, float distort) {
	if (type == 0.) {
		return sinusoidBumps(p, distort);
	}
	if (type == 1.) {
		return sinbumps2(p, distort);
	}
	if (type == 2.) {
		return distort * sin(20.*p.x)*sin(20.*p.y)*sin(20.*p.z) * (.2*zsin(time*100.)+.8);
	}
}
	

// The whole scene
float scene(in vec3 p, float type, float distort, float bumptype, float dom_mod, float mod_dom_on) {
	p = p.yxz;
	dom_mod *= 10.;
        if (type == 0.) {
        return box(p, vec3(0.2,0.2,0.2)) + bumps(p, bumptype, distort);
        }
        if (type == 1.) {
        return sphere(p, vec3(0.,0.,0.), 0.3) + bumps(p, bumptype, distort);
        }
        if (type == 2.) {
        return hex(p, vec2(0.2,0.4)) + bumps(p, bumptype, distort);
        }
	if (type == 3.) {
		// domain repetition
		vec3 pp = vec3(mod(p.x, dom_mod), mod(p.y, dom_mod) , mod(p.z, dom_mod));
		if (mod_dom_on == 0.) {
			pp = p;
		}

		// geom params
		float sph_size = m0;
		float sph1 = sphere(pp, vec3(asin(cnt(1000)),0.,0.), sph_size);
		float sph2 = sphere(pp, vec3(-.3,sin(time),0.), sph_size);
		float sph3 = sphere(pp, vec3(.0,.3,.0), sph_size);
		float ret = smin(sph1, sph2, mx21);
		ret = smin(ret, sph3, mx21);
		float hex1 = hex(pp, vec2(0.2, 0.4));
		ret = max(-hex1, ret);

		// apply some bumps
		ret += bumps(pp, bumptype, distort);

		// second geom params
		//float hex2 = hex(pp-vec3(-cnt(500),-.5,.0), vec2(.9, cnt(1000)*0.8));
		//ret = max(hex2, ret);
		return ret;
	}
}

// Gets the surface normal for p
vec3 getNormal(in vec3 original_p, float type, float rotation, float distort, float bumptype, float dom_mod, float dom_mod_on) {
        mat4 rotmat = rmat(vec3(0.,1.,0.), radians(rotation));
        vec3 p = rot(original_p, rotmat);
        return normalize(vec3(
                scene(vec3(p.x+eps,p.y,p.z),type, distort, bumptype, dom_mod, dom_mod_on)-scene(vec3(p.x-eps,p.y,p.z), type, distort, bumptype, dom_mod, dom_mod_on),
                scene(vec3(p.x,p.y+eps,p.z),type, distort, bumptype, dom_mod, dom_mod_on)-scene(vec3(p.x,p.y-eps,p.z), type, distort, bumptype, dom_mod, dom_mod_on),
                scene(vec3(p.x,p.y,p.z+eps),type, distort, bumptype, dom_mod, dom_mod_on)-scene(vec3(p.x,p.y,p.z-eps), type, distort, bumptype, dom_mod, dom_mod_on)
        ));
}

// Raymarches
float rayMarching( vec3 origin, vec3 dir, float start, float end, float type, float rotation, float distort, float bumptype, float dom_mod, float dom_mod_on) {
        float sceneDist = 1e4;
        float rayDepth = start;
        for ( int i = 0; i < maxIterations; i++ ) {
                mat4 rotmat = rmat(vec3(0.,1.,0.), radians(rotation));
                vec3 rotated_p = rot(origin + dir * rayDepth, rotmat);
                sceneDist = scene( rotated_p, type, distort, bumptype, dom_mod, dom_mod_on );
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

vec3 lights(vec3 camPos, vec3 sp, vec3 surfNormal, float bumptype, float distort, vec3 col) {
        vec3 lp = vec3(1.5*sin(time*.5), 0.75+0.25*cos(time*0.5), -1.0);
        vec3 ld = lp-sp;
        vec3 lcolor = vec3(1.,1.,1.);
        float len = length( ld );
        ld /= len;
        float lightAtten = min( 1.0 / ( 0.25*len*len ), 1.0 ); // Keeps things between 0 and 1.
        vec3 ref = reflect(-ld, surfNormal);
        vec3 sceneColor = vec3(0.0);
        vec3 objColor = vec3(col.r, col.g, col.b);
	float c_bumps = bumps(sp, bumptype, distort);
	objColor += clamp(objColor-vec3(0.3, 0.4, 0.)*c_bumps, 0.0, 1.0);
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


void main(void) {

	// vertical
	p = vec2(1.-p.y,p.x);
	pp = p;
	c=vec3(0.);

	float type = mx13; 		// type denotes the scene that is rendered
	float rotation = mx91*360.; 	// rotation in degrees
	//rotation = zsin(time)*360;
	rotation = cnt(5000)*360.*m0;
	float distort = mx16; 		// height of sinx bumps
	float bumptype = mx23;		// type of bumps
	float colormod = mx15; 		// some kind of color setting
	float aspect = mx14; 		// aspect ratio (maybe useless)
	float dom_mod = mx26;		// size of domain modulo for (X, Z)
	float dom_mod_on = mx22;	// if domain modulo enabled
	vec3 color = vec3(mx86b, mx85b, mx84b);	// RGB color

	mx11*=10.; 			// camera position
	vec3 camPos = vec3(sin(mx11), sin(mx11), cos(mx11)); //static vec3(0., 0., -1.); //rotate vec3(sin(mx11), sin(mx11), cos(mx11)); 
	vec3 rd = rd_gen(camPos, aspect);
	float dist = rayMarching(camPos, rd, clipNear, clipFar, type, rotation, distort, bumptype, dom_mod, dom_mod_on );
	float outside = 0.;
        if ( dist >= clipFar ) {
	    outside = 1.;
        }
	vec3 sp = camPos + rd*dist*3.5*colormod;
	vec3 surfNormal = getNormal(sp*sp, type, rotation, distort, bumptype, dom_mod, dom_mod_on);
	if (outside == 0.) {
		c += lights(camPos, sp, surfNormal, bumptype, distort, color);
	}

	
	// ZRPT'S
	if (mx42 == 1.) {
	
		// pretty cool
		p = dx_p1(time,int(mx34*10.*m1),mx35*m0);

		float zw = mx41;
		float zrep = mx44;
		float zamp = mx45*m0*5.;
		float zfreq = mx46*100.;
		float zmod_amp = mx31*10.;
		if (mx43 == 0.) {
			c += zrpty(0.5, zw, zrep, zamp, zfreq, zmod_amp);
		} 
		if (mx43 == 1.) {
                        c += zrpty(0.1, zw, zrep, zamp, zfreq, zmod_amp);
			c += zrpty(0.9, zw, zrep, zamp, zfreq, zmod_amp);
                }
		if (mx43 == 2.) {
                        c += zrpty(0.5, zw, zrep, zamp, zfreq, zmod_amp);
                        c -= zrptx(0.5, zw, zrep, zamp, zfreq, zmod_amp);
                }
		if (mx43 == 3.) {
                        c += zrpty(0.1, zw, zrep, zamp, zfreq, zmod_amp);
			c -= zrpty(0.9, zw, zrep, zamp, zfreq, zmod_amp);
                        c -= zrptx(0.5, zw, zrep, zamp, zfreq, zmod_amp);
                }

	}
	p = pp;


	// FLACKS
	if (mx42 == 2.) { 
		// rotate
		float cen = .5;
		p -= vec2(cen);
		p = rot2d(mx91);
		p += vec2(cen);
		p = vec2(p.x - mx46*m0, p.y);
		p *= mx45 * 10.;

		// domain distort
		p = dx_p1(time,int(mx34*10.),mx35);

		vec3 color = vec3(mx86b, mx85b, mx84b); // RGB color
		if (mx43 == 1.) {
			c+=flacks(15, color, mx44)*mx41*5.;
		}
		if (mx43 == 2.) {
			c+=flacks2(15, color, mx44)*mx41*100.;
		}
	}
	p = pp;

	// RAYS
	float RAY = 1.0;
	if (RAY == 1.) {
		vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-0.4;
        	p_ray = vec4(p, vec2(mx55,0.))-0.5;
		vec4 d = p_ray;
		if (mx52 == 1.) {
			p = p_ray.xy;
			p = dx_p1(0.5, 5, mx54*m2);
			d=p.xyxy;
		}
		vec4 t;
		p_ray.z += time;
		for(float i=mx56*5.*m0*10.; i>0.; i-=.1)
		{
			t = abs(mod(p_ray, .5)-.25);
			float x = max(t.y*2., length(t.xz)*2.);
			if (mx53 == 0.) {
				c = vec3(i)*mx51;
			}
			if (mx53 == 1.) {
				c = vec3(sin(mod(i*p.y*i*10.,2.)), tan(i*i*p.y), tan(i/p.x))*mx51;
		}
			if (mx53 == 2.) {
				c = vec3(sin(100.*abs(noise2f(vec2(i*p.x*10./p.y)))))*mx51;
			}
			if(x<.2) break;
			p_ray -= d*x;
		}
	}
	p = pp; // never forget to put this back after domain distort

	// FEEDBACK
	if (mx62 == 1.) {
		p = p.yx;
		float xp = .5;
		float yp = .5;
        	float xs = mx66;
        	float ys = mx65;
		float bs = mx64*10.;
		if (mx63 == 0.) {
			for (int i=0; i<10; i++) {
				cx = feedb_sqr(float(i)/10.+0.05, yp, xs*rand.x, ys, bs, 0., 1., c);
				c = mix(c,cx,mx61);
			}
		}
		if (mx63 == 1.) {
			ys *= 10.;
			float spacing = (1.0/ys);
			for (int i=1; i<=int(ys); i++) {
       				cx = feedb_crc(spacing*float(i)-spacing/2., yp*rand.y, xs, bs, c);
				c = mix(c,cx,mx61);
			} 
		}
		//c = mix(c,cx,mx61);
	}

	p = pp;
	// BACKBUFFER
	if (mx72 == 1.) {
		float feedpos = mx76;			// diagonal position
		//feedpos = .5;
		float feedsize = mx75 * 2. - 1.;	// feedback size	
		// normal / inverse y / inverse x
		if (mx73 == 1.) {
			p = vec2(p.x,1.-p.y);
		}
		if (mx73 == 2.) {
                	p = vec2(1.-p.x,p.y);
        	}
		if (mx73 == 3.) {
			p = vec2(p.y,p.x);
		}
		cx = texture2D(backbuffer, (p-feedpos)*(1.-feedsize)+0.5).xyz;
		// inverse colors
		if (mx74 > .5) {
			cx = 1.-cx;
			cx *= cx*cx;
		}
		c=mix(c,cx,mx71);
	}


	// COLOR STUFF
	cx = c;
	if (mx82 == 1.) { 
		cx = vec3(dot(c.rgb, vec3(0.299, 0.587, 0.114)));
		cx = mix(cx,vec3(mx86a,mx85a,mx84a),mx81);		//this is not really doing colorize
	} else {
		cx = vec3(cx.r*mx86a,cx.g*mx85a,cx.b*mx84a);
	}

	if (mx33 == 1.) {
		cx += cnt(50);
	}

	// DRAW
	gl_FragColor = vec4(cx, 1.0);

}
