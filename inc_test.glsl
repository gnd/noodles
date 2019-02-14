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
vec3 ca,cc,cd,ce,cf,cx;
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

//#include primitives.glsl

void main(void) {
	/// INIT
	p = vec2( gl_FragCoord.x / resolution.x, gl_FragCoord.y / resolution.y);
	asp = resolution.x / resolution.y;
	c= ck;

	float BLOB = mx13;
	if (BLOB == 1.0) {
		c += lights();
	}

	float RAY = 1.0;
	if (RAY == 1.) {
	 	vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-0.4;
		p_ray = vec4(p, vec2(mx55,0.))-0.5;
	 	vec4 d=p_ray;
	 	vec4 t;
	    p_ray.z += time;
	    for(float i=mx56*5.; i>0.; i-=.1) {
	        t = abs(mod(p_ray, 0.5)-.25*m0);
	        float x = max(t.y*2., length(t.xz)*2.*m1);
			if (mx52 == 0.) {
				c = vec3(i)*mx51;
			}
			if (mx52 == 1.) {
	        	c = vec3(sin(mod(i*p.y*i*10.,2.)), tan(i*i*p.y), tan(i/p.x))*mx51;
			}
			if (mx52 == 2.) {
				c = vec3(sin(100.*abs(noise2f(vec2(i*p.x*10./p.y)))))*mx51;
	        }
	        if (x<.2) break;
	        p_ray -= d*x;
	     }
	}

	float OSC=1.;
	if (OSC == 1.) {
	    float cent = 0.5;
	    float w1 = mx46;
		float w2 = mx45;

		if (mx42 >= 1.) {
	        	c+=oscx(cent, w1, m1, m2, vec3(2.,1.,.5))*mx41;
		}
		if (mx42 >= 2.) {
	        	c+=oscy(cent, w1, m1, m2, vec3(2.,1.,.5))*mx41;
		}
		if (mx42 >= 3.) {
	        	c+=oscx(m3, w2, m3, m4, vec3(.5,1.,2.))*mx41;
		}
		if (mx42 >= 4.) {
	        	c+=oscy(cent, w2, m3, m4, vec3(.5,1.,2.))*mx41;
		}
	}

	// FEEDBACK
	if (mx63 == 1.) {
		float xp = .5;
		float yp = .5;
	    float xs = mx66*m2;
	    float ys = mx65*m1;
		float bs = mx64*10.*m0;
		for (int i=0; i<10; i++) {
	    	cx= feedb(float(i)/10.+0.05, yp, xs, ys, bs, 0., 1., c);
			c = mix(c,cx,mx61);
		}
	//	c = mix(c,cx,mx61);
	}


	// BACKBUFFER
	mx75 = mx75*2.-1.;
	if (mx73 == 1.) {
		// normal / inverse y / inverse x
		if (mx72 == 1.) {
			p = vec2(p.x,1.-p.y);
		}
		if (mx72 == 2.) {
	                p = vec2(1.-p.x,p.y);
	        }
		if (mx72 == 3.) {
			p = vec2(p.y,p.x);
		}
		cx = texture2D(backbuffer, (p-mx76)*(1.-mx75)+0.5).xyz;
		//threshold
		c=mix(c,cx,mx71);
	}



	// COLOR STUFF
	cx = c;
	if (mx83 == 1.) {
		cx = vec3(dot(c.rgb, vec3(0.299, 0.587, 0.114)));
		cx = mix(cx,vec3(mx86,mx85,mx84),mx81);
	} else {
		cx = vec3(cx.r*mx86,cx.g*mx85,cx.b*mx84);
	}

	cx+=vec3(cnt(50))*mx32;


	// DRAW
	gl_FragColor = vec4(cx, 1.0);

}
