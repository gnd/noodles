//////////////////// GLOBALS /////////////////
uniform float time;
uniform vec2 resolution;
uniform vec2 rand;
uniform vec2 mouse;
uniform sampler2D backbuffer;
uniform sampler2D video;
vec2 p,pp; //relative position
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
	vec3 cl = ((p.x > xp- s/2.) && (p.x < xp + s/2.) && (p.y > yp-s/2.) && (p.y < yp + ys/2.)) ? col : vec3(0.0);
	return cl;
}


void main( void ) {
////////////////////// INIT ////////////////////////
p = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
pp = p;
asp = resolution.x / resolution.y;
vec3 cw = vec3(1.0);
vec3 cr = vec3(1.,0.,0.);
vec3 cb = vec3(0.);

////////////////////// END OF INIT /////////////////
c=vec3(0.);

// CHECK AUTONOM
c+=liney(0.1,cnt(1000)/10., cw);
c+=liney(cnt(10000),0.01, cr);

// CHECK RANDOMX/Y
float once = cnt(1000);
if (once > 0.99) {
c+=sqr(0.2,rand.x,0.1, cw);
c+=sqr(0.8,rand.y,0.1, cw);
}

// CHECK NET
c+=sqr(0.075,.5,m0, cw);
c+=sqr(0.2,.5,m1, cr);
c+=sqr(0.325,.5,m2, cw);
c+=sqr(0.45,.5,m3, cr);
c+=sqr(0.575,.5,m4, cw);
c+=sqr(0.7,.5,m5, cr);
c+=sqr(0.825,.5,m6, cw);
c+=sqr(0.95,.5,m7, cr);

// CHECK BACKBUFFER
p = vec2(p.x,1.-p.y);
c += texture2D(backbuffer, (p-0.5)*(1.+0.1)+0.5).xyz*.9;
p = pp;

// CHECK VIDEO
//float pixelx = gl_TexCoord[0].x;
//float pixely = gl_TexCoord[0].y;
vec3 vc = texture2D(video, p);
c+=sqr(.5,.5,.8,vc);

// DRAW
gl_FragColor = vec4(c, 1.0);

}
