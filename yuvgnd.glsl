uniform sampler2D video;
uniform sampler2D backbuffer;
uniform vec2 resolution;
uniform vec2 rand;
uniform float time;
vec2 p, pp;
float asp = resolution.x / resolution.y;

//#include primitives_sandbox.glsl

vec3 yuv2rgb(in float y, in float u, in float v) {
	float r,g,b;
        u -= 0.5;
        v -= 0.5;
        y = (y - 0.0625) * 1.1643;
        r = clamp(y + 1.0958 * v, 0., 1.);
        g = clamp(y - 0.8129 * v - 0.39173 * u, 0., 1.);
        b = clamp(y + 2.017 * u, 0., 1.);
        return vec3(r,g,b);
}

vec3 vid(in float posx, in float posy, in float sx, in float sy) {
	vec3 rgb = vec3(0.);
        float y, u, v;
	float px = gl_FragCoord.x - posx * resolution.x;
	float py = sy*resolution.y - (gl_FragCoord.y - posy * resolution.y);
	sx *= resolution.x; 
	sy *= resolution.y;
	vec2 size = vec2(sx,sy);
	if ((px > 0. && px < sx) && (py > 0. && py < sy)) {
        	if (mod(px, 2.) == 0.01) {
                	y = texture2D(video, vec2(px,py) / size).r;
                	u = texture2D(video, vec2(px,py) / size).g;
                	v = texture2D(video, vec2(px+1.,py) / size).g;
        	} else {
                	y = texture2D(video, vec2(px,py) / size).r;
                	u = texture2D(video, vec2(px-1.,py) / size).g;
                	v = texture2D(video, vec2(px,py) / size).g;
        	} 
        	rgb = yuv2rgb(y,u,v);
	}
	return rgb;
}
	

void main()
{
	p = vec2( gl_FragCoord.x / resolution.x, 1. - gl_FragCoord.y / resolution.y);
	vec3 cr = vec3(1.,0.,0.);
	// ADD position
	//y0u0, y1v0, y2u2, y3v2

	float y, u, v;
	float px = gl_FragCoord.x;
	float py = resolution.y - gl_FragCoord.y;
	if (mod(px, 2.) == 0.01) {
		y = texture2D(video, vec2(px,py) / resolution).r;
		u = texture2D(video, vec2(px,py) / resolution).g;	
		v = texture2D(video, vec2(px+1.,py) / resolution).g;
	} else {
		y = texture2D(video, vec2(px,py) / resolution).r;
                u = texture2D(video, vec2(px-1.,py) / resolution).g; 
                v = texture2D(video, vec2(px,py) / resolution).g;
	}	

	vec3 c = vec3(0.);
	c += yuv2rgb(y,u,v);

	c += vid(0.1, 0.8, 0.1, .1);
	c += circ_del(0.1, 0.9, 0.1, c, cr);

	gl_FragColor = vec4(c, 1.);

}
