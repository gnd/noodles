#ifdef GL_ES
precision mediump float;
#endif
#define GX 6.0
#define GY 5.0

uniform float time;
uniform vec2 resolution;
vec2 position;
vec3 c;
uniform vec2 uv;
const int IMAX = int(GX * GY);
float a[IMAX];
vec4 rules;

float zsin(in float a) {
	return (sin(a) + 1.0) / 2.0;
}

float zcos(in float a) {
        return (cos(a) + 1.0) / 2.0;
}

float mdl(in float a, in float b)
{
        return a - b * floor(a/b + 0.001);
}

float cnt(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

vec3 linex(in float xpos, in vec3 col)
{
	return int(mod(position.x*100.0,100.0)) == int(xpos*100.0) ? col : vec3(0.0);
}

vec3 liney(in float ypos, in vec3 col)
{
        return int(mod(position.y*100.0,100.0)) == int(ypos*100.0) ? col : vec3(0.0);
}

vec3 linexx(in float xpos, in vec3 col)
{
        return int(mod(position.x*5.0,5.0)) == int(xpos*5.0) ? col : vec3(0.0);
}

vec3 linexxx(in float xpos, in vec3 col, float siz)
{
        return int(mod(position.x*siz,siz)) == int(xpos*siz) ? col : vec3(0.0);
}

vec3 sqr(in float xpos, in float ypos, in float size, in vec3 col)
{
	float ysize = size * 1.0;
//	vec3 cl = ((position.x > xpos) && (position.x < xpos + size) && (position.y > ypos) && (position.y < ypos + ysize)) ? col : vec3(0.0);
	vec3 cl = ((position.x*zcos(time) > xpos*zsin(time)) && (position.x < xpos + size) && (position.y > ypos*zsin(time*10.0)) && (position.y*zcos(time*10.0) < ypos + ysize)) ? col : vec3(0.0);
	return cl;
}

vec3 crc(in float xpos, in float ypos, in float rad, in vec3 col)
{
	float dist = sqrt( (position.x - xpos) * (position.x - xpos) + (position.y - ypos) * (position.y - ypos) );
	//float dist = sqrt (8.0);
	vec3 cl = dist < rad ? col : vec3(0.0);
	return cl;
} 

vec3 crcx(in float xpos, in float ypos, in float rad, in vec3 col)
{
	ypos = ypos * 1.33;
        float dist = sqrt( (position.x - xpos) * (position.x - xpos) + (position.y - ypos) * (position.y - ypos) );
        //float dist = sqrt (8.0);
        vec3 cl = dist < rad ? col : vec3(0.0);
        return cl;
}


void main( void ) {
vec3 col;
position = ( gl_FragCoord.xy / resolution.xy );

float kak = 0.1;
float kok = 0.91;
float kek = 0.5;
c += vec3(sin(time*(position.y-500.0+(500.0*cnt(100)/2.0))/sqrt(position.x/10.0*cnt(1000)))) * cnt(200);
c+= vec3(cos(sqrt(position.x*40.0)/sqrt(position.y)*position.x*1000.0*cnt(100)*cnt(100)));
c+=vec3(sin(position.y*10.0*cnt(1000)*cnt(100)));
c-=vec3(cnt(500),0.0,0.0);
vec3 kik = c;
c=linexx(cnt(100),c);
c+=linexx(cnt(200),kik);
c+=linexx(cnt(500),kik);
c-=linexxx(cnt(1000),kik,cnt(1000));
c*=linexxx(cnt(44),kik,sin(time*1000.0));
c-=vec3(sin(time*(position.y-100.0+(500.0*cnt(200)/2.0))/sqrt(position.y/10.0*cnt(1000)))) * cnt(200);

gl_FragColor = vec4(c, 1.0 );

}

