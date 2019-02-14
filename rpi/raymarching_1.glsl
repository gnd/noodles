#ifdef GL_ES
precision lowp float;
#endif

uniform float time;
uniform vec2 resolution;
uniform float m0;

void main(void) {

/// RED SHITE - uses m0
vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-0.5;
vec4 d=p_ray; vec4 t; p_ray.z += time;
for(float i=1.8; i>0.; i-=.1) {
	t = abs(mod(cos(p_ray),.3)-m0);
	float x = max(t.y*2., length(t.xz)*2.);
	gl_FragColor = vec4(i*2.,cos(i),i, 2.0);
	if(x<.5) break;
	p_ray -= d*x;
}
for(float i=1.8; i>0.; i-=.1) {
        t = abs(mod(p_ray,.3)-0.34);
        float x = max(t.x*3., length(t.yz)*3.);
        gl_FragColor -= vec4(i*0.1,i,i, 2.0);
        if(x<.5) break;
        p_ray -= d*x*time;
}
/// END OF RED SHITE


	
}
