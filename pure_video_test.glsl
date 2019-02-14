uniform vec2 resolution;
uniform sampler2D video;

void main(void) {
	vec2 p = vec2( gl_FragCoord.x / resolution.x, 1.0 - gl_FragCoord.y / resolution.y);
	gl_FragColor = texture2D(video, p);
}
