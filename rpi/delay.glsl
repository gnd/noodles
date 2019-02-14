/* -*- Mode: c; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- */

#ifdef GL_ES
precision mediump float;
#endif


uniform vec2 resolution;
uniform sampler2D backbuffer;
uniform sampler2D prev_layer;
uniform float time;
uniform float sa;
uniform float sb;

// counter 0 to 1
float cnt(in int m)
{
        float divider = 1000.0 / float(m)*10.0;
        return mod(time*divider, 10.0) / 10.0;
}


//float del = cnt(10000)*1.;
//float del = cnt(500)*1.;
//base
float del = 0.5;
void main(void) {
    //vec2 uv = (gl_FragCoord.xy - vec2(sin(time*10.)*100., cos(time*10.)*100.)) / resolution.xy;
    //base
    vec2 uv = (gl_FragCoord.xy - vec2(0.5,0.5)) / resolution.xy;
    gl_FragColor = mix(texture2D(prev_layer, uv), texture2D(backbuffer, uv), del);
}

