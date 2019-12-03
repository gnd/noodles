// sony height mapping
vec3 texcol(vec3 p) {
    vec3 image = texture2D(sony, (vec2(9.-p.x,4.-p.z))*vec2(.05,.1)).xyz;
    return image;
}
