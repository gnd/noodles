// "Volumetric explosion" by Duke
// https://www.shadertoy.com/view/lsySzd
//-------------------------------------------------------------------------------------
// Based on "Supernova remnant" (https://www.shadertoy.com/view/MdKXzc)
// and other previous shaders
// otaviogood's "Alien Beacon" (https://www.shadertoy.com/view/ld2SzK)
// and Shane's "Cheap Cloud Flythrough" (https://www.shadertoy.com/view/Xsc3R4) shaders
// Some ideas came from other shaders from this wonderful site
// Press 1-2-3 to zoom in and zoom out.
// License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
//-------------------------------------------------------------------------------------
#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define eps 0.00001
#define MAXSTEPS 128
#define MAXDIST 12
#define pi 3.14159265
#define R(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)

// Random number generator
float rnd(vec2 co){
        // implementation found at: lumina.sourceforge.net/Tutorials/Noise.html
        float k = sin(dot(co.xy ,vec2(12.9898,78.233)));
        return fract(k + k);
}

float noise( in vec3 p )
{
        vec2 ip = vec2(floor(p.xy));
        vec2 u = fract(p.xy);
        // http://www.iquilezles.org/www/articles/morenoise/morenoise.htm
        //u = u*u;
        u = u*u*u*((6.0*u-15.0)*u+10.0);

        float res = mix(
                                        mix(rnd(ip), rnd(ip+vec2(1.0,0.0)), u.x),
                                        mix(rnd(ip+vec2(0.0,1.0)), rnd(ip+vec2(1.0,1.0)), u.x),
                                        u.y);

        //return res - 0.25;
        return 2.0* (res-0.5);
}

float fbm( vec3 p ) {
   return noise(p*.06125)*.5 + noise(p*.125)*.25 + noise(p*.25)*.125 + noise(p*.4)*.2;
}

float Sphere( vec3 p, float r ) {
    return length(p)-r;
}

//==============================================================
// otaviogood's noise from https://www.shadertoy.com/view/ld2SzK
//--------------------------------------------------------------
// This spiral noise works by successively adding and rotating sin waves while increasing frequency.
// It should work the same on all computers since it's not based on a hash function like some other noises.
// It can be much faster than other noise functions if you're ok with some repetition.
const float nudge = 4.;	// size of perpendicular vector
float normalizer = 1.0 / sqrt(1.0 + nudge*nudge);	// pythagorean theorem on that perpendicular to maintain scale
float SpiralNoiseC(vec3 p) {
    float n = -mod(time * 0.2,-2.); // noise amount
    float iter = 2.0;
    for (int i = 0; i < 8; i++)
    {
        // add sin and cos scaled inverse with the frequency
        n += -abs(sin(p.y*iter) + cos(p.x*iter)) / iter;	// abs for a ridged look
        // rotate by adding perpendicular and scaling down
        p.xy += vec2(p.y, -p.x) * nudge;
        p.xy *= normalizer;
        // rotate on other axis
        p.xz += vec2(p.z, -p.x) * nudge;
        p.xz *= normalizer;
        // increase the frequency
        iter *= 1.733733;
    }
    return n;
}

float VolumetricExplosion(vec3 p) {
    float final = Sphere(p,4.);
    final += fbm(p*50.);
    final += SpiralNoiseC(p.zxy*0.4132+333.)*3.0; //1.25;

    return final;
}

float map(vec3 p) {
	float VolExplosion = VolumetricExplosion(p/0.5)*0.5; // scale
	return VolExplosion;
}

// assign color to the media
vec3 computeColor( float density, float radius ) {
	// color based on density alone, gives impression of occlusion within
	// the media
	vec3 result = mix( vec3(1.0,0.9,0.8), vec3(0.4,0.15,0.1), density );

	// color added to the media
	vec3 colCenter = 7.*vec3(0.8,1.0,1.0);
	vec3 colEdge = 1.5*vec3(0.48,0.53,0.5);
	result *= mix( colCenter, colEdge, min( (radius+.05)/.9, 1.15 ) );

	return result;
}

bool RaySphereIntersect(vec3 org, vec3 dir, out float near, out float far) {
	float b = dot(dir, org);
	float c = dot(org, org) - 8.;
	float delta = b*b - c;
	if( delta < 0.0)
		return false;
	float deltasqrt = sqrt(delta);
	near = -b - deltasqrt;
	far = -b + deltasqrt;
	return far > 0.0;
}

void main()
{
    float key = -1.;

    vec2 uv = gl_FragCoord.xy/resolution.xy;

	// ro: ray origin
	// rd: direction of the ray
	vec3 rd = normalize(vec3((gl_FragCoord.xy-0.5*resolution.xy)/resolution.y, 1.));
	vec3 ro = vec3(0., 0., -6.+key*1.6);

	// ld, td: local, total density
	// w: weighting factor
	float ld=0., td=0., w=0.;

	// t: length of the ray
	// d: distance function
	float d=1., t=0.;

    const float h = 0.1;

	vec4 sum = vec4(0.0);

    float min_dist=0.0, max_dist=0.0;

    if (RaySphereIntersect(ro, rd, min_dist, max_dist)) {
		t = min_dist*step(t,min_dist);

		// raymarch loop
    	for (int i=0; i<86; i++) {
			vec3 pos = ro + t*rd;

			// Loop break conditions.
	    	if(td>0.9 || d<0.12*t || t>10. || sum.a > 0.99 || t>max_dist) break;

        	// evaluate distance function
        	float d = map(pos);

        	//d = uv.x < 0.5 ? abs(d)+0.07 : d;
        	d = abs(d)+0.0007;
       
			// change this string to control density
			d = max(d,0.003);

        	// point light calculations
        	vec3 ldst = vec3(0.0)-pos;
        	float lDist = max(length(ldst), 0.001);

        	// the color of light
        	vec3 lightColor=vec3(1.0,0.5,0.25);

        	sum.rgb+=(lightColor/exp(lDist*lDist*lDist*.08)/30.); // bloom

			if (d<h) {
				// compute local density
				ld = h - d;
	            // compute weighting factor
				w = (1. - td) * ld;
				// accumulate density
				td += w + 1./200.;
				vec4 col = vec4( computeColor(td,lDist), td );
	            // emission
    	        sum += sum.a * vec4(sum.rgb, 0.0) * 0.2 / lDist;
				// uniform scale density
				col.a *= 0.2;
				// colour by alpha
				col.rgb *= col.a;
				// alpha blend in contribution
				sum = sum + col*(1.0 - sum.a);
			}

		td += 1./70.;

        // trying to optimize step size
        t += max(d * 0.08 * max(min(length(ldst),d),2.0), 0.01);
 
	}

    // simple scattering
    sum *= 1. / exp( ld * 0.2 ) * 0.8;

   	sum = clamp( sum, 0.0, 1.0 );

    sum.xyz = sum.xyz*sum.xyz*(3.0-2.0*sum.xyz);

	}

    PixelColor = vec4(sum.xyz,1.0);
}
