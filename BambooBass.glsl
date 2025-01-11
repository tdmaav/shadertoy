/*
 * "Bamboo bass" by Alexander Alekseev aka TDM - 2025
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com
 */

const float PI = 3.141592;

vec2 hash21(float p) {
    return fract(sin(p * vec2(110.082, 175.025)) * vec2(19871.8972, 11571.72));
}

float noteToFreq(float n) {
	return 440.0 * exp2((n - 69.0) / 12.0);
}

float tx_wave(float x) {    
    float y = sin(x * PI);
    return y * y * -sign(mod(x,2.0)-1.0);
}

float bamboo(float t, float h) {
    float ph1 = t * h;
    float ph2 = ph1 * 2.0;
    float mod1 = pow(max(1.0 - t * 3.0, 0.0), 2.0) * 3.0;
    float mod2 = max(1.0 - t * 1000.0, 0.0) * 80.0; // click
    float wave = tx_wave(ph1 + tx_wave(ph2) * (mod1 + mod2));
    float amp = smoothstep(0.0, 1.0, 1.0 - t);
    return wave * amp * step(0.0,t);
}

// from https://www.shadertoy.com/view/XldcRr
vec2 bamboo_reverb(float t, float h) {
    const float invRate = 1./1000.;
    const float dur = 2.;    
    vec2 reverb = vec2(0.0);
    vec2 sum = vec2(0.);
    for ( float i = 0.; i < dur; i += invRate ) {
        vec2 rand = hash21(i);
        vec2 t2 = i + invRate*rand*5.0;
        vec2 amp = exp2(-t2*3.0);
        reverb.x += bamboo(t - t2.x, h ) * amp.x;
        reverb.y += bamboo(t - t2.y, h ) * amp.y;
        sum += amp;
    }
    reverb /= sum;
    return bamboo(t,h) + reverb * 12.0;
}

vec2 mainSound( in int samp,float time) {
    float note = mod(floor(time), 16.0);
    float btime = fract(time);
    return bamboo_reverb(btime, noteToFreq(note + 40.0));
}

#define EPSILON_SMOOTH(p) length(fwidth(p))

vec2 rotateZ(vec2 v, float a) {
    vec2 s = vec2(sin(a), cos(a));
    return vec2(v.x*s.y - v.y*s.x, v.x*s.x + v.y*s.y);
}

float circle(vec2 p, vec2 c, float r) {
    float ret = 1.0 - length(p-c) + r;
    return clamp(smoothstep(1.0-EPSILON_SMOOTH(p),1.0,ret),0.0,1.0);
}

vec3 mainImage(vec2 uv, float t) {
    float woof = sin(t*80.0) * pow(1.0 - fract(t),2.0);
    uv = rotateZ(uv,0.7);
    uv.x *= 1.2;
    uv *= 1.0 + woof * 0.02;
    
    vec2 dif = vec2(0.15 + woof * 0.08,0.0);
    float r = circle(uv,vec2(0.),1.0);
    r = max(r, circle(uv, vec2(0.8,0.0), 0.5));
    r = min(r, 1.0-circle(uv,vec2(0.),0.95));
    r = max(r, circle(uv, vec2(woof*0.03,0.0),0.9));
    r = min(r, 1.0-circle(uv, dif, 0.42));
    r = max(r, circle(uv,dif, 0.35)*(1.0-circle(uv,dif+vec2(0.15,0.0),0.4)));
    
	return mix(vec3(0.0,0.2,0.3),vec3(1.0,0.1,0.8),r);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
	vec2 uv = fragCoord.xy / iResolution.xy * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y;
    
    float n = texture(iChannel0, fragCoord / 1024.0).x;
    vec3 col = mainImage(uv,iTime - iTimeDelta * n * 3.0);
	fragColor = vec4(col,1.0);
}
