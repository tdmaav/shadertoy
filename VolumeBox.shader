/*
 * "Volume box" by Alexander Alekseev aka TDM - 2017
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com
 */

mat3 fromEuler(vec3 ang) {
    vec2 a1 = vec2(sin(ang.x),cos(ang.x));
    vec2 a2 = vec2(sin(ang.y),cos(ang.y));
    vec2 a3 = vec2(sin(ang.z),cos(ang.z));
    mat3 m;
    m[0] = vec3(a1.y*a3.y+a1.x*a2.x*a3.x,a1.y*a2.x*a3.x+a3.y*a1.x,-a2.y*a3.x);
    m[1] = vec3(-a2.y*a1.x,a1.y*a2.y,a2.x);
    m[2] = vec3(a3.y*a1.x*a2.x+a1.y*a3.x,a1.x*a3.x-a1.y*a3.y*a2.x,a2.y*a3.y);
    return m;
}
bool intersectionRayBox(vec3 o, vec3 d, vec3 ext, out vec3 r0, out vec3 r1) {
    vec3 t0 = (-o - ext) / d; 
    vec3 t1 = (-o + ext) / d;    
    vec3 n = min(t0,t1); n.x = max(max(n.x,n.y),n.z);
    vec3 f = max(t0,t1); f.x = min(min(f.x,f.y),f.z);
    r0 = o + d * n.x;
    r1 = o + d * f.x;
    return bool(step(n.x,f.x));
}

float integrationFunc(float x) {
    //return 2.0 * pow(max(x,1e-6),1.5) / 3.0;
    //return x*x/2.0;
    //return x*x*x/3.0;
    return x*0.5-cos(x * 30.0) / (2.0 * 30.0);
}

float functionMean(float a, float b) {
    a = -a * 0.5 + 0.5;
    b = -b * 0.5 + 0.5;
    float Fa = integrationFunc(a);
    float Fb = integrationFunc(b);
    return (Fb - Fa) / (b - a);
}

// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
vec3 ACESFilm( vec3 x ) {
    float a = 2.51;
    float b = 0.03;
    float c = 2.43;
    float d = 0.59;
    float e = 0.14;
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), vec3(0.0),vec3(1.0));
}

// main
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 iuv = fragCoord.xy / iResolution.xy * 2.0 - 1.0;
    vec2 uv = iuv;
    uv.x *= iResolution.x / iResolution.y;    
    vec2 mouse = iMouse.xy / iResolution.xy * 4.0 - 2.0;
        
    // ray   
    vec3 ang;
    if(iMouse.z > 0.0) {
        ang = vec3(0.0,-mouse.y,mouse.x);
    } else {
        ang = vec3(0.0,sin(iGlobalTime)*0.75,cos(iGlobalTime*1.5)*0.75);
    }
    mat3 rot = fromEuler(ang);    
    vec3 ori = vec3(0.0,0.0,5.0);
    vec3 dir = normalize(vec3(uv.xy,-3.0));      
    ori = ori * rot;
    dir = dir * rot;
             
    // color
    vec3 p, r0, r1;
    vec3 color = vec3(0.0);
    if(intersectionRayBox(ori,dir,vec3(1.0),r0,r1)) {
        
        // modulate density
        float density = length(r1 - r0) * functionMean(r0.y,r1.y);
        color = vec3(density);
        
        // modulate color
        color.yz *= functionMean(r0.x,r1.x);
        
        // tonemapping
        color = ACESFilm(color);
    }
               
    fragColor = vec4(color,1.0);
}
