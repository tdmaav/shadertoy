const int NUM_SHAPES = 29;
const float R = 1.0;
const float R_MAGNETIC = R * 0.85;
const float R_POLE = R_MAGNETIC * 0.1;
const vec2 pole = vec2(0.0, R_MAGNETIC*0.5);
const float SHAPE_TO_ANGLE = 3.1415 / float(NUM_SHAPES+1);
const float EPSILON = 1e-6;
#define EPSILON_SMOOTH (3.0 / iResolution.x)

float cross2(vec2 a, vec2 b) {
    return a.x * b.y - a.y * b.x;
}
float triangle(vec2 p, vec2 a, vec2 b, vec2 c) {
    float ret;
    ret = cross2(normalize(b-a),p-a);
    ret = max(ret,cross2(normalize(c-b),p-b));
    ret = max(ret,cross2(normalize(a-c),p-c));
    return clamp(smoothstep(0.0,EPSILON_SMOOTH,-ret),0.0,1.0);
}
float circle(vec2 p, vec2 c, float r) {
    float ret = 1.0 - length(p-c) + r;
    return clamp(smoothstep(1.0-EPSILON_SMOOTH,1.0,ret),0.0,1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 uv = fragCoord.xy / iResolution.xy * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y;
     
    // initial circle   
    float sum = circle(uv,vec2(0.0,0.0),R);
        
    // magnetic field shape
    for(int i = 1; i <= NUM_SHAPES; i++) {
        float a = float(i) * SHAPE_TO_ANGLE;
        vec2 c = vec2(sin(a),cos(a));
        sum += triangle(uv, pole,-pole,-c*R_MAGNETIC);
        sum += triangle(uv,-pole, pole, c*R_MAGNETIC);
    }   
        
    // checkerize
    float f = fract(max(sum-EPSILON,0.0));    
    float c = mod(floor(max(sum-EPSILON,0.0)),2.0);    
    sum = mix(f,1.0-f,step(0.1,c));
        
    // poles
    sum = mix(sum,1.0,circle(uv,pole,R_POLE)+circle(uv,-pole,R_POLE));
    
    // color
    vec3 col = mix(vec3(0.475,0.482,0.263),vec3(0.9453125,0.953125,0.8125),sum);    
    
    float noize = texture(iChannel0,uv*0.5).x + texture(iChannel0,uv*0.704).x;
    col *= 1.0 - noize * (1.0-col.x) * 0.5;    
    noize = texture(iChannel0,vec2(atan(uv.y/uv.x) * 0.1,length(uv))).x;
    col *= 1.0 - noize * col.x * 0.2;
    
    // print it to field!
    fragColor = vec4(col,1.0);
}
