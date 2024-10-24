/*
 * "2D Physics (analytic springs)" by Alexander Alekseev aka TDM - 2023
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com 
 */
 
const int NUM_OBJECTS = 64;
const float OBJECTS_GAP = 0.5 * (8.0 / float(NUM_OBJECTS));
const float BALL_SIZE = OBJECTS_GAP / 2.0;

/*
 * Dynamics
 */
 
const float k = 100.0;
const float damp = 0.45;

vec2 springSystemModel(float t, vec2 b0, vec2 b1, vec2 b2) {
    float rh = sqrt(3.*k);
    float rp = sqrt(k);
    
    vec2 CH = vec2(b0.x - (b1.x+b2.x)/2., 
                  (b0.y - (b1.y+b2.y)/2.) / rh);
    vec2 CP = vec2((b1.x+b2.x)/2.,
                   (b1.y+b2.y)/(2.*rp));
                   
    vec2 exph = vec2(cos(rh*t), sin(rh*t));
    vec2 expp = vec2(cos(rp*t), sin(rp*t));
       
    float x  = dot(CH, exph) + 
               dot(CP, expp);
    float dx = dot(CH, rh*exph.yx*vec2(-1.0,1.0)) + 
               dot(CP, rp*expp.yx*vec2(-1.0,1.0));
    
    return vec2(x,dx) * exp(-t * damp); // new state
}

/*
 * body
 */

vec2 getBody(sampler2D buf, vec2 ires, int i) {
    return texture(buf, (vec2(float(i),0.0) + 0.5) * ires).xy;
}

void initBody(int id, inout vec2 body) {
    body = vec2(0.0);
}

/*
 * store
 */

float isInside( vec2 p, vec2 c ) { 
    vec2 d = abs(p-0.5-c) - 0.5;
    return -max(d.x,d.y);
}
void storeBody(in int id, in vec2 b, inout vec4 col, in vec2 uv) {
    col = isInside(uv,vec2(float(id),0)) > 0.0 ? vec4(b,0.0,0.0) : col;
}

/**
 * coords
 */
 
vec2 toScreenspace(int id, float y) {
    float x = -float(NUM_OBJECTS-1) * 0.5 * OBJECTS_GAP + 
    float(id) * OBJECTS_GAP;
    return vec2(x,y);
}

/*
 * "2D Physics (analytic springs)" by Alexander Alekseev aka TDM - 2023
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com 
 */

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    int id = int(fragCoord.x);
    if(id >= NUM_OBJECTS) discard;
    vec2 ires = 1.0 / iChannelResolution[0].xy;
    
    // load    
    vec2 body = getBody(iChannel0, ires, id);
    if(iFrame == 0) {
        initBody(id, body); // init
    } else {
    
        float dt = min(iTimeDelta, MAX_DT);

        // mouse impulse
        if(iMouse.z > 0.5) {
            vec2 mouse = iMouse.xy / iResolution.xy * 2.0 - 1.0;
            mouse.x *= iResolution.x / iResolution.y;
            vec2 pos_ss = toScreenspace(id,body.x);
            vec2 dir = pos_ss.xy - mouse;

            float t = abs(dir.x);
            body.y += -k * 0.5 * dt * smoothstep(0.2,0.0,t);
        }

        // spring
        vec2 bl = getBody(iChannel0, ires, id-1);
        vec2 br = getBody(iChannel0, ires, id+1);
        body = springSystemModel(dt, body, bl, br);
    }
    
    // store
    storeBody(id, body, fragColor, fragCoord);
}

/*
 * "2D Physics (analytic springs)" by Alexander Alekseev aka TDM - 2023
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com 
 */

const vec3 COLOR = vec3(72, 143, 240) / 255.;

float line(vec2 p, vec2 p0, vec2 p1, float w) {
    vec2 d = p1 - p0;
    float t = clamp(dot(d,p-p0) / dot(d,d), 0.0,1.0);
    vec2 proj = p0 + d * t;
    float dist = length(p - proj);
    return step(dist-w,0.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
	vec2 uv = fragCoord.xy / iResolution.xy * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y;
    
    vec2 mouse = iMouse.xy / iResolution.xy * 2.0 - 1.0;
    mouse.x *= iResolution.x / iResolution.y;
    mouse.y = 0.0;
    
    vec3 c = vec3(1.0);
    vec2 ires = 1.0 / iChannelResolution[0].xy;
        
    // objects
    for(int i = 0; i < NUM_OBJECTS-1; i++) {
        vec2 body = getBody(iChannel0, ires, i);
        vec2 b1 = getBody(iChannel0, ires, i+1);
        vec2 pos_ss0 = toScreenspace(i,body.x);
        vec2 pos_ss1 = toScreenspace(i+1,b1.x);
        
        float ba = line(uv,pos_ss0,pos_ss1,BALL_SIZE);
        c = mix(c,COLOR,ba);
        
        vec2 mid = (pos_ss0 + pos_ss1) * 0.5;
        ba = line(uv,
            mid,
            mid+vec2(0.0,-2.0),
            BALL_SIZE);
        c = mix(c,COLOR,ba * 0.2 - min(uv.y,0.0)*0.005);
    }
    
    // final
	fragColor = vec4(c,1.0);
}
