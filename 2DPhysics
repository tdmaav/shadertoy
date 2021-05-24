/*
 * "2D Physics (balls)" by Alexander Alekseev aka TDM - 2021
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com 
 */
 
const int NUM_OBJECTS = 20;
const float BALL_SIZE = 0.15;
const float MAX_VELOCITY = 4.0;
const float ELASTICITY = 0.5;
const vec2 FRAME_SIZE = vec2(1.2,0.8);
const vec2 GRAVITY = vec2(0.0,-1.0);

const float PI = 3.141592;
const float DEG2RAD = PI / 180.0;

/*
 * math
 */
 
vec3 hash3( uint n ) {
	n = (n << 13U) ^ n;
    n = n * (n * n * 15731U + 789221U) + 1376312589U;
    uvec3 k = n * uvec3(n,n*16807U,n*48271U);
    return vec3( k & uvec3(0x7fffffffU))/float(0x7fffffff);
}
float cross2(vec2 a, vec2 b) {
    return a.x * b.y - a.y * b.x;
}
vec2 cross2(vec2 a, float b) {
    return vec2(a.y * b, a.x * -b);
}
vec2 cross2(float a, vec2 b) {
    return vec2(b.y * -a, b.x * a);
}
vec2 rotateZ(vec2 v, float a) {
    lowp vec2 sc = vec2(sin(a),cos(a));
    lowp vec2 ret = v;
    ret.x = v.x * sc.y - v.y * sc.x;
    ret.y = v.x * sc.x + v.y * sc.y;
    return ret;
}

/*
 * body
 */

struct Body {
    vec2 pos;
    vec2 vel;
    float ang;
    float ang_vel;
    float inv_mass;
    float inv_momentum;
};

Body getBody(sampler2D buf, vec2 ires, int i) {
    vec4 data0 = texture(buf, (vec2(float(i),0.0) + 0.5) * ires);
    vec4 data1 = texture(buf, (vec2(float(i),1.0) + 0.5) * ires);
    
    Body body;
    body.pos = data0.xy;
    body.vel = data0.zw;
    body.ang = data1.x;
    body.ang_vel = data1.y;
    body.inv_mass = data1.z;
    body.inv_momentum = data1.w;
    return body;
}

void initBody(int id, inout Body body) {
    vec3 rnd = hash3(uint(id));
    body.pos = (rnd.xy * 2.0 - 1.0) * 0.4;
    body.vel = cross2(1.0,body.pos) * 2.0;
    body.ang_vel = length(body.pos) * -8.0;
    body.inv_mass = 1.0;
    body.inv_momentum = 1.0 / (0.5 * (1.0/body.inv_mass) * BALL_SIZE * BALL_SIZE);
}

/*
 * solver
 */

vec2 collisionWithPlane(inout Body b0, vec3 plane) {
    vec2 displace = vec2(0.0);
    vec2 normal = normalize(plane.xy);
    float dist = dot(b0.pos,normal) + plane.z;
    float penetration = BALL_SIZE - dist;
    if(penetration > 0.0) {
        displace += normal * penetration;

        vec2 r0 = -normal * BALL_SIZE;        

        // normal
        vec2 vel0 = b0.vel + cross2(b0.ang_vel,r0);
        vec2 rel_vel = vel0;  
        
        float w1 = cross2(r0,normal);

        float a = (1.0 + ELASTICITY) * dot(normal,rel_vel);
        float b = b0.inv_mass + w1 * w1 * b0.inv_momentum;
        float lambda = max(-a / b, 0.0);

        b0.vel += normal * (lambda * b0.inv_mass);
        b0.ang_vel += cross2(r0, normal) * lambda * b0.inv_momentum;

        // friction
        vel0 = b0.vel + cross2(b0.ang_vel,r0);
        rel_vel = vel0;  

        vec2 tangent = cross2(normal,1.0);
        w1 = cross2(r0,tangent);

        a = (1.0 + ELASTICITY) * dot(tangent,rel_vel);
        b = b0.inv_mass + w1 * w1 * b0.inv_momentum;
        float lambdaF = clamp(-a / b, -lambda, lambda);

        b0.vel += tangent * (lambdaF * b0.inv_mass);
        b0.ang_vel += cross2(r0, tangent) * lambdaF * b0.inv_momentum;
    }
    return displace;
}

vec2 collisionWithBody(inout Body b0, in Body b1) {
    vec2 displace = vec2(0.0);
    vec2 normal = b0.pos - b1.pos;
    float dist = length(normal);
    float penetration = 2.0 * BALL_SIZE - dist;
    if(penetration > 0.0) {
        normal /= dist;
        displace += normal * penetration * 0.5;

        vec2 r0 = -normal * BALL_SIZE;
        vec2 r1 = normal * BALL_SIZE;
        
        // normal
        vec2 vel0 = b0.vel + cross2(b0.ang_vel,r0);
        vec2 vel1 = b1.vel + cross2(b1.ang_vel,r1);
        vec2 rel_vel = vel0 - vel1;
        
        float w1 = cross2(r0,normal);
        float w2 = cross2(r1,normal);

        float a = (1.0 + ELASTICITY) * dot(normal,rel_vel);
        float b = b0.inv_mass + b1.inv_mass +
            w1 * w1 * b0.inv_momentum +
            w2 * w2 * b1.inv_momentum;
        float lambda = max(-a / b, 0.0);

        b0.vel += normal * (lambda * b0.inv_mass);
        b0.ang_vel += cross2(r0, normal) * lambda * b0.inv_momentum;
        b1.vel -= normal * (lambda * b1.inv_mass);
        b1.ang_vel -= cross2(r1, normal) * lambda * b1.inv_momentum;

        // friction
        vel0 = b0.vel + cross2(b0.ang_vel,r0);
        vel1 = b1.vel + cross2(b1.ang_vel,r1);
        rel_vel = vel0 - vel1;  

        vec2 tangent = cross2(normal,1.0);
        w1 = cross2(r0,tangent);
        w2 = cross2(r1,tangent);

        a = (1.0 + ELASTICITY) * dot(tangent,rel_vel);
        b = b0.inv_mass + b1.inv_mass +
            w1 * w1 * b0.inv_momentum +
            w2 * w2 * b1.inv_momentum;
        float lambdaF = clamp(-a / b, -lambda, lambda);

        b0.vel += tangent * (lambdaF * b0.inv_mass);
        b0.ang_vel += cross2(r0, tangent) * lambdaF * b0.inv_momentum;
    }
    return displace;
}

void solve(sampler2D data, inout Body b0, int id, vec2 ires) {
    vec2 displace = vec2(0.0);
    
    // collision detection
    for(int i = 0; i < NUM_OBJECTS; i++) {
        if(i == id) continue;
        
        Body b1 = getBody(data, ires, i);
        displace += collisionWithBody(b0,b1);
    }
    
    // walls
    displace += collisionWithPlane(b0, vec3(0.0,1.0,FRAME_SIZE.y));
    displace += collisionWithPlane(b0, vec3(0.0,-1.0,FRAME_SIZE.y));
    displace += collisionWithPlane(b0, vec3(1.0,0.0,FRAME_SIZE.x));
    displace += collisionWithPlane(b0, vec3(-1.0,-.0,FRAME_SIZE.x));

    b0.pos += displace;
}

/*
 * store
 */

float isInside( vec2 p, vec2 c ) { 
    vec2 d = abs(p-0.5-c) - 0.5;
    return -max(d.x,d.y);
}
void storeBody(in int id, in Body b, inout vec4 col, in vec2 uv) {
    col = isInside(uv,vec2(float(id),0)) > 0.0 ? vec4(b.pos,b.vel) : col;
    col = isInside(uv,vec2(float(id),1)) > 0.0 ? vec4(b.ang,b.ang_vel,b.inv_mass,b.inv_momentum) : col;
}

/*
 * Dynamics
 */
 
vec2 getForce(vec2 x, vec2 v) {
    vec2 force = vec2(0.0);
    
    if(iMouse.z > 0.5) {
        vec2 mouse = iMouse.xy / iResolution.xy * 2.0 - 1.0;
        mouse.x *= iResolution.x / iResolution.y;
        vec2 dir = x.xy - mouse;
        float p = length(dir);        
        force += 5.0 * normalize(dir) / p;
    }
    
    return force;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    int id = int(fragCoord.x);
    if(id >= NUM_OBJECTS) discard;
    vec2 ires = 1.0 / iChannelResolution[0].xy;
    
    // load    
    Body body = getBody(iChannel0, ires, id);
    if(iFrame == 0) initBody(id, body); // init
    
    float dt = min(iTimeDelta, 0.07);
   
    // semi-implicit Euler
    
    // integrate forces
    vec2 force = getForce(body.pos, body.vel);
    body.vel += (force * body.inv_mass + GRAVITY) * dt;
    
    // limit max velocity
    float len2 = dot(body.vel,body.vel);
    if(len2 > MAX_VELOCITY * MAX_VELOCITY)
        body.vel *= inversesqrt(len2) * MAX_VELOCITY;
    
    // integrate velocity
    body.pos += body.vel * dt;
    body.ang += body.ang_vel * dt;
    
    // store
    fragColor = vec4(0.0);
    storeBody(id, body, fragColor, fragCoord);
}

/*
 * Collision solver (1st iteration)
 */

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    int id = int(fragCoord.x);
    if(id >= NUM_OBJECTS) discard;
    vec2 ires = 1.0 / iChannelResolution[0].xy;
    
    // solve collisions    
    Body body = getBody(iChannel0, ires, id);
    solve(iChannel0,body,id,ires);
    
    // store
    fragColor = vec4(0.0);
    storeBody(id, body, fragColor, fragCoord);
}

/*
 * "2D Physics (balls)" by Alexander Alekseev aka TDM - 2021
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com
 *
 * Render 
 */
 
#define PIX length(fwidth(p))
const vec3[] COLORS = vec3[4] (
    vec3(234.,67.,53.) / 255.,
    vec3(66, 133, 244) / 255.,
    vec3(251, 188, 5) / 255.,
    vec3(52, 168, 83) / 255.
);

float circle(vec2 p, vec2 c, float w) {
    float dist = length(p - c) - w;
    return smoothstep(PIX,0.0,dist);
}

float frame(vec2 p, vec2 size, float w) {
    const float SMOOTH = 0.2;
    size -= SMOOTH;
	p = abs(p)-size;
    float dist = length(p-min(p,0.0)) - SMOOTH;
    float shad = 1.0 - dist * 2.0;
    shad = 1.0 - shad * shad * shad;
    shad = 1.0 - (1.0 - shad) * smoothstep(0.0,PIX,dist);
    return shad * 0.1 + 0.9;
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
    for(int i = 0; i < NUM_OBJECTS; i++) {
        Body body = getBody(iChannel0, ires, i);
        float ba = circle(uv,body.pos,BALL_SIZE*0.98);
        ba *= 1.0-circle(uv,body.pos,BALL_SIZE*0.3);
                  
        for(int j = 0; j < 5; j++) {
            float ang = body.ang + float(j) * (360./5.) * DEG2RAD;
            vec2 o = rotateZ(vec2(0.0,BALL_SIZE*1.25), ang);
            ba *= 1.0 - circle(uv, body.pos + o, BALL_SIZE * 0.4);
        }  
        c = mix(c,COLORS[i%4],ba);
    }
    c *= frame(uv,vec2(FRAME_SIZE*1.08),0.01);
    
    // final
	fragColor = vec4(c,1.0);
}
