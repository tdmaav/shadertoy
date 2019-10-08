/*
 * "Volume Box" by Alexander Alekseev aka TDM - 2019
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 */

const int SLICES = 4;
const float DENSITY = 3.0;

float integrationFunc(float x, float a) {
    return x * 0.5 - cos(x * a) / (2.0 * a);
}

float functionMean(float a, float b, float f) {
    a = -a * 0.5 + 0.5;
    b = -b * 0.5 + 0.5;
    float Fa = integrationFunc(a,f);
    float Fb = integrationFunc(b,f);
    return (Fb - Fa) / (b - a);
}

// main
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 uv = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;   
    vec2 mouse = iMouse.xy / iResolution.xy * 4.0 - 2.0;
        
    // ray   
    vec3 ang;
    if(iMouse.z > 0.0) {
        ang = vec3(0.0,-mouse.y,mouse.x);
    } else {
        ang = vec3(sin(iTime*0.4)*2.0,0.0,cos(iTime*0.35)*3.0);
    }
	mat3 rot = fromEuler(ang);    
    vec3 ori = vec3(0.0,0.0,4.5) * rot;
    vec3 dir = normalize(vec3(uv.xy,-3.0)) * rot;      
             
    // color
    vec3 p, rp0, rp1;
    vec3 rcolor = vec3(0.0);
    if(intersectionRayBox(ori,dir,vec3(1.0),rp0,rp1)) {       
        for(int i = 0; i < SLICES; i++) {
            vec3 r0 = rp0;
            vec3 r1 = rp1;
            
            r0 = mix(r0,r1,float(SLICES-i-1)/float(SLICES));
            r1 = r0 + (rp1-rp0)/float(SLICES);
            r0 += 1.7;
            r1 += 1.7;

            // modulate color
            float fm = 0.6;
            vec3 color;
            color.x = functionMean(r0.x,r1.x,7.0*fm);
            color.y = functionMean(r0.x,r1.x,11.0*fm);
            color.z = functionMean(r0.x,r1.x,13.0*fm);
            color.yz *= functionMean(r0.x,r1.x,27.0*fm);
            color = 1.0 - color;

            color.z *= functionMean(r0.y,r1.y,11.0*fm);
            color.y *= functionMean(r0.y,r1.y,13.0*fm);
            color.x *= functionMean(r0.y,r1.y,17.0*fm);
            color = 1.0 - color;

            color.z *= functionMean(r0.z,r1.z,5.0*fm);
            color.y *= functionMean(r0.z,r1.z,7.0*fm);
            color.x *= functionMean(r0.z,r1.z,11.0*fm);
            color = 1.0 - color;
            color = pow(color,vec3(8.0));
            color += 0.02;
            color = log(1.0+color*1.5);

            // modulate density
            float d = length(r1 - r0) * functionMean(r0.z+r0.y,r1.z+r1.y,22.0);            
            rcolor = mix(rcolor, color, clamp(log(1.0+d*DENSITY),0.0,1.0));
        }
        
        // sRGB
        rcolor = pow(rcolor, vec3(1.0/2.2));
        rcolor = saturation(rcolor, 2.0);
    }
               
	fragColor = vec4(rcolor,1.0);
}
