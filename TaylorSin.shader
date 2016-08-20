float taylor_sin(float x) {
    float x2 = x * x;    
    float xi = x2 * x;   
    
    x -= xi * 0.16666666666;
    xi *= x2; x += xi * 0.00833333333;
    xi *= x2; x -= xi * 0.00019841269;
    xi *= x2; x += xi * 0.00000275573;
    xi *= x2; x -= xi * 2.50521084e-8;
    xi *= x2; x += xi * 1.6059044e-10;
    xi *= x2; x -= xi * 7.6471635e-13;
    xi *= x2; x += xi * 2.8114572e-15;
    
    return x;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {    
	vec2 uv = fragCoord.xy / iResolution.xy;
    uv = (uv * 2.0 - 1.0);
    uv.x *= iResolution.x / iResolution.y;
    
    vec2 fuv = uv * 10.0;
    float lum;
    if(uv.x < 0.0) lum = abs(sin(fuv.x) * sin(fuv.y));
    else lum = abs(taylor_sin(fuv.x) * taylor_sin(fuv.y));
   
    fragColor = vec4(lum,lum,lum,1.0);
}
