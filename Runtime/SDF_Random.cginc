float hash11(float co) { return frac(sin(co * (91.3458)) * 47453.5453); }

float hash11(uint n)
{
    n = (n << 13U) ^ n;
    n = n * (n * n * 15731U + 789221U) + 1376312589U;
    // floating point conversion from http://iquilezles.org/www/articles/sfrand/sfrand.htm
    return asfloat((n >> 9U) | 0x3f800000U) - 1.0;
}

float3 hash13(float p)
{
    float3 p3 = frac(float3(p,p,p) * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.xxy + p3.yzz) * p3.zyx);
}

float hash31(float3 p3)
{
    p3 = frac(p3 * float3(.1031, .11369, .13787));
    p3 += dot(p3, p3.yzx + 19.19);
    return -1.0 + 2.0 * frac((p3.x + p3.y) * p3.z);
}

float3 hash33(float3 p3)
{
    p3 = frac(p3 * float3(.1031, .11369, .13787));
    p3 += dot(p3, p3.yxz + 19.19);
    return -1.0 + 2.0 * frac(float3((p3.x + p3.y) * p3.z, (p3.x + p3.z) * p3.y, (p3.y + p3.z) * p3.x));
}
