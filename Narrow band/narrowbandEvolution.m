function phi=narrowbandEvolution(phi,g,front,band,delta_t,bandIterNum)

for n=1:bandIterNum
    speed = calcspeed( phi, g, front );%计算边界点上的速度值
    speed = extendspeedband( speed, front, band );%将边界上点的速度值扩展到整个窄带域
    phi = phi - delta_t .* speed;
end