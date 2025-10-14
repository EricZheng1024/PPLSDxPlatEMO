function g = g_cn(objs, z, W, W_, rho)
g = rho*sum((z-objs).*W, 2) + (1-rho)*min((z-objs).*W_, [], 2);
end