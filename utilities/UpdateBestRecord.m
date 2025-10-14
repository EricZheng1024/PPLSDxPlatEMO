function [solB, flag_updated] = UpdateBestRecord(solPrime, solB, z, W, W_, rho)
g_new = g_cn(solPrime.obj, z, W, W_, rho);
g_old = g_cn(solB.objs, z, W, W_, rho);
flag_updated = g_new>g_old;
solB(flag_updated) = solPrime;
end