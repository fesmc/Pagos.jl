Asp, u0, b = calc_matrix_ssa(ux,uy,N,N_ab,taud_acx,taud_acy,β_acx,β_acy,dx)
A == Asp
lsd.u0 == u0
lsd.b == b

@b calc_matrix_ssa($ux,$uy,$N,$N_ab,$taud_acx,$taud_acy,$β_acx,$β_acy,$dx)
