mutable struct LinSolveSSA{T<:AbstractFloat}
    ux_old::Matrix{T}
    uy_old::Matrix{T}
    N::Matrix{T}
    N_ab::Matrix{T}
    u::Vector{T}
    b::Vector{T}
    Ai::Vector{T}
    Aj::Vector{T}
    Av::Vector{T}
    Asp::SparseMatrixCSC{T}
    prob::LinearProblem
    cache::LinearSolve.LinearCache
    prob_is_defined::Bool
end

function LinSolveSSA(T::Type{<:AbstractFloat}, nx::Int, ny::Int)
    n_terms = 9;
    n_u =2*nx*ny;
    n_sprs = n_u*n_terms;

    ux_old = zeros(T, nx, ny)
    uy_old = zeros(T, nx, ny)
    N = zeros(T, nx, ny)
    N_ab = zeros(T, nx, ny)

    # Dense arrays to hold u, b
    u = zeros(T, n_u)
    b = zeros(T, n_u)
    
    # A component vectors (I, J, V)
    Ai = zeros(T, n_sprs)
    Aj = zeros(T, n_sprs)
    Av = zeros(T, n_sprs)

    # Sparse array
    # This array will be reallocated during the first run of CalculateVelocitySSA!
    # to define the sparsity pattern correctly.
    Asp = spzeros(T, n_sprs)
    
    # Initialize LinearProblem and cache
    prob = LinearProblem(Asp, b; u0=u)
    cache = init(prob)

    # Set to false for now, so that we know Asp needs to be defined with
    # correct sparsity pattern and LinearProblem and solution are needed
    prob_is_defined = false

    return LinSolveSSA{T}(ux_old, uy_old, N, N_ab, u, b, Ai, Aj, Av, Asp, prob, cache, prob_is_defined)
end

function CalculateVelocitySSA!(ux,uy,ls,H,μ,taud_acx,taud_acy,β_acx,β_acy,dx)
    # Calculate the diagnostic SSA velocity solution 
    # given ice thickness, viscosity, driving stress and basal friction coefficient

    # Make local references to LinSolveSSA struct matrices
    # (these will be updated in the struct since they are mutable arrays)
    ux_old = ls.ux_old
    uy_old = ls.uy_old
    N = ls.N 
    N_ab = ls.N_ab
    u = ls.u
    b = ls.b
    Ai = ls.Ai
    Aj = ls.Aj
    Av = ls.Av
    Asp = ls.Asp

    # Get some constants 
    nx, ny = size(H);

    dy = dx;
    
    dxdx = (dx*dx);
    dydy = (dy*dy);
    dxdy = (dx*dy);

    # Update previous velocity solution
    ux_old .= ux
    uy_old .= uy

    # Update vertically-integrated viscosity (aa-nodes)
    N .= H .* μ;

    # Stagger N to ab-nodes
    N_ab .= N;

    for i in 1:nx 
        for j in 1:ny 

            # BC: Periodic boundary conditions in x and y 
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1 
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1 
            end
            
            N_ab[i,j] = 0.25 * (N[i,j]+N[ip1,j]+N[i,jp1]+N[ip1,jp1]);

        end 
    end 

    #
    #
    # Populate SSA stress balance matrix equation Ax = b 
    # [ A_ux   A_vx ]  [ u ]  = [ b_x ]
    # [ A_uy   A_vy ]  [ v ]    [ b_y ]
    #
    #
    # Recall:
    #   n_terms = 9;
    #   n_u     = 2*nx*ny;
    #   n_sprs  = n_u*n_terms;

    # Initialize values to zero
    u .= 0.0;       # [n_u]
    b .= 0.0;       # [n_u]
    Ai .= 0;        # [n_sprs]
    Aj .= 0;        # [n_sprs]
    Av .= 0.0;      # [n_sprs]

    # Equation is being defined for acx-nodes (x-direction equation)

    k = 0 

    for i in 1:nx
        for j in 1:ny 

            # BC: Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            # Set the row in matrix A that the equation is being defined for:
            nr = (i-1)*ny + j

            # -- vx terms -- 
            
            # ux(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(ip1,j,nx,ny);
            Av[k] = ( 4.0/dxdx*N[ip1,j]);

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-4.0/dxdx * (N[ip1,j]+N[i,j])
                     -1.0/dydy * (N_ab[i,j]+N_ab[i,jm1])
                    -β_acx[i,j]);

            # ux(i-1,j)  
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 4.0/dxdx*N[i,j]);

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 1.0/dydy*N_ab[i,j]);

            # ux(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jm1,nx,ny);
            Av[k] = ( 1.0/dydy*N_ab[i,jm1]);
            
            # -- vy terms -- 

            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,j]
                     -1.0/dxdy*N_ab[i,j]);
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = (-2.0/dxdy*N[ip1,j]
                     +1.0/dxdy*N_ab[i,j]);
            
            # uy(i+1,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,jm1,nx,ny);
            Av[k] = (-2.0/dxdy*N[ip1,j]
                     -1.0/dxdy*N_ab[i,jm1]);
            
            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,j]
                     +1.0/dxdy*N_ab[i,jm1]);

            # [u] value
            u[nr] = ux[i,j];

            # [b] value 
            b[nr] = taud_acx[i,j];

        end
    end

    # Equation is being defined for acy-nodes (y-direction equation)

    for i in 1:nx 
        for j in 1:ny

            # BC: Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            # Set the row in matrix A that the equation is being defined for:
            nr = (i-1)*ny + j + nx*ny

            # -- uy terms -- 
            
            # uy(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jp1,nx,ny);
            Av[k] = ( 4.0/dydy*N[i,jp1]);  
            
            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-4.0/dydy * (N[i,jp1]+N[i,j])
                     -1.0/dxdx * (N_ab[i,j]+N_ab[im1,j])
                     -β_acy[i,j]);

            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 4.0/dydy*N[i,j]);    
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = ( 1.0/dxdx*N_ab[i,j]);     
            
            # uy(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(im1,j,nx,ny);
            Av[k] = ( 1.0/dxdx*N_ab[im1,j]);   
            
            # -- ux terms -- 

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,jp1]
                     +1.0/dxdy*N_ab[i,j]);

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,j]
                     -1.0/dxdy*N_ab[i,j]); 

            # ux(i-1,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,jp1,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,jp1]
                     -1.0/dxdy*N_ab[im1,j]);
            
            # ux(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,j]
                     +1.0/dxdy*N_ab[im1,j]);

            # [u] value
            u[nr] = uy[i,j];

            # [b] value 
            b[nr] = taud_acy[i,j];
            
        end
    end
    
    # Now u, b and A components (I, J, V vectors) have been defined.
    # Convert into a sparse array for solving:

    if !ls.prob_is_defined
        # Intitially allocate sparse array with current sparsity pattern
        Asp = sparse(Ai,Aj,Av);

        # Define th linear problem and the cache
        ls.prob = LinearProblem(Asp, b; u0=u)
        ls.cache = init(ls.prob, LUFactorization())  # Or other solver method, like KrylovJL_GMRES()

        # Set to true
        ls.prob_is_defined = true
    end

    # Solve linear problem: x = A/b

# ALLOCATING:
    # prob = LinearProblem(Asp, b; u0=u);
    # sol = solve(prob);                  # Is this allocating??

# NON-ALLOCATING:
    Asp.nzval .= Av
    ls.cache.b .= b
    ls.sol = solve!(ls.cache)  # In-place solve, no allocation

    # Fill output velocity arrays with new solution

    for i = 1:nx
        for j = 1:ny
            n = ij2n_ux(i,j,nx,ny);
            ux[i,j] = ls.sol.u[n];
            n = ij2n_uy(i,j,nx,ny);
            uy[i,j] = ls.sol.u[n];
        end
    end

    return
end