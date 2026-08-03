
foo(x) = (1 / x + 1) / 2
M = rand(10000, 10000)
mask = rand(size(M)...) .< 0.1
function foo1(Y, X::AbstractMatrix)
    Y .= (1 ./ X .+ 1) ./ 2
    return nothing
end
function foo2(Y, X::AbstractMatrix)
    @inbounds for I in CartesianIndices(X)
        Y[I] = foo(X[I])
    end
    return nothing
end
idx = CartesianIndices(mask)
mask_idx = idx[mask]
@b map!(foo, M, M)
@b map!(foo, M)
# @b map!(foo, view(M, mask))
@b map!(foo, view(M, mask_idx))
@b foo1(M, M)
@b foo2(M, M)
@b map!(foo, view(M, mask), view(M, mask))
@b view(M, mask)



using KernelAbstractions
@kernel function foo3(Y, @Const(X))
    I = @index(Global)
    @inbounds Y[I] = foo(X[I])
end

function myfoo3(A, B)
    backend = get_backend(A)
    kernel = foo3(backend)
    kernel(A, B, ndrange = length(A))
    KernelAbstractions.synchronize(backend)
end
@b myfoo3(M2, M)


@kernel function foo4(Y, @Const(X))
    I = @index(Global)
    @inbounds Y[I] = foo(X[I])
end

function myfoo4(Y, X)
    backend = get_backend(Y)
    kernel = foo4(backend, 32, size(Y)) # if size(A) varies this will cause recompilation
    kernel(Y, X, ndrange = size(Y))
    KernelAbstractions.synchronize(backend)
    return
end

@b myfoo4(M2, M)


# Check when using KernelAbstractions, memory is allocated, even with static size.
# Maybe this just can't be avoided...
@kernel function copy_kernel!(A, @Const(B))
    I = @index(Global)
    @inbounds A[I] = B[I]
end

function mycopy_static!(A, B)
    backend = get_backend(A)
    @assert size(A) == size(B)
    @assert get_backend(B) == backend

    kernel = copy_kernel!(backend, 32, size(A)) # if size(A) varies this will cause recompilation
    kernel(A, B, ndrange = size(A))
    return
end

A = rand(128, 128)
B = rand(128, 128)
mycopy_static!(A, B)
@b mycopy_static!(A, B)