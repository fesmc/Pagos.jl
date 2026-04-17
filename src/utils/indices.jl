abstract type AbstractIndexing end

struct StrictIndexing <: AbstractIndexing
    i1::Int
    i2::Int
end
struct FlatIndexing <: AbstractIndexing
    i1::Int
    i2::Int
end
struct ReflectiveIndexing <: AbstractIndexing
    i1::Int
    i2::Int
end
struct PeriodicIndexing <: AbstractIndexing
    i1::Int
    i2::Int
end

function index(i, d, idx::StrictIndexing)
    j = i + d
    if j < idx.i1 || j > idx.i2
        error("Index $j out of bounds [$(idx.i1), $(idx.i2)]")
    end
    return j
end
function index(i, d, idx::FlatIndexing)
    j = i + d
    if j < idx.i1
        j = idx.i1
    elseif j > idx.i2
        j = idx.i2
    end
    return j
end
function index(i, d, idx::ReflectiveIndexing)
    j = i + d
    if j < idx.i1
        j = idx.i1 + (idx.i1 - j)
    elseif j > idx.i2
        j = idx.i2 - (j - idx.i2)
    end
    return j
end
function index(i, d, idx::PeriodicIndexing)
    j = i + d
    if j < idx.i1
        j = idx.i2 - (idx.i1 - j - 1)
    elseif j > idx.i2
        j = idx.i1 + (j - idx.i2 - 1)
    end
    return j
end

function stencil(i, idx)
    return index(i, -1, idx), index(i, 1, idx)
end
function stencil(i, j, i_idx, j_idx)
    return stencil(i, i_idx)..., stencil(j, j_idx)...
end

function ij2n_ux(i,j,nx,ny)

    n = (i-1)*ny + j

    return n 
end

function ij2n_uy(i,j,nx,ny)

    n = (i-1)*ny + j + nx*ny
    
    return n 
end


function von_neumann_neighbours(I::CartesianIndex{N}) where N   #, bc::FlatBC
    neighbours = CartesianIndex[]
    for d in 1:N
        push!(neighbours, I + CartesianIndex(ntuple(i -> i == d ? -1 : 0, N)))
        push!(neighbours, I + CartesianIndex(ntuple(i -> i == d ? 1 : 0, N)))
    end
    return neighbours
end

function moore_neighbours(I::CartesianIndex{N}) where N   #, bc::FlatBC
    neighbours = CartesianIndex[]
    for offset in Iterators.product(ntuple(_ -> -1:1, N)...)
        if any(x -> x != 0, offset)
            push!(neighbours, I + CartesianIndex(offset))
        end
    end
    return neighbours
end