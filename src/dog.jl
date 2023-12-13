export BlobLoG # from ImageFiltering
export BlobDoG

struct BlobDoG{T,S,N} <: AbstractBlob{T,S,N}
    location::CartesianIndex{N}
    σ::S
    amplitude::T
end

function blob_DoG(img::AbstractArray{T,N}, σscales_low, σscales_high;
    edges::Union{Bool,Tuple{Bool,Vararg{Bool,N}}}=(true, ntuple(d -> false, Val(N))...),
    σshape::NTuple{N,Real}=ntuple(d -> 1, Val(N)),
    rthresh::Real=1 // 1000
) where {T<:Union{Real,AbstractGray},N}
    if edges isa Bool
        edges = (edges, ntuple(d -> edges, Val(N))...)
    end
    sigmas_low = sort(σscales_low)
    sigmas_high = sort(σscales_high)
    img_DoG = multiDoG(img, sigmas_low, sigmas_high, σshape)
    maxima = findlocalmaxima(img_DoG; edges)
    if !iszero(rthresh)
        imgmax = maximum(abs, img)
        [
            BlobDoG(CartesianIndex(Base.tail(x.I)), sigmas_low[x[1]] .* σshape, img_DoG[x])
            for x in maxima if img_DoG[x] > rthresh * imgmax
        ]
    else
        [
            BlobDoG(CartesianIndex(Base.tail(x.I)), sigmas_low[x[1]] .* σshape, img_DoG[x])
            for x in maxima
        ]
    end
end

function multiDoG(img::AbstractArray{T,N}, sigmas_low, sigmas_high, σshape;
    l= 2 .* sigmas_high .- 1
) where {T,N}
    (issorted(sigmas_low) && issorted(sigmas_high)) || error("sigmas must be sorted")
    all(sigmas_high .> sigmas_low) || error("the lowpass sigmas must be larger than the highcut sigmas")
    img_DoG = similar(img, float(eltype(T)), (Base.OneTo(length(sigmas_low)), axes(img)...))
    colons = ntuple(d -> Colon(), Val(N))
    @inbounds for (isigma, σs) in enumerate(zip(sigmas_low, sigmas_high))
        σ1, σ2 = σs
        DoG_slice = @view img_DoG[isigma, colons...]
        imfilter!(DoG_slice, img,
            Kernel.DoG(σ1 .* σshape, σ2 .* σshape, l .* σshape),
            "reflect"
        )
        DoG_slice .*= (σ1 / (σ2 - σ1))^2
    end
    return img_DoG
end
