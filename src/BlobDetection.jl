module BlobDetection

using ImageFiltering
using PaddedViews
using ColorTypes

export AbstractBlob

abstract type AbstractBlob{T,S,N} end

function Base.show(io::IO, blob::AbstractBlob)
    l = location(blob)
    σ = scale(blob)
    a = amplitude(blob)
    m0 = zeroth_moment(blob)
    m2 = second_moment(blob)
    print(io, 
        "Blob(location=$l, σ=$σ, amplitude=$a, m₀=$m0, m₂=$m2)"
    )
end

include("image_moments.jl")
include("blob_api.jl")
include("blobs.jl")
include("location_refinement.jl")
#include("filtering.jl")

end
