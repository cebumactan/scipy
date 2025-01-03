using LinearAlgebra

cd("/home/leem/Dropbox/workspace/scipy/3dsolve/")
#include("./generate_complex.jl")

# numbering = load("./numbering.jld2")
# Lof = numbering["Lof"]
# Kof = numbering["Kof"]
# Jof = numbering["Jof"]
# KeyofL = numbering["KeyofL"]
# KeyofK = numbering["KeyofK"]
# KeyofJ = numbering["KeyofJ"]


D₀ = zeros(n₁,n₀)
for key ∈ Key₁
    mass1 = C₁[key].mass
    j = Jof[key]
    (i₁, i₂) = C₁[key].boundary
    D₀[j,i₂] =  1.0/mass1
    D₀[j,i₁] = -1.0/mass1
end

println("The rank of D₀ is $(rank(D₀)) and this must be n₀-1=$(n₀-1).")

D₁ = zeros(n₂,n₁)
for key ∈ Key₂
    mass2 = C₂[key].mass
    k = Kof[key]
    signfactor = 1
    for bkey ∈ C₂[key].boundary
        mass1 = C₁[bkey].mass
        j = Jof[bkey]
        D₁[k,j] =  signfactor*mass1/mass2
        signfactor *= -1
    end
end

println("The rank of D₁ is $(rank(D₁)) and this must be n₁-n₀+1=$(n₁-n₀+1)=$(n₂-n₃)=n₂-n₃.")

D₂ = zeros(n₃,n₂)
for key ∈ Key₃
    mass3 = C₃[key].mass
    l = Lof[key]
    signfactor = 1
    for bkey ∈ C₃[key].boundary
        mass2 = C₂[bkey].mass
        k = Kof[bkey]
        D₂[l,k] =  signfactor*mass2/mass3
        signfactor *= -1
    end
end

println("The rank of D₂ is $(rank(D₂)) and this must be n₃=$n₃.")
