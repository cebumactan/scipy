cd("/home/leem/Dropbox/workspace/scipy/3dsolve/")

using JLD2, LinearAlgebra

points = load("./points.jld2")["points"]
Key = load("./complex_keys.jld2")
Key₀ = Key["Key0"]
Key₁ = Key["Key1"]
Key₂ = Key["Key2"]
Key₃ = Key["Key3"]

n₀ = length(Key₀)
n₁ = length(Key₁)
n₂ = length(Key₂)
n₃ = length(Key₃)

println("컴플렉스에는 $n₀ 개의 points, $n₁ 개의 lines, $n₂ 개의 triangles, $n₃ 개의 tetrahedrons가 담기게 됩니다.")
println("Euler Characteristic n₀-n₁+n₂-n₃=$(n₀-n₁+n₂-n₃)입니다. 값이 1이 맞는지 확인하십시오.")

mutable struct ThreeSimplex
    key::NTuple{4,Int64}
    mass::Float64
    boundary::NTuple{4,NTuple{3,Int64}}
end

mutable struct TwoSimplex
    key::NTuple{3,Int64}
    mass::Float64
    boardof::Vector{NTuple{4,Int64}}
    faceof::Vector{NTuple{4,Int64}}
    multiplicity::Int64 
    boundary::NTuple{3,NTuple{2,Int64}}
end

mutable struct OneSimplex
    key::Tuple{Int64,Int64}
    mass::Float64
    boardof::Vector{NTuple{3,Int64}}
    faceof::Vector{NTuple{4,Int64}}
    multiplicity::Int64 
    boundary::NTuple{2,Int64}
end

mutable struct ZeroSimplex
    key::Int64
    boardof::Vector{NTuple{2,Int64}}
    faceof::Vector{NTuple{4,Int64}}
    multiplicity::Int64 
end

#C₃ constructJtoKey
C₃ = Dict{NTuple{4,Int64},ThreeSimplex}()
for t ∈ Key₃
    A = zeros(3,3)
    A[:,1]= points[t[2]] - points[t[1]]
    A[:,2]= points[t[3]] - points[t[1]]
    A[:,3]= points[t[4]] - points[t[1]]
    mass3 = abs(det(A)) / 6.0
    boundary = ( (t[2],t[3],t[4]), (t[1],t[3],t[4]), (t[1],t[2],t[4]), (t[1],t[2],t[3]) )
    C₃[t] = ThreeSimplex(t,mass3,boundary)
end

#C₂ construct
C₂ = Dict{NTuple{3,Int64},TwoSimplex}()
for t ∈ Key₂
    A = zeros(3,2)
    A[:,1]= points[t[2]] - points[t[1]]
    A[:,2]= points[t[3]] - points[t[1]]
    mass2 = sqrt(det( transpose(A)*A )) / 2.0
    boundary = ( (t[2],t[3]), (t[1],t[3]), (t[1],t[2]) )
    C₂[t] = TwoSimplex(t,mass2,[],[],0,boundary)
end

#C₁ construct
C₁ = Dict{NTuple{2,Int64},OneSimplex}()
for t ∈ Key₁
    A = points[t[2]] - points[t[1]]
    mass1 = sqrt(transpose(A)*A)
    C₁[t] = OneSimplex(t,mass1,[],[],0,(t[2],t[1]))
end

#C₀ construct
C₀ = Vector{ZeroSimplex}()
for i ∈ 1:n₀
    push!(C₀,ZeroSimplex(i,[],[],0))
end

for t ∈ Key₃
    for j ∈ (1,2,3,4)
        bkey = C₃[t].boundary[j]
        push!(C₂[bkey].boardof,t)
        push!(C₂[bkey].faceof,t)
        C₂[bkey].multiplicity+=1
    end

    for fkey ∈ ( (t[1],t[2]) , (t[1],t[3]) , (t[1],t[4]) , (t[2],t[3]), (t[2],t[4]) , (t[3],t[4]) )
        push!(C₁[fkey].faceof,t)
        C₁[fkey].multiplicity+=1
    end

    for fkey ∈ t
        push!(C₀[fkey].faceof,t)
        C₀[fkey].multiplicity+=1
    end
end

for t ∈ Key₂
    for j ∈ (1,2,3)
        b = C₂[t].boundary[j]
        push!(C₁[b].boardof,t)
    end
end

for t ∈ Key₁
    for j ∈ (1,2)
        b = C₁[t].boundary[j]
        push!(C₀[b].boardof,t)
    end
end

#save("./complex.jld2","C0",C₀,"C1",C₁,"C2",C₂,"C3",C₃)

total = 0
for t ∈ values(C₂)
    global total
    total += t.multiplicity
end
println(total)
println(4*n₃)

total = 0
for t ∈ values(C₁)
    global total
    total += t.multiplicity
end
println(total)
println(6*n₃)

total = 0
for t ∈ values(C₀)
    global total
    total += t.multiplicity
end
println(total)
println(4*n₃)






