Lof = Dict{NTuple{4,Int64},Int64}()
Kof = Dict{NTuple{3,Int64},Int64}()
Jof = Dict{NTuple{2,Int64},Int64}()

KeyofL = Vector{NTuple{4,Int64}}()
KeyofK = Vector{NTuple{3,Int64}}()
KeyofJ = Vector{NTuple{2,Int64}}()

i = 1
for key ∈ Key₃
    global i
    Lof[key] = i
    push!(KeyofL,key)
    i+=1
end

i = 1
for key ∈ Key₂
    global i
    Kof[key] = i
    push!(KeyofK,key)
    i+=1
end

i = 1
for key ∈ Key₁
    global i
    Jof[key] = i
    push!(KeyofJ,key)
    i+=1
end

#save("numbering.jld2","Lof",L,"Kof",K,"Jof",J,"KeyofL",KeyofL,"KeyofK",KeyofK,"KeyofJ",KeyofJ)