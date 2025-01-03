cd("/home/leem/Dropbox/workspace/scipy/3dsolve/")

using JLD2

#points.jld2 contains 35 points' xyz-coordinates
points = load("./points.jld2")["points"]

#now we collect 140 3-simplices
function addtetra(pts)
    push!(Key₃,Tuple(sort([pts[1],pts[2],pts[3],pts[6]])))
    push!(Key₃,Tuple(sort([pts[4],pts[5],pts[6],pts[2]])))
    push!(Key₃,Tuple(sort([pts[1],pts[4],pts[6],pts[2]])))
end

# 꼭지점 부터 1층 별 10개 simplices
Key₃ = [Tuple(sort([1,2,4,6]))]
push!(Key₃,Tuple(sort([1,2,6,8])))
push!(Key₃,Tuple(sort([1,2,8,10])))
push!(Key₃,Tuple(sort([1,2,10,12])))
push!(Key₃,Tuple(sort([1,2,12,4])))
push!(Key₃,Tuple(sort([1,4,5,6])))
push!(Key₃,Tuple(sort([1,6,7,8])))
push!(Key₃,Tuple(sort([1,8,9,10])))
push!(Key₃,Tuple(sort([1,10,11,12])))
push!(Key₃,Tuple(sort([1,12,3,4])))

#별1층-별2층 기둥 바깥쪽
AtoF = [(4,5,6,15,16,17)]
push!(AtoF, (6,7,8,17,18,19))
push!(AtoF, (8,9,10,19,20,21))
push!(AtoF, (10,11,12,21,22,23))
push!(AtoF, (12,3,4,23,14,15))
for pts ∈ (AtoF)
    addtetra(pts)
end

#별1층-별2층 기둥 안쪽
AtoF = [(4,2,6,15,13,17)]
push!(AtoF, (6,2,8,17,13,19))
push!(AtoF, (8,2,10,19,13,21))
push!(AtoF, (10,2,12,21,13,23))
push!(AtoF, (12,2,4,23,13,15))
for pts ∈ (AtoF)
    addtetra(pts)
end

#별2층-별3층 기둥 바깥쪽
AtoF = [(15,16,17,26,27,28)]
push!(AtoF, (17,18,19,28,29,30))
push!(AtoF, (19,20,21,30,31,32))
push!(AtoF, (21,22,23,32,33,34))
push!(AtoF, (23,14,15,34,25,26))
for pts ∈ (AtoF)
    addtetra(pts)
end

#별2층-별3층 기둥 안쪽
AtoF = [(15,13,17,26,24,28)]
push!(AtoF, (17,13,19,28,24,30))
push!(AtoF, (19,13,21,30,24,32))
push!(AtoF, (21,13,23,32,24,34))
push!(AtoF, (23,13,15,34,24,26))
for pts ∈ (AtoF)
    addtetra(pts)
end

push!(Key₃,Tuple(sort([24,26,28,35])))
push!(Key₃,Tuple(sort([24,28,30,35])))
push!(Key₃,Tuple(sort([24,30,32,35])))
push!(Key₃,Tuple(sort([24,32,34,35])))
push!(Key₃,Tuple(sort([24,34,26,35])))

push!(Key₃,Tuple(sort([26,27,28,35])))
push!(Key₃,Tuple(sort([28,29,30,35])))
push!(Key₃,Tuple(sort([30,31,32,35])))
push!(Key₃,Tuple(sort([32,33,34,35])))
push!(Key₃,Tuple(sort([34,25,26,35])))

#save("./C3.jld2","C3",Key₃)

#now, we collect 2-simplices
Key₂set = Key₂set = Set{Set{Int64}}()
for tetra ∈ Key₃
    push!(Key₂set,Set(  [  tetra[2],tetra[3],tetra[4]  ]  )  )
    push!(Key₂set,Set(  [  tetra[1],tetra[3],tetra[4]  ]  )  )
    push!(Key₂set,Set(  [  tetra[1],tetra[2],tetra[4]  ]  )  )
    push!(Key₂set,Set(  [  tetra[1],tetra[2],tetra[3]  ]  )  )
end
Key₂ = Vector{NTuple{3, Int64}}()
for triangle ∈ Key₂set
    push!(Key₂, Tuple(sort(collect(triangle))) )
end
#save("./C2.jld2","C2",Key₂)

#now, we collect 1-simplices
Key₁set = Set{Set{Int64}}()
for triangle ∈ Key₂
    push!(Key₁set,Set(  [  triangle[2],triangle[3]  ]  )  )
    push!(Key₁set,Set(  [  triangle[1],triangle[3]  ]  )  )
    push!(Key₁set,Set(  [  triangle[1],triangle[2]  ]  )  )
end
Key₁ = Vector{NTuple{2, Int64}}()
for line ∈ Key₁set
    push!(Key₁, Tuple(sort(collect(line))) )
end
#save("./C1.jld2","C1",Key₁)


Key₀ = collect(1:35)
#save("./C0.jld2","C0",Key₀)

Key₀ = Tuple(Key₀)
Key₁ = Tuple(sort(Key₁))
Key₂ = Tuple(sort(Key₂))
Key₃ = Tuple(Key₃) 

#save("./complex_keys.jld2","Key0",Key₀,"Key1",Key₁,"Key2",Key₂,"Key3",Key₃)

