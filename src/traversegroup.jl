
#import Base: length, getindex


struct GeneratorsSet{F <: AbstractComplexInvertibleFunction}
  functions::Vector{F}
  numregular::Int64

  function GeneratorsSet{F}(gs::Vector{F}) where F <: AbstractComplexInvertibleFunction
    funs = (F)[] # Vector{F}
    autoinverses = (F)[] # Vector{F}

    for f in gs
      if f isa AbstractMobiusTransformation && !(f isa InversionReflection)
        push!(funs, f)
      else
        push!(autoinverses, f)
      end # if
    end # for

    funsinv = [inverse(f) for f in funs]
    funs = append!(funs, funsinv)
    nfirstsfuns = length(funs)
    append!(funs, autoinverses)

    new(funs, nfirstsfuns)
  end # constructor GeneratorsSet
end # struct GeneratorsSet


GeneratorsSet(gs::Vector{F}) where F <: AbstractComplexInvertibleFunction =
  GeneratorsSet{F}(gs)


function inverseindex(gs::GeneratorsSet, n::Int)
  if n > gs.numregular
    return n
  end
  ((n + (gs.numregular ÷ 2) - 1) % gs.numregular) + 1
end

Base.getindex(gs::GeneratorsSet, i::Int) = gs.functions[i]
Base.firstindex(gs::GeneratorsSet) = gs.functions[begin]
Base.lastindex(gs::GeneratorsSet) = gs.functions[end]

size(gs::GeneratorsSet) = length(gs.functions)
numregulartransformations(gs::GeneratorsSet) = gs.numregular
numregulartransformationsnoinverses(gs::GeneratorsSet) = gs.numregular÷2
numautoinverses(gs::GeneratorsSet) = length(gs.functions) - gs.numregular
numgenerators(gs::GeneratorsSet) = (gs.numregular÷2) + numautoinverses(gs)


function applypaired(genset::GeneratorsSet, level::Int, levels::Int, objects::Vector, procedure::Function)
  newobjects = []
  newtags = Int[]

  for n in 1:size(genset)
    newset = []
    invindex = inverseindex(genset, n)

    for oindex in 1:length(objects)
      if oindex != invindex # if object not paired with inverse
        o = genset[n](objects[oindex])
        push!(newset, o)
      end
    end # for levelobjects

    procedure(newset, n, level, levels) # procedure on objects Vector
    push!(newobjects, newset)
    push!(newtags, n)
  end # for genset

  newobjects, newtags
end


function traversegroupfirstdepth(generators::Vector{F}, levels::Int, object,
  procedure::Function, paired::Bool=false) where F <: AbstractComplexInvertibleFunction

  #println("- Level 1, gen 0")
  procedure(object, 0, 1, levels) # Level 1: Identity transformation ("generator index = 0")

  if levels <= 1
    return
  end

  genset = GeneratorsSet(generators)

  if paired # object is a vector of elements paired with the transformations
    objects, tags = applypaired(genset, 2, levels, object, procedure)
    if levels <= 2
      return
    end
    for n in 1:size(genset) # Level 3
      invindex = inverseindex(genset,n)
      for k in 1:length(objects)
        if k != invindex
          traversedepth(genset, n, 3, levels, genset[n](objects[k]), procedure)
        end
      end
    end
  else # else paired
    for n in 1:size(genset) # Level 2: Apply transformations first time
      traversedepth(genset, n, 2, levels, genset[n](object), procedure)
    end
  end # if else paired
end


function traversedepth(genset::GeneratorsSet, genindex::Int, level::Int,
  levels::Int, object, procedure::Function)

  #println("-"^level * " Level $level, gen $genindex")
  procedure(object, genindex, level, levels)

  if level == levels
    return
  end

  N = size(genset)
  invindex = inverseindex(genset, genindex)

  # Traverse group (Cayley graph) to lower nodes,
  # but avoiding the corresponding inverse node

  for n in (invindex+1):N
    traversedepth(genset, n, level+1, levels, genset[n](object), procedure)
  end

  for n in 1:(invindex-1)
    traversedepth(genset, n, level+1, levels, genset[n](object), procedure)
  end
end


function traversegroupfirstdepthcomp(generators::Vector{F}, levels::Int,
  procedure::Function) where F <: AbstractMobiusTransformation

  #for g in generators
  #  procedure(g, 0, 1, levels) # Level 1: Identity transformation ("generator index = 0")
  #end

  if levels <= 0
    return
  end

  genset = GeneratorsSet(generators)

  for n in 1:length(generators)
    traversedepthcomp(genset, n, 1, levels, genset[n], procedure)
  end
end


function traversedepthcomp(genset::GeneratorsSet{F}, genindex::Int, level::Int,
  levels::Int, T::F, procedure::Function) where F <: AbstractMobiusTransformation

  procedure(T, genindex, level, levels)

  if level == levels
    return
  end

  N = size(genset)
  invindex = inverseindex(genset, genindex)

  # Traverse group (Cayley graph) to lower nodes,
  # but avoiding the corresponding inverse node and auto iteratives g^n

  for n in (invindex+1):N
    #if n != genindex
      traversedepthcomp(genset, n, level+1, levels, genset[n]∘T, procedure)
    #end
  end

  for n in 1:(invindex-1)
    #if n != genindex
      traversedepthcomp(genset, n, level+1, levels, genset[n]∘T, procedure)
    #end
  end
end



function traversesemigroup(generators::Vector{F}, levels::Int, object,
  procedure::Function) where F <: AbstractMobiusTransformation

  procedure(object, 0, 1, levels)

  for n in 1:length(generators)
    traversesemigroup(generators, n, 1, levels, generators[n](object), procedure)
  end
end


function traversesemigroup(generators::Vector{F}, genindex::Int, level::Int,
  levels::Int, object, procedure::Function) where F <: AbstractMobiusTransformation

  procedure(object, genindex, level, levels)

  if level == levels
    return
  end

  for n in length(generators)
    traversesemigroup(generators, n, level+1, levels, generators[n](object), procedure)
  end
end



function traversegroupfirstbreath(generators::Vector{F}, levels::Int, object,
  procedure::Function, paired::Bool=false) where F <: AbstractComplexInvertibleFunction

  procedure(object, 0, 1, levels) # Level 1: Identity transformation ("generator index = 0")

  if levels <= 1
    return
  end

  genset = GeneratorsSet(generators)
  N = size(genset)

  currentlevel = 2

  newobjects = [object]
  newtags = [0]

  if paired
    newobjects, newtags = applypaired(genset, currentlevel, levels, object, procedure)
    #newobjects = [object]
    currentlevel = 3
  end # if paired

  if levels < currentlevel
    return
  end

  for l in currentlevel:levels
    levelobjects = newobjects
    leveltags = newtags
    newobjects = []
    newtags = Int[]

    for oindex in 1:length(levelobjects)
      for n in 1:N
        invindex = inverseindex(genset, n)

        if leveltags[oindex] == invindex
          continue
        end

        o = genset[n](levelobjects[oindex])
        procedure(o, n, l, levels)
        push!(newobjects, o)
        push!(newtags, n)
      end # for genset
    end # for levelobjects

    empty!(levelobjects)
    empty!(leveltags)
  end # for levels

  empty!(newobjects)
  empty!(newtags)

end # function traversebreath
