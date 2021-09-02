

function drawlimitset(generators::Vector{F}, z0::Number=0;
  iterations::Int=5) where F <: AbstractComplexInvertibleFunction

  SDDGraphics.newdrawing()

  #SDDGraphics.updatecolorarray(iterations)

  genset = GeneratorsSet(generators)

  SDDGraphics.updatecolorarray(size(genset))

  function procedure(zn::Number, genindex::Int, level::Int, levels::Int)
    if level == levels && SDDGraphics.insiderectregion(zn)
      SDDGraphics.color(genindex)
      SDDGraphics.drawpoint(zn)
    end
  end

  traversegroupfirstdepth(generators, iterations, z0, procedure)

  SDDGraphics.drawing()
end


const drawΛ = drawlimitset


function drawfixedpointslimitset(generators::Vector{F};
  iterations::Int=5) where F <: AbstractMobiusTransformation

  SDDGraphics.newdrawing()

  SDDGraphics.updatecolorarray(iterations)

  #genset = GeneratorSet(generators)

  #SDDGraphics.updatecolorarray(length(genset))

  function procedure(f::AbstractMobiusTransformation, genindex::Int, level::Int, levels::Int)
    #z1, z2 = fixedpoints(f)
    for z in fixedpoints(f) #[z1,z2]
      if SDDGraphics.insiderectregion(z)
        SDDGraphics.color(level)
        SDDGraphics.drawpoint(z)
      end
    end
  end

  traversegroupfirstdepthcomp(generators, iterations, procedure)

  SDDGraphics.drawing()
end


const drawfixedpointsΛ = drawfixedpointslimitset


function drawchaosgamelimitset(generators::Vector{F}, z0::Number=0;
  preiterations::Int=100, iterations::Int=100) where F <: AbstractComplexInvertibleFunction

  SDDGraphics.newdrawing()

  #SDDGraphics.updatecolorarray(iterations)

  genset = GeneratorsSet(generators)

  SDDGraphics.updatecolorarray(size(genset))

  previndex = 0

  function randT(z::Number)
    k = rand(1:size(genset))
    invindex = inverseindex(genset,k)
    while previndex == invindex
      k = rand(1:size(genset))
      invindex = inverseindex(genset,k)
    end
    previndex = k
    genset[k](z)
  end

  zn = z0

  for n in 1:preiterations
    zn = randT(zn)
  end

  for n in 1:iterations
    zn = randT(zn)
    if SDDGraphics.insiderectregion(zn)
      SDDGraphics.color(previndex)
      SDDGraphics.drawpoint(zn)
    end
  end

  SDDGraphics.drawing()
end


const drawchaosgameΛ = drawchaosgamelimitset


function drawtrappedpointslimitset(generators::Vector{F};
  distance::Function=(z1::Number,z2::Number) -> abs2(z1-z2),
  maxiterations::Int=50) where F <: AbstractMobiusTransformation

  # Veryfying if graphics backend supports functions
  SDDGraphics.supported(:drawpixel)

  # Verifying functions
  @assert typeof(distance(1., im)) <: Real

  SDDGraphics.newdrawing()

  SDDGraphics.updatecolorarray(maxiterations)

  finitefixedpoints = ComplexF64[]
  mobiusT = AbstractMobiusTransformation[]
  parabolic = Bool[]

  R = 0

  # ToDo: Handle parabolic transformations with fixed point at Inf.

  for f in generators
    z1, z2 = fixedpoints(f)
    if kind(f) == :parabolic
      if isfinite(z1)
        push!(finitefixedpoints, z1)
        push!(mobiusT, f)
        push!(parabolic, true)
        if R < abs2(z1)
          R = abs2(z1)
        end
      end
    else
      fder = derivative(f)
      for z in [z1,z2]
        if isfinite(z)
          push!(finitefixedpoints, z)
          push!(parabolic, false)
          if abs2(fder(z)) >= 1 # loxodromic or elliptic
            push!(mobiusT, f)
          else # loxodromic
            push!(mobiusT, inverse(f))
          end # if |f'(z)|>=1

          if R < abs2(z)
            R = abs2(z)
          end
        end # if |z|<Inf
      end # for z in [z1,z2]
    end # if else
  end # for generators

  #R *= √2

  # Quasi-PieceWise Transformation, using a Voronoi partition over fixed points
  function voronoiqpwt(z::Number)
    k = 1
    dist = Inf

    for n in 1:length(finitefixedpoints)
      distzfix = distance(z,finitefixedpoints[n])
      if dist > distzfix
        dist = distzfix
        k = n
      end
    end

    if parabolic[k]
      z1 = mobiusT[k](z)
      z2 = inverse(mobiusT[k])(z)
      if distance(z1,finitefixedpoints[k]) >= distance(z2,finitefixedpoints[k])
        return z1
      else
        return z2
      end
    end

    mobiusT[k](z)
  end

  @sweeprectregion SDDGraphics.xlims() SDDGraphics.ylims() SDDGraphics.canvassize() begin
      zn = complex(x,y)
      escapetime = maxiterations
      for n in 1:maxiterations
          if abs2(zn) > R
              escapetime = n
              break
          end # if hasescaped
          zn = voronoiqpwt(zn)
      end # for n maxiterations
      SDDGraphics.color(escapetime)
      SDDGraphics.drawpixel(i,j)
  end # Implemented algorithm

  SDDGraphics.drawing()

end


const drawtrappedpointsΛ = drawtrappedpointslimitset
