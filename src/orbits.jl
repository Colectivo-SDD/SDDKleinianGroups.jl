
###
### Groups with 1 generators, cyclic groups
###

"""
    drawpointorbit(f, z0 [; preiterations, iterations])

Return the drawing of the orbit of \$z_0\$ under a complex invertible function.

The orbit of \$D\$ under \$f\$ is defined as
  \$\\big<f\\big>z_0=\\{\\dots,f^{-n}(z_0),\\dots,f^{-1}(z_0),z_0,f^(z_0),\\dots,f^n(z_0),\\dots\\}\$

#### Arguments
- `f::AbstractComplexInvertibleFunction`: A Möbius transformation or inversion.
- `z0::Number`: A complex number.
- `preiterations::Integer`: Number of first iterations to calculate but not to draw.
- `iterations::Integer`: Number of iterations to calculate (after `preiterations`).
"""
function drawpointorbit(f::AbstractComplexInvertibleFunction, z0::Number;
  preiterations::Int=0, iterations::Int=20)

  finv = inverse(f)

  SDDGraphics.newdrawing()
  SDDGraphics.updatecolorarray(iterations)

  zn = z0
  zninv = z0

  for n in 1:preiterations
      zn = f(zn)
      zninv = finv(zninv)
  end # for n preiterations

  for n in 1:iterations
      SDDGraphics.color(n)

      if SDDGraphics.insiderectregion(zn)
          SDDGraphics.drawpoint(zn)
      end # if
      zn = f(zn)

      if SDDGraphics.insiderectregion(zninv)
          SDDGraphics.drawpoint(zninv)
      end # if
      zninv = finv(zninv)
  end # for n iterations

  SDDGraphics.drawing()
end # function drawpointorbit


function drawpointssetorbit(f::AbstractComplexInvertibleFunction,
  ps::Array{N,1}; preiterations::Int=0, iterations::Int=20) where N <: Number

  finv = inverse(f)

  SDDGraphics.newdrawing()
  SDDGraphics.updatecolorarray(iterations)

  pns = deepcopy(ps)
  pnsinv = deepcopy(ps)

  for n in 1:preiterations
      pns = f.(pns)
      pnsinv = finv.(pnsinv)
  end # for n preiterations

  npns = length(pns)

  for n in 1:iterations
      SDDGraphics.color(n)

      for k in 1:npns
          if SDDGraphics.insiderectregion(pns[k])
              SDDGraphics.drawpoint(pns[k])
          end # if
          pns[k] = f(pns[k])

          if SDDGraphics.insiderectregion(pnsinv[k])
              SDDGraphics.drawpoint(pnsinv[k])
          end # if
          pnsinv[k] = finv(pnsinv[k])
      end # for k pns
  end # for n iterations

  SDDGraphics.drawing()
end # function drawpointssetorbit


function drawcircleorbit(f::AbstractComplexInvertibleFunction,
  c0::AbstractCircularLinearCurve; preiterations::Int=0, iterations::Int=20)

  finv = inverse(f)

  SDDGraphics.newdrawing()
  SDDGraphics.updatecolorarray(iterations)

  cn = c0
  cninv = c0

  for n in 1:preiterations
      cn = f(cn)
      cninv = finv(cninv)
  end # for n preiterations

  for n in 1:iterations
      SDDGraphics.color(n)

      if cn isa Circle
          SDDGraphics.drawcircle(cn)
      else
          SDDGraphics.drawline(cn)
      end # if
      cn = f(cn)

      if cninv isa Circle
          SDDGraphics.drawcircle(cninv)
      else
          SDDGraphics.drawline(cninv)
      end # if
      cninv = finv(cninv)
  end # for n iterations

  SDDGraphics.drawing()
end # function drawcircleorbit


function drawcirclessetorbit(f::AbstractComplexInvertibleFunction, cs::Vector{C};
  preiterations::Int=0, iterations::Int=20) where C <: AbstractCircularLinearCurve

  finv = inverse(f)

  SDDGraphics.newdrawing()
  SDDGraphics.updatecolorarray(iterations)

  cns = deepcopy(cs)
  cnsinv = deepcopy(cs)

  for n in 1:preiterations
    cns = f.(cns)
    cnsinv = finv.(cnsinv)
  end # for n preiterations

  ncns = length(cns)

  for n in 1:iterations
    SDDGraphics.color(n)

    for k in 1:ncns
      if cns[k] isa Circle
        SDDGraphics.drawcircle(cns[k])
      else
        SDDGraphics.drawline(cnk[k])
      end # if
      cns[k] = f(cns[k])

      if cnsinv[k] isa Circle
        SDDGraphics.drawcircle(cnsinv[k])
      else
        SDDGraphics.drawline(cnsinv[k])
      end # if
      cnsinv[k] = finv(cnsinv[k])
    end
  end # for n iterations

  SDDGraphics.drawing()
end # function drawcircleorbit


#=
"""
    drawdiscorbit(f, disc [; preiterations, iterations])

Return the drawing of the orbit of a disc \$D\$ under a complex invertible function.

The orbit of \$D\$ under \$f\$ is defined as
  \$\\big<f\\big>D=\\{\\dots,f^{-n}(D),\\dots,f^{-1}(D),D,f^(D),\\dots,f^n(D),\\dots\\}\$

#### Arguments
- `f::AbstractComplexInvertibleFunction`: A Möbius transformation or inversion.
- `disc::AbstractDisc`: A disc in the Riemann sphere
- `preiterations::Integer`: Number of first iterations to calculate but not to draw.
- `iterations::Integer`: Number of iterations to calculate (after `preiterations`).
"""
function drawdiscorbit(f::AbstractComplexInvertibleFunction, d0::AbstractDisc;
  preiterations::Int=0, iterations::Int=20)

  finv = inverse(f)

  SDDGraphics.newdrawing()
  SDDGraphics.updatecolorarray(iterations)

  dn = d0
  dninv = d0

  for n in 1:preiterations
      dn = f(dn)
      dninv = finv(dninv)
  end # for n preiterations

  for n in 1:iterations
      SDDGraphics.color(n)

      if iscircle(cn)
          SDDGraphics.drawcircle(cn)
      else
          SDDGraphics.drawline(cn)
      end # if
      cn = f(cn)

      if iscircle(cninv)
          SDDGraphics.drawcircle(cninv)
      else
          SDDGraphics.drawline(cninv)
      end # if
      cninv = finv(cninv)
  end # for n iterations

  SDDGraphics.drawing()
end
=#


###
### Groups with n generators
###

function drawpointorbit(generators::Vector{F}, z0::Number=0;
  preiterations::Int=0, iterations::Int=5) where F <: AbstractComplexInvertibleFunction

  SDDGraphics.newdrawing()

  levels = preiterations + iterations + 1
  SDDGraphics.updatecolorarray(levels)

  function procedure(zn::Number, genindex::Int, level::Int, levels::Int)
    if SDDGraphics.insiderectregion(zn) && level > preiterations
      SDDGraphics.color(level)
      SDDGraphics.drawpoint(zn)
    end
  end

  traversegroupfirstdepth(generators, levels, z0, procedure)

  SDDGraphics.drawing()
end


function drawpointssetorbit(generators::Vector{F}, ps::Vector{N};
  preiterations::Int=0, iterations::Int=5)  where {F <: AbstractMobiusTransformation, N <: Number}

  SDDGraphics.newdrawing()

  levels = preiterations + iterations + 1
  SDDGraphics.updatecolorarray(levels)

  function procedure(psn::Vector{N}, genindex::Int, level::Int, levels::Int) where N <: Number
    if level > preiterations
      SDDGraphics.color(level)
      for z in psn
        if SDDGraphics.insiderectregion(z)
          SDDGraphics.drawpoint(z)
        end
      end
    end
  end

  traversegroupfirstdepth(generators, levels, ps, procedure)

  SDDGraphics.drawing()
end

function drawcircleorbit(generators::Vector{F}, c0::AbstractCircularLinearCurve;
  preiterations::Int=0, iterations::Int=5, ε::Float64=0.0001) where F <: AbstractComplexInvertibleFunction

  SDDGraphics.newdrawing()

  levels = preiterations + iterations + 1
  SDDGraphics.updatecolorarray(levels)

  function procedure(cn, genindex::Int, level::Int, levels::Int)
    if cn isa Circle && radius(cn) < ε
      return
    end
    if level > preiterations
      SDDGraphics.color(level)
      if cn isa Circle
        SDDGraphics.drawcircle(cn)
      else cn isa Line
        SDDGraphics.drawline(cn)
      end
    end
  end

  traversegroupfirstdepth(generators, levels, c0, procedure)

  SDDGraphics.drawing()
end


function drawcirclessetorbit(generators::Vector{F}, cs::Vector{C};
  preiterations::Int=0, iterations::Int=5, ε::Float64=0.0001, color=:levels,
  traverse=:firstdepth, paired::Bool=false) where {F <: AbstractComplexInvertibleFunction,
  C <: AbstractCircularLinearCurve}

  SDDGraphics.newdrawing()

  levels = preiterations + iterations + 1

  if color == :generators
    SDDGraphics.updatecolorarray(2length(generators))
  elseif color == :circles
    SDDGraphics.updatecolorarray(length(cs))
  else
    SDDGraphics.updatecolorarray(levels)
  end

  function procedure(cns, genindex::Int, level::Int, levels::Int)
    if level > preiterations
      if color == :generators
        if genindex == 0
          SDDGraphics.color(SDDGraphics.fgcolor())
        else
          SDDGraphics.color(genindex)
        end
      elseif color == :levels
        SDDGraphics.color(level)
      end
      for cn in cns
        if cn isa Circle
          if radius(cn) < ε
            continue
          end
          SDDGraphics.drawcircle(cn)
        else cn isa Line
          SDDGraphics.drawline(cn)
        end
      end
    end
  end

  if traverse == :firstbreath
    traversegroupfirstbreath(generators, levels, cs, procedure, paired)
  elseif traverse == :semigroup
    traversesemigroup(generators, levels, cs, procedure)
  else
    traversegroupfirstdepth(generators, levels, cs, procedure, paired)
  end

  SDDGraphics.drawing()
end
