module MagGrad

function dvdr(r::Real, θ::Real) :: Real
	return r^(-1.5)*cos(θ)^2/√(4-3*sin(θ)^2)
end 

function dvdt(r::Real, θ::Real) :: Real
	return -2*r^(-0.5)*sin(θ)*cos(θ)/√(4-3*sin(θ)^2)*(3*cos(θ)^2/(4-3*sin(θ)^2)-2)
end

function dudt(r::Real, θ::Real) :: Real
	return -r^(-1.5)/√(4-3*sin(θ)^2)*(cos(θ)^2-sin(θ)^2+3*sin(θ)^2*cos(θ)^2/(4-3*sin(θ)^2))
end

function dudr(r::Real, θ::Real) :: Real
	return 3*r^(-2.5)*cos(θ)*sin(θ)/2/√(4-3*sin(θ)^2)
end

function u(r::Real, θ::Real) :: Real
	return -√(1/r)*cos(θ)^2*2/√(4-3*sin(θ)^2)
end

function v(r::Real, θ::Real) :: Real
	return -√(1/r^3)*cos(θ)*sin(θ)/√(4-3*sin(θ)^2)
end
 
function grad(r::Real, θ::Real, α::Real, β::Real) :: Real
	n = [cos(α), 1/r*sin(α)*cos(β), 1/r/sin(θ)*sin(α)*sin(β)]
	dv = [dvdr(r,θ) dvdt(r,θ); dudr(r,θ) dudt(r,θ)]
	return (dv[1,1]*n[1]^2 + r*sin(θ)^2*n[3]^2*v(r,θ) + r^2*dv[2,1]*n[2]*n[1]+
                    dv[1,2]*n[1]*n[2] + r*v(r,θ)*n[2]^2 + r^2*dv[2,2]*n[2]^2+
                    r^2*sin(θ)*cos(θ)*u(r,θ)*n[3]^2)
end

function grad_integrated_all_directions(r :: Real, θ :: Real; n = 100 :: Integer) :: Real
	αarray = [0:π/n:π;]
	βarray = [0:2*π/n:2*π;]
	s = 0
	for α ∈ αarray, β ∈ βarray
		s = s + abs(grad(r, θ, α, β)) * sin(α)
	end
	return s
end

function grad_integrated_over_star(r :: Real, θ :: Real; n = 100 :: Integer) :: Real
	ϕ = 0

	try
		ϕ = asin(1/r)
	catch e
		if isa(e, DomainError)
			ϕ = π/2
		end
	end

	αarray = [π-ϕ:ϕ/n:π;]
	βarray = [0:2*π/n:2*π;]
	s = 0
	for α ∈ αarray, β ∈ βarray
		s = s + abs(grad.(r, θ, α, β)) * sin(α)
	end
	return s
end

function get_source(rmi :: Real, rmo :: Real; grid = (5,20)) :: Array{Float64, 2}
	mr, mθ = grid

	rmags = [rmi:(rmo-rmi)/(mr-1):rmo;]
	source = zeros(Float64, 5, 20)
  	i = 1
	for rmag ∈ rmags
		j = 1
		θstart = asin(√(1/rmag))
		θend = π/2
		θspan = θend - θstart
		θstep = θspan/mθ
		θend = θend - θstep
		θs = [θstart:θstep:θend;]
		for θ ∈ θs
			r = rmag*sin(θ)^2
			source[i,j] = grad_integrated_over_star(r,θ)/grad_integrated_all_directions(r,θ)
			j+=1
		end 
		i += 1
	end 
	return source
end


end